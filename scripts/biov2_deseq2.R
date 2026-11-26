#!/usr/bin/env Rscript
# BIO-V2: DESeq2 treatment-control comparisons of original ENCODE RSEM gene-level
# estimated counts using exact matched RSEM TPM and effective-length offsets.
# The design has small isogenic biological replicate numbers (2-3 per arm).
# Statistical response != direct TF binding, absence of DE != biological negative.
suppressPackageStartupMessages({library(DESeq2); library(tximport)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('Usage: Rscript biov2_deseq2.R <qc_root> <output_root>')
qc_root <- normalizePath(args[1]); out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
all_summary <- list()
for (tf in c('TFDP1','KLF2','NR2F2','RELA','SP1','NR2C2')) {
  location <- file.path(qc_root, tf)
  read_mat <- function(name) {
    # read.csv() manages the lifecycle of an unopened gzfile connection.
    # Closing it a second time caused "invalid connection" before DESeq2.
    path <- file.path(location, name)
    if (!file.exists(path) || file.info(path)$size <= 0) {
      stop(tf, ': missing or empty original QC matrix: ', name)
    }
    read.csv(gzfile(path, open='rt'), header=TRUE,
             check.names=FALSE, row.names=1)
  }
  counts <- as.matrix(read_mat('rsem_expected_counts.csv.gz'))
  abundance <- as.matrix(read_mat('rsem_TPM.csv.gz'))
  length_mat <- as.matrix(read_mat('rsem_effective_length.csv.gz'))
  meta <- read.csv(file.path(location,'coldata.csv'),stringsAsFactors=FALSE)
  if (!identical(rownames(counts),rownames(abundance)) ||
      !identical(rownames(counts),rownames(length_mat)) ||
      !identical(colnames(counts),meta$sample_id) ||
      !identical(colnames(abundance),meta$sample_id) ||
      !identical(colnames(length_mat),meta$sample_id)) stop(tf,': mismatched gene/sample identities')
  if (length(unique(meta$original_biosample_id))!=nrow(meta) ||
      sum(meta$condition=='control') < 2 || sum(meta$condition=='perturbed') < 2)
      stop(tf,': replicate independence / contrast invalid')
  if (any(!is.finite(counts)) || any(counts<0) ||
      any(!is.finite(abundance)) || any(abundance<0) ||
      any(!is.finite(length_mat)) || any(length_mat<0)) stop(tf,': invalid expression matrices')
  # Predeclared low expression and RSEM gene-length filter. Never impute a zero length.
  # Threshold only uses original estimated counts; no gold or source rankings.
  keep <- rowSums(counts)>=10 & apply(length_mat,1,min)>0
  n_before <- nrow(counts)
  if (sum(keep)<1000) stop(tf,': insufficient genes after fixed count/length filtering')
  counts <- counts[keep,,drop=FALSE]
  abundance <- abundance[keep,,drop=FALSE]
  length_mat <- length_mat[keep,,drop=FALSE]
  coldata <- data.frame(condition=factor(meta$condition,levels=c('control','perturbed')),
                        row.names=meta$sample_id)
  # tximport expects RSEM-style abundance, effective lengths and estimated counts.
  # DESeq2 uses gene-length offsets from tximport, not TPM as counts.
  txi <- list(counts=counts, abundance=abundance, length=length_mat,
              countsFromAbundance='no')
  dds <- DESeqDataSetFromTximport(txi,colData=coldata,design=~condition)
  dds <- DESeq(dds,fitType='parametric',minReplicatesForReplace=Inf,quiet=TRUE)
  res <- results(dds,contrast=c('condition','perturbed','control'),alpha=0.05,
                 independentFiltering=TRUE)
  result_df <- as.data.frame(res)
  result_df$gene_id <- rownames(result_df)
  result_df <- result_df[,c('gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
  result_df$bh_fdr_005 <- !is.na(result_df$padj) & result_df$padj < 0.05
  tfdir <- file.path(out,tf)
  dir.create(tfdir,recursive=TRUE,showWarnings=FALSE)
  write.csv(result_df,gzfile(file.path(tfdir,'deseq2_genewise.csv.gz')),row.names=FALSE,na='')
  summary <- data.frame(tf=tf, n_controls=sum(meta$condition=='control'),
      n_perturbed=sum(meta$condition=='perturbed'), genes_input=n_before,
      genes_retained=sum(keep), tested_pvalues=sum(!is.na(result_df$pvalue)),
      tested_padj=sum(!is.na(result_df$padj)),
      fdr_005=sum(result_df$bh_fdr_005),
      fdr_005_up=sum(result_df$bh_fdr_005 & result_df$log2FoldChange>0,na.rm=TRUE),
      fdr_005_down=sum(result_df$bh_fdr_005 & result_df$log2FoldChange<0,na.rm=TRUE),
      comparison='perturbed_vs_shared_control',
      method='DESeq2_rounded_RSEM_estimated_count_with_tximport_length_offsets',
      is_direct_TF_binding_truth=FALSE,
      is_independent_validation_of_GCEV=FALSE,
      stringsAsFactors=FALSE)
  write.csv(summary,file.path(tfdir,'summary.csv'),row.names=FALSE)
  all_summary[[tf]] <- summary
  print(summary)
}
write.csv(do.call(rbind,all_summary),file.path(out,'deseq2_summary_all_six.csv'),row.names=FALSE)
writeLines(paste0('R_version=',R.version.string,'\nDESeq2_version=',
                  as.character(packageVersion('DESeq2')), '\ntximport_version=',
                  as.character(packageVersion('tximport'))),
           file.path(out,'deseq2_environment.txt'))
