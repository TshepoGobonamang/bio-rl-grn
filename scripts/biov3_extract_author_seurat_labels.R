#!/usr/bin/env Rscript
# Extract only verified per-cell identities from ORIGINAL authors' Seurat RDS.
# The RDS source must be independently obtained, checksum-pinned and referenced
# by the original MAESTER paper. No marker-based relabeling or numeric recoding.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop('Usage: Rscript biov3_extract_author_seurat_labels.R ORIGINAL_Seurat_Keep.rds author_labels.csv')
source <- args[1]; target <- args[2]
if (!file.exists(source) || file.info(source)$size==0) stop('Original author RDS not supplied')
seu <- readRDS(source)
if (!methods::isS4(seu) || !('meta.data' %in% methods::slotNames(seu))) stop('File is not an expected Seurat S4 object')
meta <- methods::slot(seu,'meta.data')
if (!'CellType' %in% names(meta)) stop('Original author Seurat CellType metadata missing')
if (!identical(sort(unique(as.character(meta$CellType))), c('BT142','K562'))) stop('Original author CellType is not categorical BT142/K562')
if (anyDuplicated(rownames(meta)) || anyNA(rownames(meta)) || any(rownames(meta)=='')) stop('Missing or duplicated original author cell IDs')
out <- data.frame(cell_id=rownames(meta),cell_line=as.character(meta$CellType))
write.csv(out,target,row.names=FALSE,quote=TRUE)
cat('Original author Seurat cells:',nrow(out),'\n')
print(table(out$cell_line))
