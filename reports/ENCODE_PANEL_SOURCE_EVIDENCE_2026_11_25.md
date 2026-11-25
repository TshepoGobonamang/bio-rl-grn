# Evidence review: pre-selected K562 CRISPRi panel

Date: 2026-11-25. **This is a public record-level design review, not execution of RNA-seq differential expression.**

| TF | ENCODE experiment | Listed isogenic biological replicates | Shared control link | Original record |
|---|---|---:|---|---|
| TFDP1 | ENCSR552HTC | 2 | ENCSR095PIC | https://www.encodeproject.org/experiments/ENCSR552HTC/ |
| KLF2 | ENCSR620DOP | 2 | ENCSR095PIC | https://www.encodeproject.org/experiments/ENCSR620DOP/ |
| NR2F2 | ENCSR688UTN | 3 | ENCSR095PIC | https://www.encodeproject.org/experiments/ENCSR688UTN/ |
| RELA | ENCSR052ZPX | 2 | ENCSR095PIC | https://www.encodeproject.org/experiments/ENCSR052ZPX/ |
| SP1 | ENCSR650QAK | 3 | Control linkage to verify live | https://www.encodeproject.org/experiments/ENCSR650QAK/ |
| NR2C2 | ENCSR808HQN | 3 | ENCSR095PIC | https://www.encodeproject.org/gene-silencing-series/ENCSR004ZBD/ |
| Non-targeting control | ENCSR095PIC | 2 | Not applicable | https://www.encodeproject.org/experiments/ENCSR095PIC/ |

ENCODE records describe all as K562 isogenic CRISPRi RNA-seq with released metadata. **SP1 control linkage is not visible in the summary retrieved on 2026-11-25** and must pass the live original JSON audit; do not infer it solely from the shared lab. Both CRISPRi TFDP1 and its control records carry portal audit warnings, requiring QC examination.

A live source audit must confirm exact target, control link, distinct replicate biosample IDs, per-replicate gene quantification provenance and matching GRCh38 annotation. ENCODE's bulk RNA standards require separate biological replicates; replicate numbers alone and merged quantification files do not establish independently analyzable effects.

**Interpretation:** perturbing a TF and observing downstream expression change supports a transcriptional response under the assay's conditions. It does not distinguish direct TF-DNA binding from downstream responses, nor does an unresponsive gene become a proven regulatory negative.

Original panel purpose is a small, deliberately disclosed **pilot biological reference**. Six TFs may be too few for a strong across-TF population conclusion, even if all controls are verified. Extend only by a newly frozen, provenance-compatible confirmatory panel; never add only favorable outcomes to this panel.
