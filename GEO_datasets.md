# scHiCAR Processed Data

The processed data deposited under [GSE305889](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305889) for the scHiCAR datasets used in our [Nature Biotechnology paper](https://www.nature.com/articles/s41587-026-03013-7) are described below.

## 1. cellline_scHiCAR

Subseries [GSE267117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267117) corresponds to human cell line dataset presented in **Fig. 1** of paper.

The following processed files are available:

* **`.pairsam` file**: Read 2 represents the ATAC-seq reads, while paired Read 1 and Read 2 represent the two ends of each chromatin contact.
* **RNA matrix**: barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
* **Metadata files**: cell-level information for each cell barcode.

> **Note:** For the cell line datasets, the RNA and DNA data were generated from **two different batches**. The two libraries were evaluated separately to assess the quality and advantages of the RNA and DNA libraries of scHiCAR, and were not matched at the individual cell level.

## 2. mouse_brain_scHiCAR_2

Subseries [GSE267126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267126) corresponds to mouse brain dataset (1.62 million cells) presented in **Figs. 3-5** of paper.

For each cell type, the following processed files are available:

* **`.pairsam` file**
* **RNA matrix**
* **Metadata files**

### Extracting ATAC Fragments and Chromatin Contact Pairs

For both the cell line and mouse brain datasets described above, ATAC fragments and chromatin contact pairs can be extracted from the `.pairsam` files by following **Steps 3-6** in our processing pipeline:

[scHiCAR processing pipeline – Step 3: Extract Read 2 from PETs to call open chromatin peaks](https://github.com/DiaoLab/scHiCAR/tree/main/3_create_pseudo-bulk_fastq/paper_scripts.md#3-extract-read2-from-pets-to-call-open-chromatin-peaks-with-q-value-cutoff-001)

## 3. mouse_brain_scHiCAR_1

Subseries [GSE305439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305439) was generated using our **latest protocol** and corresponds to mouse brain dataset (5,313 cells) presented in **Fig. 2** of paper.

Unlike the earlier protocol, the latest protocol captures **ATAC-seq information in Read 1**.

The following processed files are available:

* **ATAC fragment files**
* **Chromatin contact pair files**
* **RNA matrix files**
* **Metadata files**

## 4. mouse_skeletal_muscle_scHiCAR

Subseries [GSE304674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304674) was also generated using our **latest protocol** and corresponds to mouse skeletal muscle dataset presented in **Fig. 6** of paper.

The following processed files are available:

* **ATAC fragment files**
* **Chromatin contact pair files**
* **RNA matrix files**
* **Metadata files**

## Summary of scHiCAR Datasets

| Dataset                       | GEO Accession                                                             | Figure    | Protocol         | Processed Files                                     |
| ----------------------------- | ------------------------------------------------------------------------- | --------- | ---------------- | --------------------------------------------------- |
| Cell line scHiCAR                         | [GSE267117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267117) | Fig. 1    | Earlier protocol | `.pairsam`, RNA matrix, metadata                    |
| Mouse brain scHiCAR (1.62 million cells)  | [GSE267126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267126) | Figs. 3-5 | Earlier protocol | `.pairsam`, RNA matrix, metadata                    |
| Mouse brain scHiCAR (5,313 cells)         | [GSE305439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305439) | Fig. 2    | Latest protocol  | ATAC fragments, contact pairs, RNA matrix, metadata |
| Mouse skeletal muscle scHiCAR             | [GSE304674](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304674) | Fig. 6    | Latest protocol  | ATAC fragments, contact pairs, RNA matrix, metadata |

