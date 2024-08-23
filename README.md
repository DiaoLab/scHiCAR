# scHiCAR
## Pipelines for scHiCAR data processing

### Tools requirement

- **Python**: python3.7 or later
- **snakemake**:  `pip install snakemake==5.13.0`
- **cutadapt**: `pip install cutadapt==3.3`
- **STAR**: [v2.7.5c](https://github.com/alexdobin/STAR/releases/tag/2.7.5c)

### 1. Process FASTQ files of RNA library with Snakemake ([1_RNA](https://github.com/monnneee/scHiCAR/tree/main/1_RNA))
a. Extract ***RNA barcodes*** from the read sequence and add them to the beginning of read 1 (*_R1_001.fastq). Remove the adaptors from the read sequence. If a read sequence does not contain any RNA barcodes, remove the entire read.

b. generate filtered matrix (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) for use in standard scRNA-seq downstream analysis.

### 2. Process FASTQ files of DNA library with Snakemake ([2_DNA](https://github.com/monnneee/scHiCAR/tree/main/2_DNA))
a. Extract ***DNA barcodes*** from the read sequence and add them to the read name. Remove the adaptors from the read sequence. If a read sequence does not contain any DNA barcodes, remove the entire read.

b. Generate a filtered list of DNA barcodes by removing background noise

### 3. Generate pseudo-bulk FASTQ files of DNA library ([3_create_pseudo-bulk_fastq](https://github.com/monnneee/scHiCAR/tree/main/3_create_pseudo-bulk_fastq))
Match the ***DNA barcodes*** corresponding to the cells in each cluster or cell type identified from the RNA library. Extract reads with read names containing these DNA barcodes from the processed FASTQ files of the DNA library, and generate a new FASTQ files for each cluster or cell type.

The downstream scripts for pseudo-bulk FASTQ files used in scHiCAR paper refer to the [nf-core/hicar](https://github.com/jianhong/hicar/tree/dev2rc) pipeline, the details can be found in [paper_scripts.md](https://github.com/monnneee/scHiCAR/blob/main/3_create_pseudo-bulk_fastq/paper_scripts.md)

**Note: the \*.fastq.gz files uploaded in [GSE267126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267126) and [GSE267117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267117) were already processed with Snakemake.**
