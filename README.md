# scHiCAR
![alt text](scHiCAR.png)
## Pipelines for scHiCAR data processing

### Tools requirement

- **Python**: python3.7 or later
- **snakemake**:  `pip install snakemake==5.13.0`
- **cutadapt**: `pip install cutadapt==3.3`
- **pairtools**: `pip install pairtools==1.1.0`
- **snaptools**: `pip install snaptools==1.4.8`
- **STAR**: [v2.7.5c](https://github.com/alexdobin/STAR/releases/tag/2.7.5c)
- **BWA**: [v0.7.17](https://github.com/lh3/bwa/releases/tag/v0.7.17)
- **SAMTOOLS**: [v1.12](https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2)

### 1. Preprocess raw FASTQ files of the RNA library with Snakemake and align sequences to the genome ([code](https://github.com/DiaoLab/scHiCAR/tree/dev/1_RNA_preprocess))

#### a. Snakemake procedures:
- Trim specific sequences at the 5′ end of Read 1.  
- Extract RNA barcodes from Read 1 and append them to its 5′ end.  
- Remove the template-switching oligo (TSO) from the 5′ end of Read 2.  
- Remove poly(A) tails and adaptor sequences from the 3′ end of Read 2.  
- Split FASTQ files into two subsets based on primer strategy: oligo-dT vs. random hexamer.  
- Extract and count all barcodes from the dataset.  
- Correct barcodes that have only one mismatch relative to the whitelist.  
- Compress and merge filtered FASTQ files from the two primer strategies.

The resulting files (`03_corrected_fq/*_all_L001_R*_001.fastq.gz`) are ready for alignment using STAR.

#### b. Generate filtered gene expression matrices (`barcodes.tsv`, `features.tsv`, and `matrix.mtx`) with STAR.

### 2. Preprocess raw FASTQ files of the DNA library with Snakemake ([code](https://github.com/DiaoLab/scHiCAR/tree/dev/2_DNA_preprocess))

#### Snakemake procedures:
- Trim specific sequences at the 5′ end of both Read 1 and Read 2.  
- Extract barcodes from the read sequences and append them to the read names following the `@` symbol.  
- Generate a list of all extracted barcodes and count their occurrences.  
- Compare extracted barcodes against a provided whitelist.  
- Correct barcodes that have only one mismatch relative to the whitelist.  
- Compress the filtered FASTQ files.  
- Remove ME (mosaic end) sequences from the reads.  

The resulting files (`05_cutME_fq/*_cutME_L001_R*_001.fastq.gz`) are ready for generating ATAC fragment files and chromatin contact pair files.

### 3. Generate ATAC fragment files with Snakemake (`*.tsv.gz`)([code](https://github.com/DiaoLab/scHiCAR/tree/dev/3_ATAC_fragment))

#### Snakemake procedures:
- Align R1 reads to the reference genome using Snaptools with BWA, and sort BAM files by read name.
- Convert BAM files to fragment-level BED format.
- Extract high-quality cell barcodes based on the knee point of the barcode rank curve.
- Filter fragments using the high-quality barcode list.

The resulting files (`03/filtered/*.filtered.tsv.gz`) can be used in standard scATAC-seq downstream analysis.

### 4. Generate chromatin contact pair files with Snakemake (`*.dedup.pairs.gz`)([code](https://github.com/DiaoLab/scHiCAR/tree/dev/4_chromatin_contact))

#### Snakemake procedures:
- Align reads using BWA-MEM with -SP mode and convert SAM output to BAM format.
- Parse mapped reads, select valid contact pairs, and retain 28-bp cell barcodes.
- Flip contact pairs to generate an upper-triangular matrix.
- Sort contact pairs by cell barcode and genomic position.
- Remove PCR and optical duplicates.
- Summarize read pair statistics.
- Filter contact pairs based on the knee point of the barcode rank curve.


The resulting files (`05_filtered/*.dedup.filtered.pairs.gz`) can be used in downsteam pseudo-bulk or single-cell contact analysis.

### 5. Downstream pseudo-bulk / single-cell analysis ([code](https://github.com/DiaoLab/scHiCAR/tree/dev/5_downstream_analysis))

### Citation
Wei, X., Xu, Y., Yang, D. et al. Trimodal single-cell profiling of transcriptome, epigenome and 3D genome in complex tissues with scHiCAR. *Nat Biotechnol* (2026). https://doi.org/10.1038/s41587-026-03013-7

### Contributors
This pipeline was contributed by [Yueyuan Xu](https://github.com/yueyuanxu), [Xiaolin Wei](https://github.com/wxl123wxl), and [Yu Xiang](https://github.com/yxian9?tab=repositories).
