### 1. Download all the files to your folder
```
Your_folder
├── bind_index
├── cluster.json
├── config.yaml
├── ME_index
├── sample2json.py
├── sciHiCAR_18bp_barcode_440k.txt.gz
├── Snakefile
├── fq  # move your raw fastq files to this folder
│   ├── DNA_example_R1_001.fastq.gz
│   └── DNA_example_R2_001.fastq.gz
└── script
    ├── barcode_hash_v2.py
    ├── extract_barcode_based_on_knee.R
    ├── fq_barcode_correction.py
    ├── raw_fq_update.py
    ├── read_summary.R
```

### 2.Create samples.json file

`python3 sample2json.py --fastq_dir fq`

### 3. Run snakemake pipeline

`snakemake --latency-wait 60 -p -j 99 --cluster-config cluster.json --cluster "sbatch -p common,scavenger -J {cluster.job} --mem={cluster.mem} -N 1 -n {threads} -o {cluster.out} -e {cluster.err} " &> log &`
