# 1. Match each RNA barcode with its corresponding DNA barcode

This step is based on the RNA object generated from the filtered expression matrix using [Seurat](https://satijalab.org/seurat/articles/install_v5) R package, along with a 6-bp barcode list created during library preparation to match RNA and DNA.

#### 1.1 Load RNA expression matrix
```r
# R
library(Seurat)
data<- Read10X(data.dir = "1_RNA_preprocess/GeneFull/filtered")
rna<-CreateSeuratObject(counts = data, project = "example", min.cells = 3, min.features = 200)
```

#### 1.2 Extract RNA barcode and split into 3 DNA parts
```r
# R
df <- data.frame(RNAlibrary = rna$orig.ident)
df$DNAbarcode1<-substr(rownames(df),1,6)
df$DNAbarcode2<-substr(rownames(df),7,12)
df$DNAbarcode3<-substr(rownames(df),13,18)
```

#### 1.3 Replace the first 6bp using predefined map
```r
# R
barcode_map <- c( #predefined 6-bp RNA–DNA barcode pairing scheme used during library preparation
  "TCATCC" = "TACCCG", #left represents RNA and right represents DNA
  "AGTCAA" = "GAGTTT",
        ...
  "ATGTCA" = "CGATTT",
  "CCCTAC" = "CTGGAG"
)
df$DNAbarcode1 <- sapply(df$DNAbarcode1, function(barcode) {
  for (old_prefix in names(barcode_map)) {
    if (startsWith(barcode, old_prefix)) {
      return(sub(paste0("^", old_prefix), barcode_map[old_prefix], barcode))
    }
  }
  return(barcode)
})
```

#### 1.4 Reverse complement the last 6bp
```r
# R
reverse_complement <- function(seq) {
  complement <- setNames(c("T", "G", "C", "A"), c("A", "C", "G", "T"))
  comp_seq <- unname(complement[unlist(strsplit(seq, ""))])
  rev_seq <- paste(rev(comp_seq), collapse = "")
  return(rev_seq)
}
df$DNAbarcode3<-unname(sapply(df$DNAbarcode3, reverse_complement))
df$DNAbarcode<-paste0(df$DNAbarcode1,df$DNAbarcode2,df$DNAbarcode3) # combine full DNA barcode
df$DNAbarcode1<-NULL
df$DNAbarcode2<-NULL
df$DNAbarcode3<-NULL
```

#### 1.5 Filter cells based on the high-quliaty barcode list generate from [3_ATAC_fragment](https://github.com/DiaoLab/scHiCAR/tree/dev/3_ATAC_fragment) and [4_chromatin_contact](https://github.com/DiaoLab/scHiCAR/tree/dev/4_chromatin_contact)
```r
# R
atac_cut_rank<-read.table("3_ATAC_fragment/02_fragment/*.barcode.cut_rank")$V1
pairs_cut_rank<-read.table("4_chromatin_contact/03_dedup/*.barcode.cut_rank")$V1
valid_DNA_barcodes <- intersect(atac_cut_rank, pairs_cut_rank)
df2 <- df[df$DNAbarcode %in% valid_DNA_barcodes, ]
# Note: If you also want to filter cells based on TSS enrichment, please run the ArchR functions (`createArrowFiles` and  `ArchRProject` steps in the "3. Single-cell processing") first and export the processed barcode list for additional filtering here.
```

#### 1.6 Subset RNA object and add DNA barcode information to meta.data
```r
# R
rna_filter <- subset(rna, cells = rownames(df2))
rna_filter <- AddMetaData(rna_filter, metadata = setNames(df2$DNAbarcode, rownames(df2)), col.name = "DNAbarcode")
```

#### 1.7 Cell clustering and annotation
```r
# R
rna_filter <- NormalizeData(rna_filter)
rna_filter <- FindVariableFeatures(rna_filter, selection.method = "vst", nfeatures = 2000)
rna_filter <- ScaleData(rna_filter)
rna_filter <- RunPCA(rna_filter,npcs=30)
rna_filter <- FindNeighbors(rna_filter, reduction = "pca", dims = 1:30)
rna_filter <- FindClusters(rna_filter, resolution = 0.5)
rna_filter <- RunUMAP(rna_filter, reduction = "pca", dims = 1:30)
marker_genes <- FindAllMarkers(rna_filter, only.pos = TRUE)
new.cluster.ids <- c("celltype1", "celltype2", "celltype3",...,"celltypeN")
names(new.cluster.ids) <- levels(rna_filter)
rna_filter <- RenameIdents(rna_filter, new.cluster.ids)
rna_filter$celltype <- Idents(rna_filter) # Add cell type annotation to metadata
DimPlot(rna_filter, reduction = "umap", group.by = "celltype", label = TRUE)
write.table(rna_filter@meta.data,"single_cell_metadata.txt",quote=F,sep='\t',col.names=T,row.names=T) # For downstream single cell ATAC analysis
write.table(rna_filter@meta.data$DNAbarcode,"all_valid_DNA_barcodes.txt",quote=F,col.names=F,row.names=F) # For downstream single cell 3D genome analysis
```

#### 1.8 Extract DNA barcodes for each cell type for downstream pseudo-bulk analyses of the DNA library
```r
# R
meta <- rna_filter@meta.data
barcode_list <- split(meta$DNAbarcode, meta$celltype)
for (celltype in names(barcode_list)) {
  write.table(barcode_list[[celltype]],file=paste0(celltype, ".DNAbarcode.txt"),quote=F,col.names=F,row.names=F)
}
```

# 2. Pseudo-bulk processing

#### 2.1 Generate pseudo-bulk ATAC fragment files and contact pair files for each cell type
```bash
# Bash
for i in {celltype1,celltype2,celltype3,...,celltypeN}
do
python3 extract_ATAC_fragment.py -l 3_ATAC_fragment/03_filtered/*.filtered.tsv -s ${i}.DNAbarcode.txt -o ${i}.ATAC.fragment.tsv
python3 extract_pairs.py -l 4_chromatin_contact/05_filtered/*.dedup.filtered.pairs -s ${i}.DNAbarcode.txt -o ${i}.contact.pairs
done
```

#### 2.2 Call open chromatin peaks for each cell types
macs2: [install](https://github.com/macs3-project/MACS/wiki/Install-macs2)
```bash
# Bash
genome_size=mm # or hs
for i in {celltype1,celltype2,celltype3,...,celltypeN}
do
macs2 callpeak --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits --qval 0.01 --gsize $genome_size --format BED --name ${i} --treatment ${i}.ATAC.fragment.tsv
done
```

#### 2.3 Generate pvalue BIGWIG file for open chromatin visualization
slopBed: [install](https://github.com/arq5x/bedtools2/releases/tag/v2.31.0)

bedClip and bedGraphToBigWig: [install](https://github.com/ENCODE-DCC/kentUtils)
```bash
# Bash
chrom_size=The/PATH/of/chrom/sizes/file
for i in {celltype1,celltype2,celltype3,...,celltypeN}
do
sval=$(wc -l ${i}.ATAC.fragment.tsv | awk '{printf "%f", $1/1000000}') #counts the number of tags per million in the (compressed) BED file
macs2 bdgcmp -t ${i}_treat_pileup.bdg -c ${i}_control_lambda.bdg --o-prefix ${i} -m ppois -S $sval
slopBed -i ${i}_ppois.bdg -g $chrom_size -b 0 | bedClip stdin $chrom_size ${i}.pval.signal.bedgraph
grep -E 'chrX|chrY' ${i}.pval.signal.bedgraph > temp
grep -E -v 'chrX|chrY' ${i}.pval.signal.bedgraph|sort -k1.4n -k 2,2n|cat - temp > ${i}.pval.signal.srt.bedgraph
bedGraphToBigWig ${i}.pval.signal.srt.bedgraph $chrom_size ${i}_sig.pval.signal.bigwig
rm -f ${i}.pval.signal.bedgraph ${i}.pval.signal.srt.bedgraph
rm temp
done
```

#### 2.4 Aggregate read pairs into contact matrix in the cooler format (5kb resolution)
cooler: [install](https://cooler.readthedocs.io/en/latest/quickstart.html)
```bash
# Bash
chrom_size=The/PATH/of/chrom/sizes/file
for j in {celltype1,celltype2,celltype3,...,celltypeN}
do
cat 4_chromatin_contact/05_filtered/*.dedup.pairs.head ${j}.contact.pairs|bgzip -c > ${j}.contact.pairs.gz
pairtools sort -o ${j}.sort.pairs.gz ${j}.contact.pairs.gz
pairix -f ${j}.sort.pairs.gz
cooler cload pairix --max-split 2 --nproc 12 ${chrom_size}:5000 ${j}.sort.pairs.gz ${j}.5000.cool
cooler zoomify --balance -r 5000N -n 12 -o ${j}.5000.mcool ${j}.5000.cool
done
```
The *.cool file can be used to call A/B compartment, TAD, and chromatin loops for each cell types. The *.mcool file can be used to visualization with Higlass.

# 3. Single-cell processing

#### 3.1 single cell ATAC analysis
The `03_filtered/*.filtered.tsv.gz` files can be used  for TSS enrichment score and gene score calculation, as well as cell clustering, using the [ArchR](https://www.archrproject.com/index.html) R package.
```r
# R
library(ArchR)
addArchRThreads(threads = 1)
addArchRGenome("mm10")
createArrowFiles(inputFiles = "03_filtered/*.filtered.tsv.gz", sampleNames = "example", minTSS = 1, minFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE) # Set a small cutoff for minTSS and minFrags if you just want to retain the same cells as in rna_filter@meta.data.
dna <- ArchRProject(ArrowFiles = "example.arrow",outputDirectory = "obj.name", copyArrows = TRUE)
rna_metadata<-read.table("single_cell_metadata.txt",header=T,sep='\t',row.names=1)
dna_filter<-dna[paste("obj.name#",rna_metadata$DNAbarcode,sep=""), ] #filter cells based on the cell barcodes of rna_filter@meta.data
GSmatrix <- getMatrixFromProject(dna_filter,"GeneScoreMatrix")
rownames(GSmatrix)<-rowData(GSmatrix)$name
genescore<-assays(GSmatrix)$GeneScoreMatrix
colnames(genescore)<-gsub("#","_",colnames(genescore))
saveRDS(genescore,"genescore_matrix.rds")
dna_filter <- addIterativeLSI(ArchRProj = dna_filter,useMatrix = "TileMatrix",name = "IterativeLSI",iterations = 4,clusterParams = list(resolution = c(0.2,0.6,1.2),sampleCells = 20000, n.start = 10), varFeatures = 100000,dimsToUse = 1:50,outlierQuantiles = NULL)
dna_filter <- addClusters(input = dna_filter,reducedDims = "IterativeLSI",method = "Seurat",name = "Clusters",resolution = 0.8)
dna_filter <- addUMAP(ArchRProj = dna_filter, reducedDims = "IterativeLSI", name = "UMAP", minDist = 0.2)
plotEmbedding(ArchRProj = dna_filter, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
saveRDS(dna_filter,"dna_filter.rds")
```

#### 3.2 single cell 3D genome analysis
Generate a list of contact pair files, each representing a single cell in tab-delimited format:
```bash
# Bash
python3 extract_pairs.py -l 4_chromatin_contact/05_filtered/*.dedup.filtered.pairs -s all_valid_DNA_barcodes.txt -o all_valid_DNA_contact.pairs
cut -d':' -f1,3 all_valid_DNA_contact.pairs | perl -p -e "s/://g" | cut -f1-5 > all_valid_DNA_contact.pairs.5columns
python3 split_per_cell_pairs.py --input_file all_valid_DNA_contact.pairs.5columns --output_dir per_cell/
```
The `per_cell` folder can be used for cell embedding and imputing single-cell chromatin contact maps with [Fast-Higashi and Higashi](https://github.com/ma-compbio/Higashi).
```python
# Python
### Process the input data：
from higashi.Higashi_wrapper import *
config = ./higashi.JSON"
higashi_model = Higashi(config)
higashi_model.process_data()

### Initialize Fast-Higashi model and turn sparse matrices into sparse tensors:
from fasthigashi.FastHigashi_Wrapper import *
from umap import UMAP
import seaborn as sns
import matplotlib.pyplot as plt
config = "./Fasthigashi.JSON"
fastHigashi_model = FastHigashi(config_path=config,
                     path2input_cache="./Fasthigashi",
                     path2result_dir="./Fasthigashi",
                     off_diag=100,
                     filter=False,
                     do_conv=False,
                     do_rwr=False,
                     do_col=False,
                     no_col=False)
fastHigashi_model.prep_dataset(batch_norm=True)

### Run the model:
fastHigashi_model.run_model(dim1=.6,
                rank=256,
                n_iter_parafac=1,
                extra="")

### Visualize the embedding results:
clusters = fastHigashi_model.label_info['category']
pal1 = {"GLUT":"#ff595e", "GABA":"#1982c4", "NON":"#8ac926"}
embed = fastHigashi_model.fetch_cell_embedding(final_dim=256,restore_order=False)
embedding = embed['embed_l2_norm_correct_coverage_fh']
vec = UMAP(n_components=2, n_neighbors=25, random_state=0).fit_transform(embedding)
fig = plt.figure(figsize=(5, 3))
ax = plt.subplot(1, 1, 1)
sns.scatterplot(x=vec[:, 0], y=vec[:, 1], hue=clusters, ax=ax, s=5, alpha=0.8, linewidth=0, palette=pal1)
handles, labels = ax.get_legend_handles_labels()
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
plt.tight_layout()
plt.savefig("example_fastHigashi_UMAP.pdf", format="pdf")

### train higashi for contact map imputation：
from higashi.Higashi_wrapper import *
config = "./higashi.JSON"
higashi_model = Higashi(config)
higashi_model.prep_model()
higashi_model.train_for_imputation_nbr_0()
higashi_model.impute_no_nbr()
higashi_model.train_for_imputation_with_nbr()
higashi_model.impute_with_nbr()
with open('higashi_model.impute.pkl', 'wb') as f:
    pickle.dump(higashi_model, f)
```
