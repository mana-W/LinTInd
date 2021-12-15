# LinTInd

## Installation via GitHub
```
devtools::install_github("mana-W/LinTInd")
```
* Depends:<br />
  + ggplot2<br />
  + parallel<br />
  + stats<br />
  + S4Vectors<br />
  + data.tree<br />
  + reshape2<br />
  + networkD3<br />
  + stringdist<br />
  + purrr<br />
  + ape<br />
  + cowplot<br />
  + ggnewscale<br />
  + stringr<br />
  + dplyr<br />
  + rlist<br />
  + pheatmap<br />
  + Biostrings<br />
  + IRanges<br />
  + BiocGenerics(>= 0.36.1)<br />
  + ggtree<br />

## Installation via Bioconductor
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("LinTInd")
```

## Usage
**Input file：**<br />
*Example files is in LinTInd/inst/extdata*<br />
data is from CB_UMI<br />
fa is ref file<br />
cutsite is a file define each sgRNA start and end positon<br />
celltype.tsv is a file include cell barcode and its' annotations, header: Cell.BC Cell.type

### Quick start
```
library(LinTInd)
data<-paste0(system.file("extdata",package = 'LinTInd'),"/CB_UMI")
fafile<-paste0(system.file("extdata",package = 'LinTInd'),"/V3.fasta")
cutsite<-paste0(system.file("extdata",package = 'LinTInd'),"/V3.cutSites")
celltype<-paste0(system.file("extdata",package = 'LinTInd'),"/celltype.tsv")
data<-read.table(data,sep="\t",header=T)
ref<-ReadFasta(fafile)
cutsite<-read.table(cutsite,col.names = c("indx","start","end"))
scarref<-ReadCutsite(cutsite)
scarref_all<-ReadCutsite(cutsite,reftype="All")
celltype<-read.table(celltype,header=T,stringsAsFactors=F)
```
Or load the example data
```
data("example_data",package = "LinTInd")
```

### Array identify<br />
Alignment

```
scarinfo<-FindIndel(data=data,scarfull=ref,scar=cutsite,indel.coverage="All",type="test",cln=8)
scarform<-IndelForm(scarinfo,cln=4)

```
Define scar pattern for each cell<br />
```
cellsinfo<-IndelIdents(scarinfo,method.use="umi.num",cln=4)
```

Pattern visualization <br />
```
IndelPlot(cellsinfo = cellsinfo)
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/Indel_pattern.png" width = "620" height = "450" align=center />
</p >
<br />


### Indel extracted
```
tag<-TagProcess(cellsinfo$info,Cells=celltype)
```

### Tree reconstruct 
```
treeinfo<-BuildTree(tag)
```

### Visualization

**Similarity of each pair of clusters**
```
tag_dist=TagDist(tag,method = "Jaccard")
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/Indel.png" width = "500" height = "300" align=center />
</p >

<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/cluster_similarity.png" width = "490" height = "450" align=center />
</p >

***Visualization for tree***
```
plotinfo<-PlotTree(treeinfo = treeinfo,data.extract = "T",annotation = "T")
plotinfo$p
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/tree_tag_pattern.png" width = "400" height = "500" align=center />
</p >

Or <br />

```
plotinfo<-PlotTree(treeinfo = treeinfo,data.extract = "T",annotation = "F")
plotinfo$p
```
<p align="center">
<img src="https://github.com/mana-W/scar-barcode/blob/main/image/tree_pattern.png" width = "400" height = "500" align=center />
</p >
