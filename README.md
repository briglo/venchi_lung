# venessa's worflow
A workflow for venessa chin

*This should look like a hack version of those fancy seurat tutorials that i told you to look at*

## bits and bobs
* youll need to make a RStudio project somewhere on your laptop
* you'll need to get the r object pandora/volumes/Cancer/brian\ gloss/190531_integrated.rdata and put it in your project path
* we'll need to install a bunch of R libraries, plus maybe some others that ive forgotten
```R
install.packages("BiocManager")
library(BiocManager)
install(c("Seurat","dplyr","biomaRt","ggplot2","clusterProfiler","ReactomePA","igraph","cowplot"))
```

## an R cheatsheet
```R
setwd("NEW/PATH") #move to a new directory
ls() # gives a list of all objects currently loaded (in case you forgot what name it was)
head(seuratObj@meta.data) # shows the metadata columns of the seurat bject you can group by
pdf('test.pdf',width=36,height=16) #makes a pdf  of set dimensions to plot into, good if hard to see...

```

## reading in files
there are two ways really
1) manually `genelist=c("TP53","SPOCK2")`
2) from files (dependent on type) ` read.table("FILENAME", header=T,stringsAsFactors=F) #pretty common, will read a text file saved from an excel doc`

## simple seurat  functions
* the heatmap `DoHeatmap(seuratObj,group.by="orig.ident",features=YOURGENE(S),disp.max=200)`
* The violin plot `VlnPlot(seuratObj,split.by="orig.ident",group.by="orig_final.ident",features=YOURGENE(S),pt.size=.1)`
* the dot plot `DotPlot(seuratObj,features=YOURGENE(S),group.by="orig_final.ident",do.return=T) + theme(axis.text.x = element_text(angle = 90))`
* the UMAP plot `UMAPPlot(seuratObj,group.by='orig_final.ident',label=T)`
* add a meta score `seuratObj<-AddModuleScore(seuratObj,features = YOURGENELIST, name="my_meta_score")`
* plot expression of something in UMAP space `FeaturePlot(seuratObj,reduction="umap",features=YOURGENE(S),split.by=NULL))`

## some of the "bespoke functions"
* CellPhoneDB e.g.
```R
ti<-subsetSeurat2cellPhone(seuratObj=integrated,annoColumn="SCT_snn_res.0.15",no.cells=50,prefix="small")
```
```bash 
     #on cluster
     module load briglo/miniconda/3
     source activate cellphone
     qsub -V -cwd -b y -j y -pe smp 8 -N cpdb_1 "cellphonedb method statistical_analysis small_meta.txt small_counts.txt --project-name small --threshold 10 --threads 8"
     EASY!!!
```
```R
#back in R
cd("PATH/TO/CELLPHONEDB/OUT/small")
cellphoneDB_data<-preprocCellphone(varval=0,pval=.05)
intgraph(cellphoneDB_data$countdat, scoreCut = 0.3, numberCut = 0, numberSplit = 35)
```


* Reactome Analysis e.g.
```
markers<-makeReactomePipe(integrated,"SCT_snn_res.0.8")
for (i in 1:length(markers$CP_result)) dotplot(markers$CP_result[[i]]) + ggtitle(names(markers$CP_result)[i])
```


* Count number of cells expressing above a cutoff  subset by grouping data
```
countTab<-countByCut(seuratObj=integrated,geneName="PFN1",expCut=2,groupBy="orig_final.ident",splitBy="orig.ident")
```


## for later
so first you want to login to the cluster

[this website](https://intranet.gimr.garvan.org.au/display/PG/Wolfpack+SGE+Cheat+Sheet) is where our group put stuff to get on the cluster
## getting started

1) once you are done, login
2) make an interactive HPC job `qrsh -pe smp 4 -l mem_requested=20G`
3) navigate to the results folder you need `cd /share/ScratchGeneral/PATH`
4) setup your R environment `module load briglo/R/3.6.0`
5) Fire up R `R`
6) Set up tools you need
``` library(Seurat)
library(scFuncs) 
library(reticulate)
use_python("/share/ClusterShare/software/contrib/briglo/miniconda3/envs/magic/bin/python")
```
7) load your seurat object ` load("PATH/TO/SEURAT.rdata")`