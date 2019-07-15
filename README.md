# venessa's worflow
A workflow for venessa chin

*This should look like a hack version of those fancy seurat tutorials that i told you to look at*

## getting onto the cluster and firing up R
```bash
cluster # mounts cluster
d1 # connect to dice
brain # get an interactive compute session
scr # currently goes to /share/ScratchGeneral/briglo/scRNA
cd venchi
R
```

## an R cheatsheet
```R
library(Seurat) #loads the library Seurat, required for basically everything
library(venchi) #loads the functions i wrote for you
setwd("NEW/PATH") #move to a new directory
ls() # gives a list of all objects currently loaded (in case you forgot what name it was)
load("PATH/TO/SEURAT.rdata") #loads the seurat object, takes a while
head(integrated@meta.data) # shows the metadata columns of the seurat bject you can group by
#this will be useful when you're playing in the cluster
pdf('test.pdf',width=36,height=16) #makes a pdf  of set dimensions to plot into, good if hard to see...
dev.off() #stops whatever graphics device R is plotting to, necessary to view pdfs
```

## reading in files
there are two ways really
1) manually `genelist=c("GZMA","GZMB","GZMH","GZMK","GZMM","PFN1","FASLG","TNFSF10")`
2) from files (dependent on type) ` read.table("FILENAME", header=T,stringsAsFactors=F) #pretty common, will read a text file saved from an excel doc`

## simple seurat  functions
* the heatmap `DoHeatmap(integrated,group.by="orig_final.ident",features=genelist,disp.max=200)`
* The violin plot `VlnPlot(integrated,split.by="orig.ident",group.by="orig_final.ident",features=genelist,pt.size=.1)`
* the dot plot `DotPlot(integrated,features=genelist,group.by="orig_final.ident",do.return=T) + theme(axis.text.x = element_text(angle = 90))`
* the UMAP plot `UMAPPlot(integrated,group.by='orig_final.ident',label=T)`
* add a meta score `integrated<-AddModuleScore(integrated,features = genelist, name="my_meta_score")`
* plot expression of something in UMAP space `FeaturePlot(integrated,reduction="umap",features=genelist,split.by=NULL))` or `FeaturePlot(integrated,reduction="umap",features='my_meta_score',split_by="orig_final.ident"`

## some of the "bespoke functions"

* The kinds of comparisons we were talking about counting
```R
#for a simple single gene, single cutoff in all cells
countTab<-countByCut(integrated=integrated,geneName="PFN1",expCut=2,groupBy="orig_final.ident",splitBy="orig.ident")

#for a more complicated mulltigene, multi cutoffs
tab<-data.frame(gene=c("PFN1","GZMA","GZMB"),cuts=c(2,0,0))
#adds new meta.data columns containing info of interest, can be used like below in Reactome Analysis
integrated<-annotateByCuts(seuratObj=integrated,geneTab=tab,groupBy="orig_final.ident",groupID="CD8",splitBy="orig.ident")

outTab #a table of counts for cells matching groupID FYE
```


* Reactome Analysis e.g.
```R
#for all groupings made by a particular metadata column
markers<-makeReactomePipe(integrated,"SCT_snn_res.0.8")
for (i in 1:length(markers$CP_result)) dotplot(markers$CP_result[[i]]) + ggtitle(names(markers$CP_result)[i])

#for combining two metadata columns and directly comparing two groups(split by something else if you want)
#first extract the relevant markers and universe
enrichment<-makeSpecificMarkers(integrated,groupBy_1="orig_final.ident",groupID_1="CD8 T-cell 1",groupBy_2="orig.ident",groupID_2=c("ACITE","BCITE"),splitBy=NULL,splitID=NULL)
out<-makeReactomeForMarkers(enrichment$markerLists,enrichment$hasEntrez)
dotoplot(out)

```



## bits and bobs
* terminal is your friend 
* we need to make sure you are running a "current" version of R `R --version`
* youll need to make a RStudio project somewhere on your laptop
* you'll need to get the r object pandora/volumes/Cancer/brian\ gloss/190531_integrated.rdata and put it in your project path
* we'll need to install a bunch of R libraries, plus maybe some others that ive forgotten
```R
install.packages("BiocManager")
library(BiocManager)
install(c("devtools","Seurat","dplyr","biomaRt","ggplot2","clusterProfiler","ReactomePA","igraph","cowplot"))
devtools::install_github("https://github.com/briglo/venchi_lung.git")
```


* CellPhoneDB e.g.
```R
ti<-subsetSeurat2cellPhone(seuratObj=integrated,annoColumn="SCT_snn_res.0.15",no.cells=50,prefix="small")
```
```bash 
#running cellphoneDB
module load briglo/miniconda/3
source activate cellphone
qsub -V -cwd -b y -j y -pe smp 8 -N cpdb_1 "cellphonedb method statistical_analysis small_meta.txt small_counts.txt --project-name small --threshold 10 --threads 8"
#EASY!!!
```
```R
#back in R
cd("PATH/TO/CELLPHONEDB/OUT/small")
cellphoneDB_data<-preprocCellphone(varval=0,pval=.05)
intgraph(cellphoneDB_data$countdat, scoreCut = 0.3, numberCut = 0, numberSplit = 35)
```




# setting up: for later
## getting onto the cluster
so first you want to login to the cluster this was originally done by doing [this really](https://intranet.gimr.garvan.org.au/display/PG/A+guide+for+getting+a+mac+onto+the+cluster)

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