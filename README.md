# venessa's worflow
A workflow for venessa chin

** This should look like a hack version of those fancy seurat tutorials that i told you to look at 

so first you want to login to the cluster

[this website](https://intranet.gimr.garvan.org.au/display/PG/Wolfpack+SGE+Cheat+Sheet) is where our group put stuff to get on the cluster
## getting started

1) once you are done, login
2) make an interactive HPC job `qrsh -pe smp 4 -l mem_requested=20G`
3) navigate to the results folder you need `cd /share/ScratchGeneral/PATH`
4) setup your R environment `module load briglo/R/3.6.0`
5) Fire up R `R`
6) Set up tools you need
`` library(Seurat)
library(scFuncs) 
library(reticulate)
use_python("/share/ClusterShare/software/contrib/briglo/miniconda3/envs/magic/bin/python")
``
7) load your seurat object ` load("PATH/TO/SEURAT.rdata")`

## an R cheatsheet
```
setwd("NEW/PATH") #move to a new directory
ls() # gives a list of all objects currently loaded (in case you forgot what name it was)
head(seuratObj@meta.data) # shows the metadata columns of the seurat bject you can group by
pdf('test.pdf',width=36,height=16) #makes a pdf  of set dimensions to plot into, good if hard to see...

```

## reading in files
there are two ways really
1) manually `genelist=c("TP53","SPOCK2")`
2) from files (dependent on type) ` read.table("FILENAME", header=T,stringsAsFactors=F) #pretty common, will read a text file saved from an excel doc

## simple seurat  functions
* the heatmap `DoHeatmap(seuratObj,group.by="orig.ident",features=YOURGENELIST,disp.max=200)
* The violin plot `VlnPlot(seuratObj,split.by="orig.ident",group.by="orig_final.ident",features=YOURGENELIST,pt.size=.1)
* the dot plot `DotPlot(seuratObj,features=YOURGENELIST,group.by="orig_final.ident",do.return=T) + theme(axis.text.x = element_text(angle = 90))`
* the UMAP plot `UMAPPlot(seuratObj,group.by='orig_final.ident',label=T)`


