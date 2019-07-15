#' makeSeuratList 
#' 
#' reads multiple 10X directories into a list of Seurat Objects
#'
#' @param data.dirs a vector of  paths to cellranger output
#' @param min.cells minimum cells expressing gene to retain (per data.dir), defaults to 3
#' @param min.features minimum genes expressed per cell to retain (per data.dir), defaults to 200
#'
#' @return a list of seurat objects
#'
#' @examples
#' snam<-dir()
#' id <- makeSeuratList(snam)
#'
#' @export
makeSeuratList<-function(data.dirs=NULL,min.cells=3,min.features=200){
return(lapply(data.dirs, function(x) {
    dd<-unlist(strsplit(x,"\\/")
    pnam<-dd[grep('outs',dd)-1]
rd <- Read10X(x) #read in 10x RNA
cd<-read.csv(file = paste0(x,"/CITE/*.csv"), sep = ",",header = TRUE, row.names = 1)
percMouse<-Matrix::colSums(rd[grepl("^mm10-",rownames(rd))])/Matrix::colSums(rd) #replaces seurat build just fot percentage featureset

colMat<-CollapseSpeciesExpressionMatrix(rd,prefix = "GRCh38_",controls = "mm10_", ncontrols =100)
SO<-CreateSeuratObject(counts=colMat,project=gsub('^VENCHI_Sample',"",pnam), min.cells=min.cells, meta.data=percMouse)
SO[["ADT"]] <- CreateAssayObject(counts = cd)
return(SO)
}))
}

#to ask gammy
#R settings
options(scipen=10000)
invisible(utils::memory.limit(64000))
options(future.globals.maxSize= 16384 * 8192^2)
# options(future.globals.maxSize= 4096 * 2048^2)
plan(strategy = "multicore", workers = 2)
scratch_or_preprocessed <- 1  #1 = scratch
No_samples <- 2 #Number of samples if scratch

In_object_storage_method <- 1    #1 = Rdata, all else is RDS
umi_min <- 400
gene_min <- 50
mito_max <- 30
gene_max <- Inf
umi_max <- Inf
human_reference <- "GRCh38_"
mouse_reference <- "mm10_"
mouse_genes_keep <- 100
get function(R)


#' buildMasterSeurat
#'
#' turns a list of Seurat objects into a master object
#'
#' @param seuratList the list of Seurat objects, usually from makeSeuratList
#' @param ndims dimensions to use in integration,SCT and UMAP, default 1:30
#' @param mt.cutoff numeric percentage cutoff for mito expression, defaullt 10
#' @param regress.cell.cycle logical, whether to regress out cell cycle along with mito and nFeature
#' @param res.vec vector of resolutions to find clusters at, defaults to c(seq(0,1,0.1),1.5,2)
#'
#' @return an integrated master seurat object
#'
#' @examples
#' NULL
#'
#' @export
buildMasterSeurat<-function(seuratList,ndims=1:30,mt.cutoff=10,regress.cell.cycle=TRUE,res.vec=c(seq(0,1,0.1),1.5,2)){
    require(Seurat)
    data(cell.cyclegenes)
anchors <- FindIntegrationAnchors(object.list = seuratList, dims = ndims)
integrated <- IntegrateData(anchorset = anchors, dims = ndims) #this takes ages
integrated@active.assay="RNA"
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")
integrated<-subset(integrated, subset = percent.mt < mt.cutoff)
integrated<-CellCycleScoring(object = integrated, s.features = toupper(s.genes), g2m.features = toupper(g2m.genes), set.ident = FALSE)

integrated<- NormalizeData(integrated, assay = "ADT", normalization.method = "CLR")
integrated<- ScaleData(integrated, assay = "ADT")

integrated@active.assay="integrated"
if(regress.cell.cycle) {
    integrated<-SCTransform(integrated, vars.to.regress = c("percent.mt","nFeature_RNA","S.Score","G2M.Score"), verbose = FALSE) } else {
        integrated<-SCTransform(integrated, vars.to.regress = c("percent.mt","nFeature_RNA"), verbose = FALSE)
    }


integrated<-RunPCA(integrated, verbose = FALSE)
integrated<-RunUMAP(integrated, dims = ndims, verbose = T)
integrated<-FindNeighbors(integrated, dims = ndims, verbose = FALSE)
for (res in resvec) integrated<-FindClusters(integrated,resolution=res, verbose = FALSE)
return(integrated)
}




#' splitMasterSeurat
#'
#' turns a Seurat objects into a list of seurat objects
#'
#' @param seuratObj the list of Seurat objects, usually from makeSeuratList
#' @param metaData the metadata column you want to split on
#' @param keep a character vector of a subset of metadata that you want to keep, default NULL
#' @param ndims dimensions to use in integration,SCT and UMAP, default 1:30
#' @param res.vec vector of resolutions to find clusters at, defaults to c(seq(0,1,0.1),1.5,2)
#'
#' @return an list of seurat objects
#'
#' @examples
#' NULL
#'
#' @export
splitMasterSeurat<-function(seuratObj,metaData,keep=NULL,ndims=1:30,resVec=c(seq(0,1,0.1),1.5,2)){
splCells<-split(rownames(seuratObj@meta.data), integrated@meta.data[,metaData])
if(is.null(keep)) cu<-splCells else cu <- splCells[names(splCells) %in% keep]
comb<-lapply(cu, function(X) {
tmp<-subset(integrated,cells=X)
tmp<-RunUMAP(tmp, dims = ndims, verbose = T)
tmp<-FindNeighbors(tmp, dims = ndims, verbose = FALSE)
 for (res in resVec) tmp<-FindClusters(tmp, verbose = FALSE,resolution=res,force.recalc=T)
return(tmp)
})
}