% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratFuncs.R
\name{subsetSeurat2cellPhone}
\alias{subsetSeurat2cellPhone}
\title{subsetSeurat2cellPhone}
\usage{
subsetSeurat2cellPhone(seuratObj, annoColumn, no.cells, prefix)
}
\arguments{
\item{seuratObj}{the Seurat object}

\item{annoColumn}{name of metadata column you want}

\item{no.cells}{the number of cells per cluster you want}

\item{prefix}{file prefix for writing out data for cellphoneDB}
}
\value{
invisible(trimmedSeurat)
}
\description{
trims a seurat object to a smaller number of cells (ranked by number of genes expressed) per user defined cluster. writes out text files for cellphoneDB. Makes cellphoneDB a little quicker to run
}
\examples{
ti<-subsetSeurat2cellPhone(seuratObj=integrated,annoColumn="SCT_snn_res.0.15",no.cells=50,prefix="small")

#on cluster
module load briglo/miniconda/3
source activate cellphone
qsub -V -cwd -b y -j y -pe smp 8 -N cpdb_1 "cellphonedb method statistical_analysis small_meta.txt small_counts.txt --project-name small --threshold 10 --threads 8"
EASY!!!

}
