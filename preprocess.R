## Script purpose: Preprocess samples, transform to SingleCellExperiment object and run clustering (example of Lung7)

library(semla)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(BayesSpace)
library(ggplot2)
library(tidyverse)
library(SingleCellExperiment)
library(SummarizedExperiment)


#create directory to store figures and RDS files
#if(!file.exists("../Lung7/figures")) dir.create("../Lung7/figures")
#if(!file.exists("../Lung7/RDS")) dir.create("../Lung7/RDS") 

# Absolute paths to my personal images and rds directories
fig_directory <- ".../GBM_experiment/Images"
rds_directory <- ".../GBM_experiment/RDS"

# Path of each file
h5_matrix <- ".../GBM_experiment/Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5"
tissue_position <- ".../GBM_experiment/spatial/tissue_positions_list.csv"
tissue_hires <- ".../GBM_experiment/spatial/tissue_hires_image.png"
scalefactors_json <- ".../GBM_experiment/spatial/scalefactors_json.json"

# Load the data into a table
infoTable <- tibble(samples=h5_matrix, imgs=tissue_hires, spotfiles=tissue_position, json=scalefactors_json, sample_id=c("GBM")) # Add additional column

# Create Seurat object through STutility, keep genes that have at least 5 counts in the whole tissue
#STobject <- ReadVisiumData(infotable = infoTable, 
#                           platform =  "Visium", 
#                           minUMICountsPerGene = 5)

# This makes a seurat object
seurat_object <- ReadVisiumData(infoTable, 
                                remove_spots_outside_tissue = TRUE, 
                                remove_spots_outside_HE = TRUE,
                                verbose = TRUE)


# This DOES seem to be working
seurat_object <- SubsetSTData(seurat_object, expression = nFeature_Spatial > 5)






# Read in table of types of genes
# Needed to clean up the annotLook.txt because it had a bunch of weird stuff in it
annotLookup <- read.table(".../Code/updated_annotLook.txt", header = T, sep = "\t")


# REMOVE MITOCHONDRIAL, RIBOSOMAL, AND NON-CODING GENES

# Gene names that code for protein
protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]

# List of genes present in sample
present_genes <- rownames(seurat_object)
length(present_genes)

# Keep protein coding genes
# Filter to keep the intersection of protein_genes and selected_genes
selected_genes <- protein_genes[protein_genes %in% present_genes]

# Remove Ribosomal genes
selected_genes <- selected_genes[!grepl("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", selected_genes)]

# Remove Mitochondrial genes
selected_genes <- selected_genes[!grepl("^MT-", selected_genes)]

# Reset seurat_object to only include information relating to the selected_genes of interest
seurat_object <- SubsetSTData(seurat_object, features = selected_genes)


# Load the image in Staffli object 
seurat_object <- LoadImages(seurat_object, verbose = TRUE)
ImagePlot(seurat_object)


# POSSIBLY NOT NECESSARY BECAUSE READVISIUMDATA
# Select spots that are out of the tissue and filter out
#seurat_object <- FeatureViewer(seurat_object)
#seurat_object <- SubsetSTData(seurat_object, features = "need_to_remove")


# Quality Control -> remove spots having a low number of Features
# It doesn't work with nFeature_RNA or nCount_RNA, but does work when I replace RNA with Spatial
# VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA"))

# Mustafa's code seems to have stored the stuff in the "data" layer. Our stuff is stored in the "counts" layer 
# in the Spatial assay

# Save Violin Plot prior to filtering
VlnPlot(seurat_object, features = c("nFeature_Spatial", "nCount_Spatial"), layer = "counts")
#ggsave(".../GBM_experiment/Images/VlnPlot_prefiltering.pdf")

# Filter 
# Save Violin Plot after filtering
seurat_object <- SubsetSTData(seurat_object, expression = nFeature_Spatial > 250)
VlnPlot(seurat_object, features = c("nFeature_Spatial", "nCount_Spatial"), layer = "counts")
#ggsave(".../GBM_experiment/Images/VlnPlot_postfiltering.pdf")





# Load image in Seurat object to use Seurat plotting functions
img <- Seurat::Read10X_Image(image.dir = '.../GBM_experiment/spatial/', image.name = "tissue_lowres_image.png")

Seurat::DefaultAssay(object = img) <- "Spatial"

rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")
img <- img[colnames(x = seurat_object)]

seurat_object[['GBM_experiment']] <- img

# SpatialFeaturePlot(STobject, features = "nFeature_RNA")
MapFeatures(seurat_object, features = "nFeature_Spatial")
#ggsave(".../GBM_experiment/Images/nFeature.pdf")





# Normalization of the data returning all genes (SCT)
seurat_object <- SCTransform(seurat_object,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "Spatial", vst.flavor = "v1")

# SCTransform makes a SCT assay with data layer.
# We have just normalized the data, and it is stored in data layer (which is supposed to hold normalized data)

saveRDS(seurat_object, ".../GBM_experiment/RDS/gbm.rds")




########## Create SingleCellExperiment object with SCT normalization ###########

spatial_dir <- '.../GBM_experiment/spatial/'

colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header = FALSE)

colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")


colData$spot <- paste0(colData$spot, "_1")
rownames(colData) <- colData$spot
colData <- colData[colData$in_tissue > 0, ]

# Create SCE object by transferring SCT assay and metadata (Spatial coordinates)

## coercion from SummarizedExperiment
# se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
# as(se, "SingleCellExperiment")

ST_sce <- SummarizedExperiment(assays = list(SCT = as.matrix(GetAssayData(seurat_object, layer = "data"))))
ST_sce <- as(ST_sce, "SingleCellExperiment")

#ST_sce <- SingleCellExperiment(assays = list(SCT = as.matrix(GetAssayData(seurat_object, layer = "data"))))
#ST_sce <- SingleCellExperiment(assays = list(SCT = as.matrix(GetAssayData(seurat_object, layer = "data"))), colData = colData[colnames(ST_sce@assays@data$SCT),])

# Attempt Dimensionality Reduction
#ST_sce <- SingleCellExperiment::reducedDim(ST_sce, colnames(colData), withDimnames=TRUE):
  

# Error in SummarizedExperiment(...) : 
# the rownames and colnames of the supplied assay(s) must be NULL or identical to 
# those of the SummarizedExperiment object (or derivative) to construct




metadata(ST_sce)$BayesSpace.data <- list()
metadata(ST_sce)$BayesSpace.data$platform <- "Visium"
metadata(ST_sce)$BayesSpace.data$is.enhanced <- FALSE


# Cluster with BayesSpace
set.seed(100)
# Compute top 2000 HVGs and 15 PCA using the SCT data

# PCA = Principle Component Analysis
ST_sce <- ProjectDim(object = ST_sce@assays@data@listData$SCT, reduction = "pca")


ST_sce <- spatialPreprocess(ST_sce, platform = "Visium", n.PCs = 15, n.HVGs = 2000, log.normalize = FALSE, assay.type = "SCT", skip.PCA = TRUE)





# In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) : collapsing to unique 'x' values

# Fin q for clusters
ST_sce <- qTune(ST_sce, qs = seq(2, 20), platform="Visium", d=15)
p <- qPlot(ST_sce)

pdf("../Lung7/figures/BayesQplot.pdf", width = 8, height = 7)
print(p)
dev.off()

#number of clusters
q = 11
#number of PCA
d = 15

# Run the clustering 
ST_sce <- spatialCluster(ST_sce, q=q, d=d, init.method = "kmeans", platform = "Visium", nrep = 50000, gamma = 3)
ST_sce$spatial.cluster <- as.factor(ST_sce$spatial.cluster)

palette <- RColorBrewer::brewer.pal(q, name = "Paired")
clusterPlot(ST_sce, palette = palette, size=0.05) + labs(title = "Spot-level clustering")
saveRDS(ST_sce, "../Lung7/RDS/Lung7_sce.rds")
