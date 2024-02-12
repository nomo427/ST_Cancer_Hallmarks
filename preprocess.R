## Script purpose: Preprocess samples, transform to SingleCellExperiment object and run clustering

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


# ======================= CREATE SEURAT_OBJECT ============================ #


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

# Create Seurat Object through Semla
seurat_object <- ReadVisiumData(infoTable, 
                                remove_spots_outside_tissue = TRUE, 
                                remove_spots_outside_HE = TRUE,
                                verbose = TRUE)

# Keep genes that have at least 5 counts in the whole tissue
seurat_object <- SubsetSTData(seurat_object, expression = nFeature_Spatial > 5)


# ============== FILTER / REMOVE MITOCHONDRIAL, RIBOSOMAL, AND NON-CODING GENES ================ #


# Read in table of types of genes
annotLookup <- read.table(".../Code/updated_annotLook.txt", header = T, sep = "\t")

# Gene names that code for protein
protein_genes <- annotLookup$gene_name[annotLookup$gene_type %in% c("protein_coding", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene", "IG_LV_gene", "IG_V_gene", "IG_J_gene", "IG_C_gene" , "IG_D_gene")]

# List of genes present in sample
present_genes <- rownames(seurat_object)

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


# ============== QUALITY CONTROL (remove spots with low number of Features) ================= #


# For reference: 
# nFeature_spatial is the number of genes detected within a cell
# nCount_spatial is the total number of molecules detected within a cell

# Save Violin Plot prior to filtering
VlnPlot(seurat_object, features = c("nFeature_Spatial", "nCount_Spatial"), layer = "counts")
ggsave(".../GBM_experiment/Images/VlnPlot_prefiltering.pdf")

# Filter 
# Save Violin Plot after filtering
seurat_object <- SubsetSTData(seurat_object, expression = nFeature_Spatial > 250)
VlnPlot(seurat_object, features = c("nFeature_Spatial", "nCount_Spatial"), layer = "counts")
ggsave(".../GBM_experiment/Images/VlnPlot_postfiltering.pdf")


# ================== LOAD IMAGE INTO SEURAT OBJECT ================== #


# Load image into seurat_object to use Seurat plotting functions

# Make img (a Visium object)
img <- Seurat::Read10X_Image(image.dir = '.../GBM_experiment/spatial/', image.name = "tissue_lowres_image.png")

# Set img's default assay to "Spatial"
Seurat::DefaultAssay(object = img) <- "Spatial"

# Adds a "_1" to the end of each rowname
rownames(img@coordinates) <- paste0(rownames(img@coordinates), "_1")

# img: the rows are the tags with _1's and the cols are "tissue", "row", "col", "imagerow", "imagecol"
# seurat_object: the rows are genes, and the cols are tags

# Adds a "_1" to the end of each colname in the seurat_object
colnames(seurat_object) <- paste0(colnames(seurat_object), "_1")

# Keeps the items that are present in seurat_object's columns
img <- img[colnames(x = seurat_object)]

# Load the img into seurat_object's images directory as GBM_experiment
seurat_object[['GBM_experiment']] <- img

# Make a plot displaying how many genes were detected in each spot
# SpatialFeaturePlot(STobject, features = "nFeature_RNA")
MapFeatures(seurat_object, features = "nFeature_Spatial")
ggsave(".../GBM_experiment/Images/nFeature.pdf")


# ============ NORMALIZE DATA AND CREATE SINGLECELLEXPERIMENT OBJECT ============ #


# Normalize the data using SCTransform. This makes an SCT assay with data layer.
seurat_object <- SCTransform(seurat_object,return.only.var.genes = FALSE, variable.features.n = NULL, variable.features.rv.th = 1.1, assay = "Spatial", vst.flavor = "v1")

# Save the RDS
saveRDS(seurat_object, ".../GBM_experiment/RDS/gbm.rds")

# Make Variable containing the path to a Spatial directory
spatial_dir <- '.../GBM_experiment/spatial/'

# Read in the tissue_position_list.csv
colData <- read.csv(file.path(spatial_dir, "tissue_positions_list.csv"), header = FALSE)

# Read in the following tibble as the column names in the colData list
colnames(colData) <- c("spot", "in_tissue", "row", "col", "imagerow", "imagecol")

# Add "_1" at the end of each spot in the spot column
colData$spot <- paste0(colData$spot, "_1")

# Without this line, rownames(colData) would just be a list of integers. This line names the rows after the spots
rownames(colData) <- colData$spot

# Remove the spots with 0's in the spots columns
colData <- colData[colData$in_tissue > 0, ]

# Make SingleCellExperiment object
# In SingleCellExperiment objects, rows should represent features (genes, transcripts, genomic regions) and columns should represent cells.
# Create SCE object by transferring SCT assay and metadata (Spatial coordinates)

# Make ST_sce, a SingleCellExperiment object. Stores genes in rows, and tags in columns
ST_sce <- SingleCellExperiment(assays = list(SCT = as.matrix(GetAssayData(seurat_object, layer = "data"))))
ST_sce <- SingleCellExperiment(assays = list(SCT = as.matrix(GetAssayData(seurat_object, layer = "data"))), colData = colData[colnames(ST_sce@assays@data$SCT),])


# ================================== RUN CLUSTERING ================================= #


# Set some parameters of the ST_sce object
metadata(ST_sce)$BayesSpace.data <- list()
metadata(ST_sce)$BayesSpace.data$platform <- "Visium"
metadata(ST_sce)$BayesSpace.data$is.enhanced <- FALSE

# Cluster with BayesSpace
# Compute top 2000 HVGs and 15 PCA using the SCT data

# This line preprocesses the data and will likely give the following error:
# "In regularize.values(x, y, ties, missing(ties), na.rm = na.rm) : collapsing to unique 'x' values"
# But it can be ignored and the rest of the file should work just fine
ST_sce <- spatialPreprocess(ST_sce, platform = "Visium", n.PCs = 15, n.HVGs = 2000, log.normalize = FALSE, assay.type = "SCT", )#skip.PCA = TRUE)

# Fin q for clusters

# Set the seed for the clustering
set.seed(100)

# Run qTune
# qTune() computes the average negative log likelihood for a range of q values over iterations 100:1000, and qPlot() displays the results.
ST_sce <- qTune(ST_sce, qs = seq(2, 20), burn.in = 10, nrep=100)

# Make qPlot
p <- qPlot(ST_sce)
pdf(".../GBM_experiment/Images/BayesQplot.pdf", width = 8, height = 7)
print(p)
dev.off()

# Number of clusters
q = 11

# Number of PCA
d = 15

# Run the clustering 
ST_sce <- spatialCluster(ST_sce, q=q, d=d, init.method = "kmeans", platform = "Visium", nrep = 50000, gamma = 3)
ST_sce$spatial.cluster <- as.factor(ST_sce$spatial.cluster)

# Adjust coloring and create ClusterPlot
palette <- RColorBrewer::brewer.pal(q, name = "Paired")
clusterPlot(ST_sce, palette = palette, size=0.05) + labs(title = "Spot-level clustering")

# Save the RDS
saveRDS(ST_sce, ".../GBM_experiment/RDS/GBM_sce.rds")
