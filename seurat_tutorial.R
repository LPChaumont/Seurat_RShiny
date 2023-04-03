library(Seurat)
library(R.utils)

# download raw data
url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
file <- basename(url)
download.file(url, destfile = file)
untar(R.utils::gunzip(file))
unlink("pbmc3k_filtered_gene_bc_matrices.tar")

# Follow tutorial on https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data)
pbmc <- CreateSeuratObject(counts = pbmc.data,
                           project = "pbmc3k",
                           min.cells = 3,
                           min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Normalizing data
pbmc <- NormalizeData(pbmc)

# Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling the data
pbmc <- ScaleData(pbmc)

# Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run non-linear dimensional reduction (UMAP/tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Add cells label
pbmc$cell_type <- Idents(pbmc)

# Add random clinical information about the cells : age, sex, sampleID
set.seed(2000)
sex <- sample(c("male", "female"), ncol(pbmc), replace = TRUE)
pbmc$sex <- sex

age <- round(runif(ncol(pbmc), 18, 65))
pbmc$age <- age

sampleID <- paste0("sample", sample(1:10, ncol(pbmc), replace = TRUE))
pbmc$sampleID <- sampleID

# save Seurat object
saveRDS(pbmc, file = "seurat_object.rds")
