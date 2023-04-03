library(shiny)
library(shinydashboard)
library(shinyjs)
library(tools)
library(Seurat)
library(ggplot2)
library(shinybusy)
library(cowplot)
library(markdown)

# load Seurat object from an rds file and validate it
load_seurat_obj <- function(path){
  errors <- c() #empty vector to store errors
  
  # check if the file has a .rds extension
  if (!tolower(tools::file_ext(path)) == "rds") { # ignore case
    errors <- c(errors, "Invalid rds file")
    return(errors)
  }
  
  # try to read the rds file and check if it's a Seurat object
  tryCatch(
    {
      obj <- readRDS(path)
    },
    error = function(e) {
      errors <- c(errors, "Invalid rds file")
      return(errors)
    }
  )
  if (!inherits(obj, "Seurat")) {
    errors <- c(errors, "File is not a Seurat object")
    return(errors)
  }
  return(obj)
}

# create a UMAP plot for a given metadata column in a Seurat object
create_metadata_umap <- function(obj, col){
  
  # use custom plot if the metadata column are the following
  if (col %in% c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
    col_df <- data.frame(obj@reductions$umap@cell.embeddings, data = obj@meta.data[,col])
    umap <- ggplot(data = col_df) +
      geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(data)), size = 0.01) +
      scale_colour_gradientn(colours = rainbow(7)) +
      labs(title = col) +
      guides(col = guide_colourbar(title = "Log10")) +
      theme_cowplot() +
      theme(plot.title = element_text(hjust = 0.5))
    
  # else use DimPlot function from Seurat 
  } else if (col %in% colnames(obj@meta.data)) {
    umap <- DimPlot(obj, pt.size = .1, label = F, label.size = 4, group.by = col, reduction = "umap")
    
    # return an empty plot with error message if the metadata column doesn't exist
  } else {
    umap <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Metadata column doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(umap)
}

# create a gene expression plot for a given gene in a Seurat object
create_feature_plot <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    FP <- FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
  
    # return an empty plot with error message if the gene doesn't exist
  } else { 
    FP <- ggplot() + 
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(FP)
}