##############################################################################
################                FeaturePlot_3D                ################
##############################################################################

# [SEURAT]

##### Description #####

# This function is a tool for 3D representation of the metadata of the Seurat object according to the dimensional reduction method used.
#
# INPUT
#  object : Seurat Object
#  reduction.method : "tsne", "umap", "pca" or "ica"
#  feature : gene to show
#  dims : each dimension wanted
# OUTPUT
#  Cells graph



##### Function #####

FeaturePlot_3D = function(object, dims, reduction.method, feature){
  
  library(rgl)
  options(rgl.printRglwidget =TRUE)
  
  color.list = c("lightgrey", "blue")
  
  # Error managment
  if(!class(object)[1] == "Seurat"){
    return(message("ERROR :  Object's class must be 'Seurat'."))
  }
  if((!reduction.method== "tsne") && (!reduction.method== "umap") && (!reduction.method== "pca") && (!reduction.method== "ica")){
    return(message("ERROR : wrong reduction.method name."))
  }
  if(!length(dims) == 3){
    return(message("ERROR :  dims must be 3 digits "))
  }
  if(!feature%in% object@assays[["RNA"]]@counts@Dimnames[[1]]){
    return(message("ERROR : wrong name for feature in feature parameter"))
  }
  
  # 3D plot according to feature
  if(reduction.method == "tsne"){
    
    tSNE_1 = object@reductions[["tsne"]]@cell.embeddings[,dims[1]]
    tSNE_2 = object@reductions[["tsne"]]@cell.embeddings[,dims[2]]
    tSNE_3 = object@reductions[["tsne"]]@cell.embeddings[,dims[3]]
    
    cell.expression = object@assays[["RNA"]]@data[feature,]
    
    fig = plot_ly(x=tSNE_1, y=tSNE_2, z=tSNE_3, type="scatter3d",mode = "markers", color = cell.expression, 
                  colors = color.list, sizes = c(10:10), size = 1)
    
    return(fig)
  }
  
  if(reduction.method == "umap"){
    
    umap_1 = object@reductions[["umap"]]@cell.embeddings[,dims[1]]
    umap_2 = object@reductions[["umap"]]@cell.embeddings[,dims[2]]
    umap_3 = object@reductions[["umap"]]@cell.embeddings[,dims[3]]
    
    cell.expression = object@assays[["RNA"]]@data[feature,]
    
    fig = plot_ly(x=umap_1, y=umap_2, z=umap_3, type="scatter3d",mode = "markers", color = cell.expression, 
                  colors = color.list, sizes = c(10:10), size = 1)
    
    return(fig)
  }
  
  if(reduction.method == "pca"){
    
    pca_1 = object@reductions[["pca"]]@cell.embeddings[,dims[1]]
    pca_2 = object@reductions[["pca"]]@cell.embeddings[,dims[2]]
    pca_3 = object@reductions[["pca"]]@cell.embeddings[,dims[3]]
    
    cell.expression = object@assays[["RNA"]]@data[feature,]
    
    fig = plot_ly(x=pca_1, y=pca_2, z=pca_3, type="scatter3d",mode = "markers", color = cell.expression, 
                  colors = color.list, sizes = c(10:10), size = 1)
    
    return(fig)
  }
  
  if(reduction.method == "ica"){
    
    ica_1 = object@reductions[["ica"]]@cell.embeddings[,dims[1]]
    ica_2 = object@reductions[["ica"]]@cell.embeddings[,dims[2]]
    ica_3 = object@reductions[["ica"]]@cell.embeddings[,dims[3]]
    
    cell.expression = object@assays[["RNA"]]@data[feature,]
    
    fig = plot_ly(x=ica_1, y=ica_2, z=ica_3, type="scatter3d",mode = "markers", color = cell.expression, 
                  colors = color.list, sizes = c(10:10), size = 1)
    
    return(fig)
  }
}

# Example : 
FeaturePlot_3D(object = pbmc, dims = c(1,2,3), reduction.method = "tsne", feature = "AT5G49800")


