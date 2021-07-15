###############################################################################
##################                Modularity                ###################
###############################################################################

# [SEURAT]

# Description
#
# For each cell, the algorithm checks that the K neighbors belong to the same 
# cluster as this cell (K is k.param given in Seurat FindNeighbors())
# Calculation of the modularity of clustering.
# INPUT 1: Seurat Object
# OUTPUT: modularity score


# Function

modularity = function(object){
  
  neighbors_k = as.list(0)
  for(i in 1:length(derybel_epidermis_seurat@meta.data[["barcodes"]])){
    nn_cells = derybel_epidermis_seurat@graphs[["RNA_snn"]][which(derybel_epidermis_seurat@graphs[["RNA_snn"]][-i,i]!=0),i]
    neighbors_k[[i]] = (names(sort(nn_cells, decreasing = T)[1:20]))
  }
  
  clusfrom = as.numeric(derybel_epidermis_seurat@active.ident)
  
  for(i in 1:length(neighbors_k)){
    neighbors_k[[i]]= as.numeric(derybel_epidermis_seurat@active.ident[neighbors_k[[i]]])
  }
  
  modularity_score = 0
  
  for (i in 1: length(neighbors_k)){
    modularity_score = modularity_score + sum(clusfrom[i] == neighbors_k[[i]])
  }
  
  modularity_score = modularity_score/(length(neighbors_k)*20)
  
  return(modularity_score)
}


# Example
modularity(derybel_epidermis_seurat)
