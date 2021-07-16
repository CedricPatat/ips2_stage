##############################################################################
#################                Specificity                 #################
##############################################################################

# [SEURAT]

# This function is a tool for calculating the specificity of one or more features for a clustering data.
#
# INPUT : 
#  object : Seurat Object
#  features : all the features to analyze
# OUTPUT : data frame
#  col 1 : feature name
#  col 2 : cluster's number
#  col 3 : feature's specificity for the cluster. The specificity is normalized. For a given feature, the sum of the specificities is 1
#  col 4 : the frequency of appearance of the feature in the cells of the cluster
#  col 5 : the p-value adjusted calculated by the function FindAllMarkers() [Seurat]
#  col 6 : a star indicates the best specificity for a given feature only if the maximum is above 0.5


# Function 

specificity = function(object, features){
  
  # Error management
  if(!class(object)[1] == "Seurat"){
    return(message("ERROR :  object's class must be 'Seurat'."))
  }
  for(nfeature in 1:length(features)){
    if(! features[nfeature] %in% derybel_epidermis_seurat@assays[["RNA"]]@counts@Dimnames[[1]]){
      return(message("ERROR : wrong name for feature in feature parameter"))
    }
  }
  
  # Computing specificity
  specificity.df = as.data.frame(matrix("" , ncol=6 , nrow = (length(features) * nlevels(object@active.ident))))
  colnames(specificity.df) = c("Feature_name","Cluster_number", "Specificity_inter_clusters", "Frequency_intra_cluster", "p-value_adjusted","")
  specificity.df[,1] = rep(features, each = nlevels(object@active.ident))
  
  myline = 1 
  
  for(nfeature in 1 : length(features)){
    p_val_adj=FindAllMarkers(object,
                             features= features[nfeature],
                             min.cells.feature = 0,
                             min.cells.group = 0,
                             min.pct = 0,
                             logfc.threshold = 0,
                             verbose = F)
    
    b_line = as.numeric(myline)
    
    for(i in 0:(nlevels(object@active.ident)-1)){
      cells.cluster = which(object@meta.data[["seurat_clusters"]]==i)
      specificity.df[myline,2] = i
      specificity.df[myline,3] = mean(t(t(object@assays[["RNA"]]@counts[features[nfeature], cells.cluster]) / 
                                        (colSums(object@assays[["RNA"]]@counts[, cells.cluster]))))
      specificity.df[myline,4] = sum(object@assays[["RNA"]]@counts[features[nfeature], cells.cluster] != 0) / length(cells.cluster)
      
      if( i %in% p_val_adj[,6]){
        specificity.df[myline,5] = p_val_adj[which(p_val_adj[,6] == i),5]
      }
      myline = myline+1
    }
    specificity.df[b_line:(myline-1),3] = as.numeric(specificity.df[b_line:(myline-1),3])/ (sum(as.numeric(specificity.df[b_line:(myline-1),3])))
    
    if(max(specificity.df[b_line:(myline-1),3]) > 0.5){
      specificity.df[b_line -1 + which.max(specificity.df[b_line:(myline-1),3]),6]="*"
    }
  }
  return(specificity.df)
}

# Example
specificity(object = derybel_epidermis_seurat, features = c("AT2G14900", "AT3G53420"))