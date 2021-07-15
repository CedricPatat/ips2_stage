###############################################################################
##################                Modularity                ###################
###############################################################################

# [MONOCLE 3]

# Description

# For each cell, the algorithm checks that the K neighbors belong to the same 
# cluster as this cell (K given in Monocle 3 cluster_cells ())
# Calculation of the modularity of clustering.
# INPUT 1: cds as Monocle's cds 
# INPUT 2: reduction_method used as red_met (do not forget quotation marks)
# OUTPUT: modularity score


#Function

modularity = function(cds, red_meth) {
  myclusters= clusters(cds, reduction_method = red_meth)
  from = cds@clusters@listData[[red_meth]][["cluster_result"]][["relations"]][["from"]]
  to = cds@clusters@listData[[red_meth]][["cluster_result"]][["relations"]][["to"]]
  clus_from = as.numeric(myclusters[from], reduction_method = red_meth)
  clus_to = as.numeric(myclusters[to], reduction_method = red_meth)
  return(sum(clus_from == clus_to)/length(from))
}

# Example
modularity(derybel_cds, "tsne")
