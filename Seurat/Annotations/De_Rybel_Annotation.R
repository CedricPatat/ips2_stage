##############################################################################
################              De Rybel Annotation             ################
##############################################################################

# [SEURAT]

##### Description
#
# This script is the tool used to annotate the cell types of the Wendrich dataset.
# The annotation is made from clustering of the entire dataset and subsets of the Nemhauser data in UMAP with Seurat. 
# Cell types are checked using marker genes provided by data from Wendrich et al and Nemhauser et al..
#
# We will save two new columns in the original file. The first column entitled "cell_type_ips2"
# will store the new cell type assignments and the "gene_marker" column stores the marker genes
# of the cluster or sub-cluster which made it possible to identify the cell type.


##### Annotation Script

### Loading libraries
library(Seurat)
library(Matrix)

### Loading data

setwd("~/Bureau/De-Rybel")

derybel_matrix = readMM("GSE141730_matrix.mtx.gz")
derybel_genes = read.table("GSE141730_genes.tsv")
derybel_cells = read.table("GSE141730_barcodes.tsv")

derybel_genes = as.data.frame(derybel_genes[,1])
colnames(derybel_genes) = "genes_names"
rownames(derybel_genes)= derybel_genes[,1]

derybel_cells = cbind(derybel_cells, as.data.frame(matrix(NA, ncol=1, nrow=dim(derybel_cells)[1])), 
                      as.data.frame(matrix(NA, ncol=1, nrow=dim(derybel_cells)[1])))
colnames(derybel_cells) = c("barcodes", "cell_type_ips2", "gene_marker")
rownames(derybel_cells) = derybel_cells[,1]

derybel_matrix = as(derybel_matrix, "dgTMatrix")
colnames(derybel_matrix)= rownames(derybel_cells)
rownames(derybel_matrix) = derybel_genes[,1]


### Clustering of the whole dataset

derybel_all = CreateSeuratObject(counts = derybel_matrix,
                                 project = "Gala Annotation",
                                 assay = "RNA",
                                 names.field = 1,
                                 names.delim = NULL,
                                 meta.data = derybel_cells,
                                 min.cells = 0,
                                 min.features = 0,
                                 row.names = NULL)

# Checking outliers
derybel_all[["percent.mt"]] <- PercentageFeatureSet(derybel_all, pattern = "^ATM")
VlnPlot(derybel_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

selected_derybel_cells = derybel_cells[(derybel_all[["percent.mt"]] < 0.10) & (derybel_all[["nCount_RNA"]] < 3e+5) 
                                       & (derybel_all[["nFeature_RNA"]] > 4000),]

# Annotating outliers cells

derybel_cells[setdiff(rownames(derybel_cells),rownames(selected_derybel_cells)), 2] = "Outliers"

# Rearrangment of the matrix

derybel_filtered_matrix = derybel_matrix[, rownames(selected_derybel_cells)]

### Clustering of the whole dataset

derybel_filtered_all = CreateSeuratObject(counts = derybel_filtered_matrix,
                                          project = "De Rybel Annotation",
                                          assay = "RNA",
                                          names.field = 1,
                                          names.delim = NULL,
                                          meta.data = selected_derybel_cells,
                                          min.cells = 0,
                                          min.features = 0,
                                          row.names = NULL)

# Checking outliers
derybel_filtered_all[["percent.mt"]] <- PercentageFeatureSet(derybel_filtered_all, pattern = "^ATM")
VlnPlot(derybel_filtered_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was done this time.

derybel_filtered_all = NormalizeData(derybel_filtered_all,
                                     normalization.method = "LogNormalize",
                                     scale.factor = 10000,
                                     margin = 1, 
                                     block.size = NULL,
                                     verbose = TRUE)

derybel_filtered_all = FindVariableFeatures(derybel_filtered_all,
                                            assay = NULL,
                                            selection.method = "vst",
                                            loess.span = 0.3,
                                            clip.max = "auto",
                                            num.bin = 20,
                                            binning.method = "equal_width",
                                            nfeatures = 2000,
                                            mean.cutoff = c(0.1, 8),
                                            dispersion.cutoff = c(1, Inf),
                                            verbose = TRUE)

derybel_filtered_all = ScaleData(derybel_filtered_all,
                                 features = NULL,
                                 assay = NULL,
                                 vars.to.regress = NULL,
                                 split.by = NULL,
                                 model.use = "linear",
                                 use.umi = FALSE,
                                 do.scale = TRUE,
                                 do.center = TRUE,
                                 scale.max = 10,
                                 block.size = 1000,
                                 min.cells.to.block = 3000,
                                 verbose = TRUE)

derybel_filtered_all = RunPCA(derybel_filtered_all,
                              assay = NULL,
                              features = NULL,
                              npcs = 120,
                              rev.pca = FALSE,
                              weight.by.var = TRUE,
                              ndims.print = 1:5,
                              nfeatures.print = 30,
                              reduction.name = "pca",
                              reduction.key = "PC_",
                              seed.use = 42,
                              verbose = TRUE)

derybel_filtered_all = RunUMAP(derybel_filtered_all,
                               dims = 1:120,
                               reduction = "pca",
                               features = NULL,
                               graph = NULL,
                               assay = DefaultAssay(object = object),
                               nn.name = NULL,
                               slot = "data",
                               umap.method = "uwot",
                               reduction.model = NULL,
                               return.model = FALSE,
                               n.neighbors = 30L,
                               n.components = 4L,
                               metric = "cosine",
                               n.epochs = NULL,
                               learning.rate = 1,
                               min.dist = 0.3,
                               spread = 1,
                               set.op.mix.ratio = 1,
                               local.connectivity = 1L,
                               repulsion.strength = 1,
                               negative.sample.rate = 5L,
                               a = NULL,
                               b = NULL,
                               uwot.sgd = FALSE,
                               seed.use = 42L,
                               metric.kwds = NULL,
                               angular.rp.forest = FALSE,
                               reduction.name = "umap",
                               reduction.key = "UMAP_",
                               verbose = TRUE)

derybel_filtered_all = FindNeighbors(derybel_filtered_all,
                                     reduction = "umap", 
                                     dims=1:3,
                                     assay = NULL,
                                     features =NULL,
                                     k.param = 20,
                                     return.neighbor = F)

derybel_filtered_all = FindClusters(derybel_filtered_all,
                                    graph.name = NULL,
                                    modularity.fxn = 1,
                                    initial.membership = NULL,
                                    node.sizes = NULL,
                                    resolution = 0.001,
                                    method = "matrix",
                                    algorithm = 3,
                                    n.start = 20,
                                    n.iter = 20,
                                    random.seed = 0,
                                    group.singletons = TRUE,
                                    temp.file.location = NULL,
                                    edge.file.name = NULL,
                                    verbose = TRUE)

DimPlot_3D(object = derybel_filtered_all, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

  # We get 7 clusters. We will study directly cluster 6. Then we will do subsets.
  # Subset1 : clusters 1, 5
  # Subset2 : clusters 2, 0
  # Subset3 : cluster 3
  # Subset4 : cluster 4
  #
  # We do not create a subset for cluster 6 because he is enough small to contain only few cells. Cluster 6 is subset 0.


# Search for marker gene expression clusters

clusters.markers0 = FindAllMarkers(derybel_filtered_all, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Whole dataset Cluster 6
clusters.markers0[which(clusters.markers0[,6] == 6),][1:5,]
specificity(derybel_filtered_all, "AT5G06630")
FeaturePlot_3D(object = derybel_filtered_all, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G06630")
  # TAIR indicates that this gene is specific for the cell type: Trichoblast
  # The calculation of the specificity gives a score of 99.8% for this cluster.
  # CONCLUSION : Trichoblast -> AT5G06630


### Subclustering of the whole dataset (clusters 1 and 5)

cells.subset1 = c(names(derybel_filtered_all@active.ident[which(derybel_filtered_all@active.ident ==1)]),
                  names(derybel_filtered_all@active.ident[which(derybel_filtered_all@active.ident ==5)]))
derybel_cells.subset1 = selected_derybel_cells[cells.subset1,]
derybel_matrix.subset1 = derybel_filtered_matrix[,cells.subset1]


derybel_subset1 = CreateSeuratObject(counts = derybel_matrix.subset1,
                                     project = "De Rybel Annotation",
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = NULL,
                                     meta.data = derybel_cells.subset1,
                                     min.cells = 0,
                                     min.features = 0,
                                     row.names = NULL)

derybel_subset1 = NormalizeData(derybel_subset1,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                margin = 1, 
                                block.size = NULL,
                                verbose = TRUE)

derybel_subset1 = FindVariableFeatures(derybel_subset1,
                                       assay = NULL,
                                       selection.method = "vst",
                                       loess.span = 0.3,
                                       clip.max = "auto",
                                       num.bin = 20,
                                       binning.method = "equal_width",
                                       nfeatures = 2000,
                                       mean.cutoff = c(0.1, 8),
                                       dispersion.cutoff = c(1, Inf),
                                       verbose = TRUE)

derybel_subset1 = ScaleData(derybel_subset1,
                            features = NULL,
                            assay = NULL,
                            vars.to.regress = NULL,
                            split.by = NULL,
                            model.use = "linear",
                            use.umi = FALSE,
                            do.scale = TRUE,
                            do.center = TRUE,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = TRUE)

derybel_subset1 = RunPCA(derybel_subset1,
                         assay = NULL,
                         features = NULL,
                         npcs = 120,
                         rev.pca = FALSE,
                         weight.by.var = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.name = "pca",
                         reduction.key = "PC_",
                         seed.use = 42,
                         verbose = TRUE)

derybel_subset1 = RunUMAP(derybel_subset1,
                          dims = 1:120,
                          reduction = "pca",
                          features = NULL,
                          graph = NULL,
                          assay = DefaultAssay(object = object),
                          nn.name = NULL,
                          slot = "data",
                          umap.method = "uwot",
                          reduction.model = NULL,
                          return.model = FALSE,
                          n.neighbors = 30L,
                          n.components = 4L,
                          metric = "cosine",
                          n.epochs = NULL,
                          learning.rate = 1,
                          min.dist = 0.3,
                          spread = 1,
                          set.op.mix.ratio = 1,
                          local.connectivity = 1L,
                          repulsion.strength = 1,
                          negative.sample.rate = 5L,
                          a = NULL,
                          b = NULL,
                          uwot.sgd = FALSE,
                          seed.use = 42L,
                          metric.kwds = NULL,
                          angular.rp.forest = FALSE,
                          reduction.name = "umap",
                          reduction.key = "UMAP_",
                          verbose = TRUE)

derybel_subset1 = FindNeighbors(derybel_subset1,
                                reduction = "umap", 
                                dims=1:4,
                                assay = NULL,
                                features =NULL,
                                k.param = 20,
                                return.neighbor = F)

derybel_subset1 = FindClusters(derybel_subset1,
                               graph.name = NULL,
                               modularity.fxn = 1,
                               initial.membership = NULL,
                               node.sizes = NULL,
                               resolution = 0.5,
                               method = "matrix",
                               algorithm = 3,
                               n.start = 20,
                               n.iter = 20,
                               random.seed = 0,
                               group.singletons = TRUE,
                               temp.file.location = NULL,
                               edge.file.name = NULL,
                               verbose = TRUE)

DimPlot_3D(object = derybel_subset1, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

  # We get 19 relatively well separated clusters. So we are going to do a biological study of each of these clusters.


# Search for marker gene expression clusters

clusters.markers1 = FindAllMarkers(derybel_subset1, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset1 Cluster 0
clusters.markers1[which(clusters.markers1[,6] == 0),][1:5,]
specificity(derybel_subset1, "AT3G02850")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G02850")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 43.0% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT3G02850

  # Subset1 Cluster 1
clusters.markers1[which(clusters.markers1[,6] == 1),][1:5,]
specificity(derybel_subset1, "AT3G45700")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G45700")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of % for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT3G45700

  # Subset1 Cluster 2
clusters.markers1[which(clusters.markers1[,6] == 2),][1:5,]
specificity(derybel_subset1, "AT1G25530")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G25530")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Late Procambium
  # The calculation of the specificity gives a score of 40.7% for this cluster.
  # CONCLUSION : Late Procambium -> AT1G25530

  # Subset1 Cluster 3
clusters.markers1[which(clusters.markers1[,6] == 3),][1:5,]
specificity(derybel_subset1, "AT3G48185")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G48185")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Early procambium
  # The calculation of the specificity gives a score of 29.9% for this cluster.
  # CONCLUSION : Early procambium -> AT3G48185

  # Subset1 Cluster 4
clusters.markers1[which(clusters.markers1[,6] == 4),][1:5,]
specificity(derybel_subset1, "AT1G14190")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G14190")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Metaxylem
  # The calculation of the specificity gives a score of 84.4% for this cluster.
  # CONCLUSION : Metaxylem -> AT1G14190

  # Subset1 Cluster 5
clusters.markers1[which(clusters.markers1[,6] == 5),][1:5,]
specificity(derybel_subset1, "AT2G02020")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G02020")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Phloem cambium
  # The calculation of the specificity gives a score of 44.5% for this cluster.
  # CONCLUSION : Phloem Cambium -> AT2G02020

  # Subset1 Cluster 6
clusters.markers1[which(clusters.markers1[,6] == 6),][1:5,]
specificity(derybel_subset1, "AT3G12730")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G12730")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Companion Cell
  # The calculation of the specificity gives a score of 55.0% for this cluster.
  # CONCLUSION : Companion Cell -> AT3G12730

  # Subset1 Cluster 7
clusters.markers1[which(clusters.markers1[,6] == 7),][1:5,]
specificity(derybel_subset1, "AT2G26040")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G26040")
  # TAIR indicates that this gene is specific for the cell type: Xylem Pole Pericycle
  # The calculation of the specificity gives a score of 16.3% for this cluster.
  # CONCLUSION : Xylem Pole Pericycle -> AT2G26040

  # Subset1 Cluster 8
clusters.markers1[which(clusters.markers1[,6] == 8),][1:5,]
specificity(derybel_subset1, "AT5G66700")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G66700")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Initials
  # The calculation of the specificity gives a score of 37.1% for this cluster.
  # CONCLUSION : Initials -> AT5G66700

  # Subset1 Cluster 9
clusters.markers1[which(clusters.markers1[,6] == 9),][1:5,]
specificity(derybel_subset1, "AT3G60140")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G60140")
  # TAIR indicates that this gene is specific for the cell type: Pericycle initials
  # The calculation of the specificity gives a score of 68.6% for this cluster.
  # CONCLUSION : Pericycle Initials -> AT3G60140

  # Subset1 Cluster 10
clusters.markers1[which(clusters.markers1[,6] == 10),][1:5,]
specificity(derybel_subset1, "AT2G26040")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G26040")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Xylem Pole Pericycle
  # The calculation of the specificity gives a score of 20.5% for this cluster.
  # CONCLUSION : Xylem Pole Pericycle -> AT2G26040

  # Subset1 Cluster 11
clusters.markers1[which(clusters.markers1[,6] == 11),][1:5,]
specificity(derybel_subset1, "AT1G29200")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G29200")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Protoxylem
  # The calculation of the specificity gives a score of 97.3% for this cluster.
  # CONCLUSION : Protoxylem -> AT1G29200

  # Subset1 Cluster 12
clusters.markers1[which(clusters.markers1[,6] == 12),][1:5,]
specificity(derybel_subset1, "AT3G05140")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G05140")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Quiescent cells
  # The calculation of the specificity gives a score of 96.0% for this cluster.
  # CONCLUSION : Quiescent cells -> AT3G05140

  # Subset1 Cluster 13
clusters.markers1[which(clusters.markers1[,6] == 13),][1:5,]
specificity(derybel_subset1, "AT5G13080")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G13080")
  # TAIR indicates that this gene is specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 97.9% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT5G13080

  # Subset1 Cluster 14
clusters.markers1[which(clusters.markers1[,6] == 14),][1:5,]
specificity(derybel_subset1, "AT3G25100")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G25100")
  # TAIR indicates that this gene is specific for the cell type: Differenciated Initials
  # The calculation of the specificity gives a score of 37.4% for this cluster.
  # CONCLUSION : Differenciated Initials -> AT3G25100

  # Subset1 Cluster 15
clusters.markers1[which(clusters.markers1[,6] == 15),][1:5,]
specificity(derybel_subset1, "AT1G77200")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G77200")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Dividing Lateral Root Cap
  # The calculation of the specificity gives a score of 88.7% for this cluster.
  # CONCLUSION : Dividing Lateral Root Cap -> AT1G77200

  # Subset1 Cluster 16
clusters.markers1[which(clusters.markers1[,6] == 16),][1:5,]
specificity(derybel_subset1, "AT2G02990")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G02990")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 95.5% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT2G02990

  # Subset1 Cluster 17
clusters.markers1[which(clusters.markers1[,6] == 17),][1:5,]
specificity(derybel_subset1, "AT1G11915")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G11915")
  # TAIR indicates that this gene is specific for the cell type: Sieve element
  # The calculation of the specificity gives a score of 87.8% for this cluster.
  # CONCLUSION : Sieve Element -> AT1G11915

  # Subset1 Cluster 18
clusters.markers1[which(clusters.markers1[,6] == 18),][1:5,]
specificity(derybel_subset1, "AT4G16270")
FeaturePlot_3D(object = derybel_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G16270")
  # De Rybel indicates that this gene is specific for the cell type: Endodermis, but we have seen in his paper that he get 2 groups of
  # endodermis cells, and one of the both was mapped with the initias cells. In his atlas, this group of endodermis cell is fusionned
  # with initials cell. So that is why we gonna to annotate this cluster : Initials
  # The calculation of the specificity gives a score of 86.8% for this cluster.
  # CONCLUSION : Initials -> AT4G16270


### Subclustering of the whole dataset (clusters 0 and 2)

cells.subset2 = c(names(derybel_filtered_all@active.ident[which(derybel_filtered_all@active.ident ==0)]),
                  names(derybel_filtered_all@active.ident[which(derybel_filtered_all@active.ident ==2)]))
derybel_cells.subset2 = selected_derybel_cells[cells.subset2,]
derybel_matrix.subset2 = derybel_filtered_matrix[,cells.subset2]


derybel_subset2 = CreateSeuratObject(counts = derybel_matrix.subset2,
                                     project = "De Rybel Annotation",
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = NULL,
                                     meta.data = derybel_cells.subset2,
                                     min.cells = 0,
                                     min.features = 0,
                                     row.names = NULL)

derybel_subset2 = NormalizeData(derybel_subset2,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                margin = 1, 
                                block.size = NULL,
                                verbose = TRUE)

derybel_subset2 = FindVariableFeatures(derybel_subset2,
                                       assay = NULL,
                                       selection.method = "vst",
                                       loess.span = 0.3,
                                       clip.max = "auto",
                                       num.bin = 20,
                                       binning.method = "equal_width",
                                       nfeatures = 2000,
                                       mean.cutoff = c(0.1, 8),
                                       dispersion.cutoff = c(1, Inf),
                                       verbose = TRUE)

derybel_subset2 = ScaleData(derybel_subset2,
                            features = NULL,
                            assay = NULL,
                            vars.to.regress = NULL,
                            split.by = NULL,
                            model.use = "linear",
                            use.umi = FALSE,
                            do.scale = TRUE,
                            do.center = TRUE,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = TRUE)

derybel_subset2 = RunPCA(derybel_subset2,
                         assay = NULL,
                         features = NULL,
                         npcs = 120,
                         rev.pca = FALSE,
                         weight.by.var = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.name = "pca",
                         reduction.key = "PC_",
                         seed.use = 42,
                         verbose = TRUE)

derybel_subset2 = RunUMAP(derybel_subset2,
                          dims = 1:120,
                          reduction = "pca",
                          features = NULL,
                          graph = NULL,
                          assay = DefaultAssay(object = object),
                          nn.name = NULL,
                          slot = "data",
                          umap.method = "uwot",
                          reduction.model = NULL,
                          return.model = FALSE,
                          n.neighbors = 30L,
                          n.components = 4L,
                          metric = "cosine",
                          n.epochs = NULL,
                          learning.rate = 1,
                          min.dist = 0.3,
                          spread = 1,
                          set.op.mix.ratio = 1,
                          local.connectivity = 1L,
                          repulsion.strength = 1,
                          negative.sample.rate = 5L,
                          a = NULL,
                          b = NULL,
                          uwot.sgd = FALSE,
                          seed.use = 42L,
                          metric.kwds = NULL,
                          angular.rp.forest = FALSE,
                          reduction.name = "umap",
                          reduction.key = "UMAP_",
                          verbose = TRUE)

derybel_subset2 = FindNeighbors(derybel_subset2,
                                reduction = "umap", 
                                dims=1:4,
                                assay = NULL,
                                features =NULL,
                                k.param = 20,
                                return.neighbor = F)

derybel_subset2 = FindClusters(derybel_subset2,
                               graph.name = NULL,
                               modularity.fxn = 1,
                               initial.membership = NULL,
                               node.sizes = NULL,
                               resolution = 0.5,
                               method = "matrix",
                               algorithm = 3,
                               n.start = 20,
                               n.iter = 20,
                               random.seed = 0,
                               group.singletons = TRUE,
                               temp.file.location = NULL,
                               edge.file.name = NULL,
                               verbose = TRUE)

DimPlot_3D(object = derybel_subset2, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

  # We get 20 relatively well separated clusters. So we are going to do a biological study of each of these clusters.


# Search for marker gene expression clusters

clusters.markers2 = FindAllMarkers(derybel_subset2, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset2 Cluster 0
clusters.markers2[which(clusters.markers2[,6] == 0),][1:5,]
specificity(derybel_subset2, "AT1G67850")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G67850")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Differenciating Lateral Root Cap
  # The calculation of the specificity gives a score of 28.2% for this cluster.
  # CONCLUSION : Differenciating Lateral Root Cap -> AT1G67850

  # Subset2 Cluster 1
clusters.markers2[which(clusters.markers2[,6] == 1),][1:5,]
specificity(derybel_subset2, "AT4G16230")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G16230")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Mature Lateral Root Cap
  # The calculation of the specificity gives a score of 52.3% for this cluster.
  # CONCLUSION : Mature Lateral Root Cap -> AT4G16230

  # Subset2 Cluster 2
clusters.markers2[which(clusters.markers2[,6] == 2),][1:5,]
specificity(derybel_subset2, "AT5G45210")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Endodermis Initials
  # The calculation of the specificity gives a score of 55.4% for this cluster.
  # CONCLUSION : Endodermis Initials -> AT5G45210

  # Subset2 Cluster 3
clusters.markers2[which(clusters.markers2[,6] == 3),][1:5,]
specificity(derybel_subset2, "AT4G00480")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G00480")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Epidermis Initials
  # The calculation of the specificity gives a score of 36.5% for this cluster.
  # CONCLUSION : Epidermis Initials -> AT4G00480

  # Subset2 Cluster 4
clusters.markers2[which(clusters.markers2[,6] == 4),][1:5,]
specificity(derybel_subset2, "AT5G40330")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G40330")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Atrichoblast
  # The calculation of the specificity gives a score of 44.9% for this cluster.
  # CONCLUSION : Atrichoblast -> AT5G40330

  # Subset2 Cluster 5
clusters.markers2[which(clusters.markers2[,6] == 5),][1:5,]
specificity(derybel_subset2, "AT4G35030")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G35030")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Young Lateral Root Cap
  # The calculation of the specificity gives a score of 23.3% for this cluster.
  # CONCLUSION : Young Lateral Root Cap -> AT4G35030

  # Subset2 Cluster 6
clusters.markers2[which(clusters.markers2[,6] == 6),][1:5,]
specificity(derybel_subset2, "AT3G53650")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G53650")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Initial LRC
  # The calculation of the specificity gives a score of 65.6% for this cluster.
  # CONCLUSION : Inital Lateral Root Cap -> AT3G53650

  # Subset2 Cluster 7
clusters.markers2[which(clusters.markers2[,6] == 7),][1:5,]
specificity(derybel_subset2, "AT4G12390")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G12390")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Initials
  # The calculation of the specificity gives a score of 47.9% for this cluster.
  # CONCLUSION : Initials -> AT4G12390

  # Subset2 Cluster 8
clusters.markers2[which(clusters.markers2[,6] == 8),][1:5,]
specificity(derybel_subset2, "AT5G50090")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G50090")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Initials
  # The calculation of the specificity gives a score of 50.4% for this cluster.
  # CONCLUSION : Initials -> AT5G50090

  # Subset2 Cluster 9
clusters.markers2[which(clusters.markers2[,6] == 9),][1:5,]
specificity(derybel_subset2, "AT5G57980")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G57980")
# De Rybel's atlas indicates that this gene is specific for the cell type: Epidermis Initial
# The calculation of the specificity gives a score of 29.1% for this cluster.
# CONCLUSION : Epidermis Initial -> AT5G57980

  # Subset2 Cluster 10
clusters.markers2[which(clusters.markers2[,6] == 10),][1:5,]
specificity(derybel_subset2, "AT1G17400")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G17400")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Differenciated Columella
  # The calculation of the specificity gives a score of 93.8% for this cluster.
  # CONCLUSION : Differenciated Columella -> AT1G17400

  # Subset2 Cluster 11
clusters.markers2[which(clusters.markers2[,6] == 11),][1:5,]
specificity(derybel_subset2, "AT2G38160")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G38160")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Dividing Lateral Root Cap
  # The calculation of the specificity gives a score of 56.7% for this cluster.
  # CONCLUSION : Dividing Lateral Root Cap -> AT2G38160

  # Subset2 Cluster 12
clusters.markers2[which(clusters.markers2[,6] == 12),][1:5,]
specificity(derybel_subset2, "AT5G44480")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G44480")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Trichoblast
  # The calculation of the specificity gives a score of 83.5% for this cluster.
  # CONCLUSION : Trichoblast -> AT5G44480

  # Subset2 Cluster 13
clusters.markers2[which(clusters.markers2[,6] == 13),][1:5,]
specificity(derybel_subset2, "AT2G42110")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G42110")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Young Dividing
  # The calculation of the specificity gives a score of 59.3% for this cluster.
  # CONCLUSION : Young Dividing -> AT2G42110

  # Subset2 Cluster 14
clusters.markers2[which(clusters.markers2[,6] == 14),][1:5,]
specificity(derybel_subset2, "AT1G22500")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G22500")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Trichoblast
  # The calculation of the specificity gives a score of 49.5% for this cluster.
  # CONCLUSION : Trichoblast -> AT1G22500

  # Subset2 Cluster 15
clusters.markers2[which(clusters.markers2[,6] == 15),][1:5,]
specificity(derybel_subset2, "AT1G35230")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G35230")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Columella Initial
  # The calculation of the specificity gives a score of 67.2% for this cluster.
  # CONCLUSION : Columella Initial -> AT1G35230

  # Subset2 Cluster 16
clusters.markers2[which(clusters.markers2[,6] == 16),][1:5,]
specificity(derybel_subset2, "AT1G09750")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G09750")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Differenciating Cortex
  # The calculation of the specificity gives a score of 97.7% for this cluster.
  # CONCLUSION : Differenciating cortex -> AT1G09750

  # Subset2 Cluster 17
clusters.markers2[which(clusters.markers2[,6] == 17),][1:5,]
specificity(derybel_subset2, "AT4G18450")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G18450")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Differenciating Lateral Root Cap
  # The calculation of the specificity gives a score of 92.1% for this cluster.
  # CONCLUSION : Differenciating Lateral Root Cap -> AT4G18450

  # Subset2 Cluster 18
clusters.markers2[which(clusters.markers2[,6] == 18),][1:5,]
specificity(derybel_subset2, "AT3G17680")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G17680")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Dividing
  # The calculation of the specificity gives a score of 78.0% for this cluster.
  # CONCLUSION : Dividing -> AT3G17680

  # Subset2 Cluster 19
clusters.markers2[which(clusters.markers2[,6] == 19),][1:5,]
specificity(derybel_subset2, "AT2G16385")
FeaturePlot_3D(object = derybel_subset2, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G16385")
  # De Rybel's atlas indicates that this gene is specific for the cell type: Phloem Cambium
  # The calculation of the specificity gives a score of 73.3% for this cluster.
  # CONCLUSION : Phloem Cambium -> AT2G16385


### Subclustering of the whole dataset (cluster 3)

cells.subset3 = names(derybel_filtered_all@active.ident[which(derybel_filtered_all@active.ident ==3)])
derybel_cells.subset3 = selected_derybel_cells[cells.subset3,]
derybel_matrix.subset3 = derybel_filtered_matrix[,cells.subset3]


derybel_subset3 = CreateSeuratObject(counts = derybel_matrix.subset3,
                                     project = "De Rybel Annotation",
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = NULL,
                                     meta.data = derybel_cells.subset3,
                                     min.cells = 0,
                                     min.features = 0,
                                     row.names = NULL)

derybel_subset3 = NormalizeData(derybel_subset3,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                margin = 1, 
                                block.size = NULL,
                                verbose = TRUE)

derybel_subset3 = FindVariableFeatures(derybel_subset3,
                                       assay = NULL,
                                       selection.method = "vst",
                                       loess.span = 0.3,
                                       clip.max = "auto",
                                       num.bin = 20,
                                       binning.method = "equal_width",
                                       nfeatures = 2000,
                                       mean.cutoff = c(0.1, 8),
                                       dispersion.cutoff = c(1, Inf),
                                       verbose = TRUE)

derybel_subset3 = ScaleData(derybel_subset3,
                            features = NULL,
                            assay = NULL,
                            vars.to.regress = NULL,
                            split.by = NULL,
                            model.use = "linear",
                            use.umi = FALSE,
                            do.scale = TRUE,
                            do.center = TRUE,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = TRUE)

derybel_subset3 = RunPCA(derybel_subset3,
                         assay = NULL,
                         features = NULL,
                         npcs = 120,
                         rev.pca = FALSE,
                         weight.by.var = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.name = "pca",
                         reduction.key = "PC_",
                         seed.use = 42,
                         verbose = TRUE)

derybel_subset3 = RunUMAP(derybel_subset3,
                          dims = 1:120,
                          reduction = "pca",
                          features = NULL,
                          graph = NULL,
                          assay = DefaultAssay(object = object),
                          nn.name = NULL,
                          slot = "data",
                          umap.method = "uwot",
                          reduction.model = NULL,
                          return.model = FALSE,
                          n.neighbors = 30L,
                          n.components = 4L,
                          metric = "cosine",
                          n.epochs = NULL,
                          learning.rate = 1,
                          min.dist = 0.3,
                          spread = 1,
                          set.op.mix.ratio = 1,
                          local.connectivity = 1L,
                          repulsion.strength = 1,
                          negative.sample.rate = 5L,
                          a = NULL,
                          b = NULL,
                          uwot.sgd = FALSE,
                          seed.use = 42L,
                          metric.kwds = NULL,
                          angular.rp.forest = FALSE,
                          reduction.name = "umap",
                          reduction.key = "UMAP_",
                          verbose = TRUE)

derybel_subset3 = FindNeighbors(derybel_subset3,
                                reduction = "umap", 
                                dims=1:4,
                                assay = NULL,
                                features =NULL,
                                k.param = 20,
                                return.neighbor = F)

derybel_subset3 = FindClusters(derybel_subset3,
                               graph.name = NULL,
                               modularity.fxn = 1,
                               initial.membership = NULL,
                               node.sizes = NULL,
                               resolution = 0.8,
                               method = "matrix",
                               algorithm = 3,
                               n.start = 20,
                               n.iter = 20,
                               random.seed = 0,
                               group.singletons = TRUE,
                               temp.file.location = NULL,
                               edge.file.name = NULL,
                               verbose = TRUE)

DimPlot_3D(object = derybel_subset3, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# We get 5 relatively well separated clusters. So we are going to do a biological study of each of these clusters.
# Clusters 0 and 1 seem to be the same cell type


# Search for marker gene expression clusters

clusters.markers3 = FindAllMarkers(derybel_subset3, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

# Subset3 Cluster 0
clusters.markers3[which(clusters.markers3[,6] == 0),][1:5,]
specificity(derybel_subset3, "AT5G06360")
FeaturePlot_3D(object = derybel_subset3, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G06360")
# De Rybel's atlas indicates that this gene is specific for the cell type: Ambiguous
# The calculation of the specificity gives a score of 27.9% for this cluster.
# CONCLUSION : Ambiguous -> AT5G06360

# Subset3 Cluster 1
clusters.markers3[which(clusters.markers3[,6] == 1),][1:5,]
specificity(derybel_subset3, "AT5G60390")
FeaturePlot_3D(object = derybel_subset3, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G60390")
# De Rybel's atlas indicates that this gene is specific for the cell type: Ambiguous
# The calculation of the specificity gives a score of 28.9% for this cluster.
# CONCLUSION :  Ambiguous-> AT5G60390

# Subset3 Cluster 2
clusters.markers3[which(clusters.markers3[,6] == 2),][1:5,]
specificity(derybel_subset3, "AT5G07990")
FeaturePlot_3D(object = derybel_subset3, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G07990")
# De Rybel's atlas indicates that this gene is specific for the cell type: Intermediate Cortex
# The calculation of the specificity gives a score of 95.7% for this cluster.
# CONCLUSION : Intermediate Cortex -> AT5G07990

# Subset3 Cluster 3
clusters.markers3[which(clusters.markers3[,6] == 3),][1:5,]
specificity(derybel_subset3, "AT2G45180")
FeaturePlot_3D(object = derybel_subset3, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G45180")
# De Rybel's atlas indicates that this gene is specific for the cell type: Intermediate Cortex
# The calculation of the specificity gives a score of 35.0% for this cluster.
# CONCLUSION : Intermediate Cortex -> AT2G45180

# Subset3 Cluster 4
clusters.markers3[which(clusters.markers3[,6] == 4),][1:5,]
specificity(derybel_subset3, "AT1G07920")
FeaturePlot_3D(object = derybel_subset3, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G07920")
# De Rybel's atlas indicates that this gene is specific for the cell type: Ambiguous
# The calculation of the specificity gives a score of 24.1% for this cluster.
# CONCLUSION : Ambiguous -> AT1G07920



#########################

### Subclustering of the whole dataset (cluster 3)

cells.subset4 = names(derybel_filtered_all@active.ident[which(derybel_filtered_all@active.ident == 4)])
derybel_cells.subset4 = selected_derybel_cells[cells.subset4,]
derybel_matrix.subset4 = derybel_filtered_matrix[,cells.subset4]


derybel_subset4 = CreateSeuratObject(counts = derybel_matrix.subset4,
                                     project = "De Rybel Annotation",
                                     assay = "RNA",
                                     names.field = 1,
                                     names.delim = NULL,
                                     meta.data = derybel_cells.subset4,
                                     min.cells = 0,
                                     min.features = 0,
                                     row.names = NULL)

derybel_subset4 = NormalizeData(derybel_subset4,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                margin = 1, 
                                block.size = NULL,
                                verbose = TRUE)

derybel_subset4 = FindVariableFeatures(derybel_subset4,
                                       assay = NULL,
                                       selection.method = "vst",
                                       loess.span = 0.3,
                                       clip.max = "auto",
                                       num.bin = 20,
                                       binning.method = "equal_width",
                                       nfeatures = 2000,
                                       mean.cutoff = c(0.1, 8),
                                       dispersion.cutoff = c(1, Inf),
                                       verbose = TRUE)

derybel_subset4 = ScaleData(derybel_subset4,
                            features = NULL,
                            assay = NULL,
                            vars.to.regress = NULL,
                            split.by = NULL,
                            model.use = "linear",
                            use.umi = FALSE,
                            do.scale = TRUE,
                            do.center = TRUE,
                            scale.max = 10,
                            block.size = 1000,
                            min.cells.to.block = 3000,
                            verbose = TRUE)

derybel_subset4 = RunPCA(derybel_subset4,
                         assay = NULL,
                         features = NULL,
                         npcs = 120,
                         rev.pca = FALSE,
                         weight.by.var = TRUE,
                         ndims.print = 1:5,
                         nfeatures.print = 30,
                         reduction.name = "pca",
                         reduction.key = "PC_",
                         seed.use = 42,
                         verbose = TRUE)

derybel_subset4 = RunUMAP(derybel_subset4,
                          dims = 1:120,
                          reduction = "pca",
                          features = NULL,
                          graph = NULL,
                          assay = DefaultAssay(object = object),
                          nn.name = NULL,
                          slot = "data",
                          umap.method = "uwot",
                          reduction.model = NULL,
                          return.model = FALSE,
                          n.neighbors = 30L,
                          n.components = 4L,
                          metric = "cosine",
                          n.epochs = NULL,
                          learning.rate = 1,
                          min.dist = 0.3,
                          spread = 1,
                          set.op.mix.ratio = 1,
                          local.connectivity = 1L,
                          repulsion.strength = 1,
                          negative.sample.rate = 5L,
                          a = NULL,
                          b = NULL,
                          uwot.sgd = FALSE,
                          seed.use = 42L,
                          metric.kwds = NULL,
                          angular.rp.forest = FALSE,
                          reduction.name = "umap",
                          reduction.key = "UMAP_",
                          verbose = TRUE)

derybel_subset4 = FindNeighbors(derybel_subset4,
                                reduction = "umap", 
                                dims=1:4,
                                assay = NULL,
                                features =NULL,
                                k.param = 20,
                                return.neighbor = F)

derybel_subset4 = FindClusters(derybel_subset4,
                               graph.name = NULL,
                               modularity.fxn = 1,
                               initial.membership = NULL,
                               node.sizes = NULL,
                               resolution = 0.8,
                               method = "matrix",
                               algorithm = 3,
                               n.start = 20,
                               n.iter = 20,
                               random.seed = 0,
                               group.singletons = TRUE,
                               temp.file.location = NULL,
                               edge.file.name = NULL,
                               verbose = TRUE)

DimPlot_3D(object = derybel_subset4, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# We get 6 relatively well separated clusters. So we are going to do a biological study of each of these clusters.


# Search for marker gene expression clusters

clusters.markers4 = FindAllMarkers(derybel_subset4, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

# Subset4 Cluster 0
clusters.markers4[which(clusters.markers4[,6] == 0),][1:5,]
specificity(derybel_subset4, "AT4G12520")
FeaturePlot_3D(object = derybel_subset4, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G12520")
# De Rybel's atlas indicates that this gene is specific for the cell type: Intermediate endodermis
# The calculation of the specificity gives a score of 32.7% for this cluster.
# CONCLUSION : Intermediate endordermis -> AT4G12520

# Subset4 Cluster 1
clusters.markers4[which(clusters.markers4[,6] == 1),][1:5,]
specificity(derybel_subset4, "AT3G19820")
FeaturePlot_3D(object = derybel_subset4, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G19820")
# De Rybel's atlas indicates that this gene is specific for the cell type: Endodermis Initials
# The calculation of the specificity gives a score of 46.4% for this cluster.
# CONCLUSION : Endodermis initials -> AT3G19820

# Subset4 Cluster 2
clusters.markers4[which(clusters.markers4[,6] == 2),][1:5,]
specificity(derybel_subset4, "AT5G44160")
FeaturePlot_3D(object = derybel_subset4, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G44160")
# De Rybel's atlas indicates that this gene is specific for the cell type: 
# The calculation of the specificity gives a score of 92.1% for this cluster.
# CONCLUSION : Ambiguous -> AT5G44160

# Subset4 Cluster 3
clusters.markers4[which(clusters.markers4[,6] == 3),][1:5,]
specificity(derybel_subset4, "AT2G41110")
FeaturePlot_3D(object = derybel_subset4, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G41110")
# De Rybel's atlas indicates that this gene is specific for the cell type: Intermediate Endodermis
# The calculation of the specificity gives a score of 34.3% for this cluster.
# CONCLUSION : Intermediate endodermis -> AT2G41110

# Subset4 Cluster 4
clusters.markers4[which(clusters.markers4[,6] == 4),][1:5,]
specificity(derybel_subset4, "AT1G65690")
FeaturePlot_3D(object = derybel_subset4, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G65690")
# De Rybel's atlas indicates that this gene is specific for the cell type: Differenciated Endodermis
# The calculation of the specificity gives a score of 96.6% for this cluster.
# CONCLUSION : Differenciated Endodermis -> AT1G65690

# Subset4 Cluster 5
clusters.markers4[which(clusters.markers4[,6] == 5),][1:5,]
specificity(derybel_subset4, "AT1G53580")
FeaturePlot_3D(object = derybel_subset4, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G53580")
# De Rybel's atlas indicates that this gene is specific for the cell type: Intermediate Endodermis
# The calculation of the specificity gives a score of 51.7% for this cluster.
# CONCLUSION : Intermediate Endodermis -> AT1G53580


### Cell Annotation

## Subset 0
cells.cluster.subset0 = derybel_filtered_all@active.ident
derybel_cells[names(cells.cluster.subset0[which(cells.cluster.subset0 == 6)]),2] = "Trichoblast"
derybel_cells[names(cells.cluster.subset0[which(cells.cluster.subset0 == 6)]),3] = "AT5G06630"
## Subset 1
cells.cluster.subset1 = derybel_subset1@active.ident
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 0)]),2] = "Phloem_Pole_Pericycle"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 0)]),3] = "AT3G02850"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 1)]),2] = "Phloem_Pole_Pericycle"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 1)]),3] = "AT3G45700"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 2)]),2] = "Late_Procambium"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 2)]),3] = "AT1G25530"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 3)]),2] = "Early_Procambium"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 3)]),3] = "AT3G48185"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 4)]),2] = "Metaxylem"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 4)]),3] = "AT1G14190"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 5)]),2] = "Phloem_Cambium"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 5)]),3] = "AT2G02020"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 6)]),2] = "Campanion_Cell"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 6)]),3] = "AT3G12730"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 7)]),2] = "Xylem_Pole_Pericycle"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 7)]),3] = "AT2G26040"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 8)]),2] = "Initial"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 8)]),3] = "AT5G66700"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 9)]),2] = "Pericycle_Initial"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 9)]),3] = "AT3G60140"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 10)]),2] = "Xylem_Pole_Pericycle"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 10)]),3] = "AT2G26040"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 11)]),2] = "Protoxylem"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 11)]),3] = "AT1G29200"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 12)]),2] = "Quiescent_Cell"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 12)]),3] = "AT3G05140"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 13)]),2] = "Phloem_Pole_Pericycle"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 13)]),3] = "AT5G13080"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 14)]),2] = "Differenciated_Initial"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 14)]),3] = "AT3G25100"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 15)]),2] = "Dividing_Lateral_Root_Cap"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 15)]),3] = "AT1G77200"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 16)]),2] = "Phloem_Pole_Pericycle"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 16)]),3] = "AT2G02990"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 17)]),2] = "Sieve_Element"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 17)]),3] = "At1G11915"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 18)]),2] = "Initial"
derybel_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 18)]),3] = "AT4G16270"

## Subset 2
cells.cluster.subset2 = derybel_subset2@active.ident
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 0)]),2] = "Differenciating_Lateral_root_Cap"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 0)]),3] = "AT1G67850"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 1)]),2] = "Mature_Lateral_Root_Cap"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 1)]),3] = "AT4G16230"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 2)]),2] = "Endodermis_Initias"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 2)]),3] = "AT5G45210"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 3)]),2] = "Epidermis_Initial"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 3)]),3] = "AT4G00480"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 4)]),2] = "Atrichoblast"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 4)]),3] = "AT5G40330"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 5)]),2] = "Young_Lateral_root_Cap"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 5)]),3] = "AT4G35030"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 6)]),2] = "Initial_Lateral_Root_Cap"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 6)]),3] = "AT3G53650"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 7)]),2] = "Initial"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 7)]),3] = "AT4G12390"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 8)]),2] = "Initial"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 8)]),3] = "At5G50090"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 9)]),2] = "Epidermis_Initial"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 9)]),3] = "AT5G57980"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 10)]),2] = "Differenciated_Columella"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 10)]),3] = "AT1G17400"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 11)]),2] = "Dividing_Lateral_Root_Cap"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 11)]),3] = "AT2G38160"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 12)]),2] = "Trichoblast"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 12)]),3] = "AT5G44480"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 13)]),2] = "Young_Dividing"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 13)]),3] = "AT2G42110"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 14)]),2] = "Trichoblast"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 14)]),3] = "AT1G22500"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 15)]),2] = "Columella_Initial"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 15)]),3] = "AT1G35230"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 16)]),2] = "Differenciating_Cortex"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 16)]),3] = "AT1G09750"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 17)]),2] = "Differenciating_Lateral_Root_Cap"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 17)]),3] = "AT4G18450"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 18)]),2] = "Dividing"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 18)]),3] = "AT3G17680"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 19)]),2] = "Phloem_Cambium"
derybel_cells[names(cells.cluster.subset2[which(cells.cluster.subset2 == 19)]),3] = "At2G16385"

## Subset 3
cells.cluster.subset3 = derybel_subset3@active.ident
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 0)]),2] = "Ambiguous"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 0)]),3] = "AT5G06360"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 1)]),2] = "Ambiguous"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 1)]),3] = "AT5G60390"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 2)]),2] = "Intermediate_Cortex"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 2)]),3] = "AT5G07990"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 3)]),2] = "Intermediate_Cortex"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 3)]),3] = "AT2G45180"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 4)]),2] = "Ambiguous"
derybel_cells[names(cells.cluster.subset3[which(cells.cluster.subset3 == 4)]),3] = "AT1G07920"

## Subset 4
cells.cluster.subset4 = derybel_subset4@active.ident
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 0)]),2] = "Intermediate_Endodermis"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 0)]),3] = "AT4G12520"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 1)]),2] = "Endodermis_Initials"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 1)]),3] = "AT3G19820"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 2)]),2] = "Ambiguous"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 2)]),3] = "AT5G44160"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 3)]),2] = "Intermediate_Endodermis"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 3)]),3] = "AT2G41110"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 4)]),2] = "Differenciated_Endodermis"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 4)]),3] = "AT1G65690"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 5)]),2] = "Intermediate_Endodermis"
derybel_cells[names(cells.cluster.subset4[which(cells.cluster.subset4 == 5)]),3] = "AT1G53580"


### Storing new annotations
write.table(derybel_cells, "GSE141730_barcodes_annotated.tsv")
  # To view the new cell type assignment, use the group.by = "cell_type_ips2" parameter.



