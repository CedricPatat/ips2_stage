##############################################################################
################               Gala Annotation                ################
##############################################################################

# [SEURAT]

##### Description
#
# This script is the tool used to annotate the cell types of the Nemhauser dataset.
# The annotation is made from clustering of the entire dataset and subsets of the Nemhauser data in UMAP with Seurat. 
# Cell types are checked using marker genes provided by data from Wendrich et al.
#
# We will save two new columns in the original file. The first column entitled "cell_type_ips2"
# will store the new cell type assignments and the "gene_marker" column stores the marker genes
# of the cluster or sub-cluster which made it possible to identify the cell type.


##### Annotation Script

### Loading libraries
library(Seurat)
library(Matrix)

### Loading data

setwd("~/Bureau/Nemhauser")

gala_matrix = readMM("GSE158761_matrix.mtx.gz")
gala_genes = read.table("GSE158761_genes.tsv")
gala_cells = read.table("GSE158761_barcodes.tsv")

gala_genes = as.data.frame(gala_genes[,1])
colnames(gala_genes) = "genes_names"
rownames(gala_genes)= gala_genes[,1]

gala_cells = cbind(gala_cells, as.data.frame(matrix(NA, ncol=1, nrow=dim(gala_cells)[1])), 
                   as.data.frame(matrix(NA, ncol=1, nrow=dim(gala_cells)[1])))
colnames(gala_cells) = c("barcodes", "hour", "cell_type_gala", "cell_type_ips2", "gene_marker")
rownames(gala_cells) = gala_cells[,1]

gala_matrix = as(gala_matrix, "dgTMatrix")
colnames(gala_matrix)= rownames(gala_cells)
rownames(gala_matrix) = gala_genes[,1]


### Clustering of the whole dataset

gala_all = CreateSeuratObject(counts = gala_matrix,
                                 project = "Gala Annotation",
                                 assay = "RNA",
                                 names.field = 1,
                                 names.delim = NULL,
                                 meta.data = gala_cells,
                                 min.cells = 0,
                                 min.features = 0,
                                 row.names = NULL)

# Checking outliers
gala_all[["percent.mt"]] <- PercentageFeatureSet(gala_all, pattern = "^ATM")
VlnPlot(gala_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_all = NormalizeData(gala_all,
                         normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                         scale.factor = 10000,
                         margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                         block.size = NULL,
                         verbose = TRUE)

gala_all = FindVariableFeatures(gala_all,
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

gala_all = ScaleData(gala_all,
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

gala_all = RunPCA(gala_all,
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

gala_all = RunUMAP(gala_all,
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

gala_all = FindNeighbors(gala_all,
                         reduction = "umap", 
                         dims=1:4,
                         assay = NULL,
                         features =NULL,
                         k.param = 20,
                         return.neighbor = F)

gala_all = FindClusters(gala_all,
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

DimPlot_3D(object = gala_all, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

  # We got 3 clusters. We have therefore divided the entire dataset into two, a first subset
  # consisting of clusters 0 and 2 (supposedly stele + cortex + endodermis) and a second corresponding
  # to cluster 1 (supposedly epidermis + atrichoblast + trichoblast + columella + lateral root cap).
  # Cell types are supposed thanks to Galaés annotations.


### Subclustering of the whole dataset (cluster 1)

cells.subset1 = names(gala_all@active.ident[which(gala_all@active.ident ==1)])
gala_cells.subset1 = gala_cells[cells.subset1,]
gala_matrix.subset1 = gala_matrix[,cells.subset1]

gala_subset1 = CreateSeuratObject(counts = gala_matrix.subset1,
                                   project = "Gala Annotation",
                                   assay = "RNA",
                                   names.field = 1,
                                   names.delim = NULL,
                                   meta.data = gala_cells.subset1,
                                   min.cells = 0,
                                   min.features = 0,
                                   row.names = NULL)

# Checking outliers
gala_subset1[["percent.mt"]] <- PercentageFeatureSet(gala_subset1, pattern = "^ATM")
VlnPlot(gala_subset1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset1 = NormalizeData(gala_subset1,
                             normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                             scale.factor = 10000,
                             margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                             block.size = NULL,
                             verbose = TRUE)

gala_subset1 = FindVariableFeatures(gala_subset1,
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

gala_subset1 = ScaleData(gala_subset1,
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

gala_subset1 = RunPCA(gala_subset1,
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

gala_subset1 = RunUMAP(gala_subset1,
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

gala_subset1 = FindNeighbors(gala_subset1,
                             reduction = "umap", 
                             dims=1:4,
                             assay = NULL,
                             features =NULL,
                             k.param = 20,
                             return.neighbor = F)

gala_subset1 = FindClusters(gala_subset1,
                            graph.name = NULL,
                            modularity.fxn = 1,
                            initial.membership = NULL,
                            node.sizes = NULL,
                            resolution = 0.05,
                            method = "matrix",
                            algorithm = 3,
                            n.start = 20,
                            n.iter = 20,
                            random.seed = 0,
                            group.singletons = TRUE,
                            temp.file.location = NULL,
                            edge.file.name = NULL,
                            verbose = TRUE)

DimPlot_3D(object = gala_subset1, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")


# Search for marker gene expression clusters

clusters.markers1 = FindAllMarkers(gala_subset1, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset1 Cluster 0 and 4
  # It's hard to conclude. A sub-clustering will be carried out.

  # Subset1 Cluster 1
clusters.markers[which(clusters.markers[,6] == 1),][1:5,]
specificity(gala_subset1, "AT2G33790")
FeaturePlot_3D(object = gala_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G33790")
  # TAIR indicates that this gene is very specific for the cell type: Atrichoblast
  # The calculation of the specificity gives a score of 72.9% for this cluster.
  # CONCLUSION : Atrichoblast -> AT2G33790

  # Subset1 Cluster 2
clusters.markers[which(clusters.markers[,6] == 2),][1:5,]
specificity(gala_subset1, "AT1G27140")
FeaturePlot_3D(object = gala_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G27140")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Epidermis to Trichoblast
  # The calculation of the specificity gives a score of 67.9% for this cluster.
  # CONCLUSION : Epidermis to Trichoblast -> AT1G27140

  # Subset1 Cluster 3
clusters.markers[which(clusters.markers[,6] == 3),][1:5,]
specificity(gala_subset1, "AT4G25630")
FeaturePlot_3D(object = gala_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G25630")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: epidermis initials
  # The calculation of the specificity gives a score of 82.7% for this cluster.
  # CONCLUSION : epidermis initials -> AT4G25630

  # Subset1 Cluster 5
clusters.markers[which(clusters.markers[,6] == 5),][1:5,]
specificity(gala_subset1, "AT3G62680") 
FeaturePlot_3D(object = gala_subset1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G62680")
  # De Rybel's atlas indicates that this gene is expressed in the cell type: atrichoblast
  # The calculation of the specificity gives a score of 98.6% for this cluster.
  # CONCLUSION : Trichoblast -> AT3G62680


### Subclustering of Subset1 (clusters 0 and 4)

cells.subset1.1 = c(names(gala_subset1@active.ident[which(gala_subset1@active.ident == 0)]), 
                    names(gala_subset1@active.ident[which(gala_subset1@active.ident == 4)]))
gala_cells.subset1.1 = gala_cells[cells.subset1.1,]
gala_matrix.subset1.1 = gala_matrix[,cells.subset1.1]


gala_subset1.1 = CreateSeuratObject(counts = gala_matrix.subset1.1,
                                  project = "Gala Annotation",
                                  assay = "RNA",
                                  names.field = 1,
                                  names.delim = NULL,
                                  meta.data = gala_cells.subset1.1,
                                  min.cells = 0,
                                  min.features = 0,
                                  row.names = NULL)

# Checking outliers
gala_subset1.1[["percent.mt"]] <- PercentageFeatureSet(gala_subset1.1, pattern = "^ATM")
VlnPlot(gala_subset1.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset1.1 = NormalizeData(gala_subset1.1,
                               normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                               scale.factor = 10000,
                               margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                               block.size = NULL,
                               verbose = TRUE)

gala_subset1.1 = FindVariableFeatures(gala_subset1.1,
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

gala_subset1.1 = ScaleData(gala_subset1.1,
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

gala_subset1.1 = RunPCA(gala_subset1.1,
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

gala_subset1.1 = RunUMAP(gala_subset1.1,
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

gala_subset1.1 = FindNeighbors(gala_subset1.1,
                               reduction = "umap", 
                               dims=1:4,
                               assay = NULL,
                               features =NULL,
                               k.param = 20,
                               return.neighbor = F)

gala_subset1.1 = FindClusters(gala_subset1.1,
                              graph.name = NULL,
                              modularity.fxn = 1,
                              initial.membership = NULL,
                              node.sizes = NULL,
                              resolution = 0.3,
                              method = "matrix",
                              algorithm = 3,
                              n.start = 20,
                              n.iter = 20,
                              random.seed = 0,
                              group.singletons = TRUE,
                              temp.file.location = NULL,
                              edge.file.name = NULL,
                              verbose = TRUE)

DimPlot_3D(object = gala_subset1.1, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# Search for marker gene expression clusters

clusters.markers1.1 = FindAllMarkers(gala_subset1.1, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset1.1 Cluster 0
clusters.markers1.1[which(clusters.markers1.1[,6] == 0),][1:5,]
specificity(gala_subset1.1, "AT1G14870")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G14870")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Lateral Root Cap
  # The calculation of the specificity gives a score of 18.9% for this cluster.
  # CONCLUSION : Lateral Root Cap -> AT1G14870

  # Subset1.1 Cluster 1
clusters.markers1.1[which(clusters.markers1.1[,6] == 1),][1:5,]
specificity(gala_subset1.1, "AT1G28290")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G28290")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Lateral Root Cap
  # The calculation of the specificity gives a score of 26.1% for this cluster.
  # CONCLUSION : Lateral Root Cap -> AT1G28290

  # Subset1.1 Cluster 2
clusters.markers1.1[which(clusters.markers1.1[,6] == 2),][1:5,]
specificity(gala_subset1.1, "AT5G02960")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G02960")
  # TAIR indicates that this gene is very specific for the cell type: Young Lateral Root Cap
  # The calculation of the specificity gives a score of 18.1% for this cluster.
  # CONCLUSION : Young Lateral Root Cap -> AT5G02960

  # Subset1.1 Cluster 3
clusters.markers1.1[which(clusters.markers1.1[,6] == 3),][1:5,]
specificity(gala_subset1.1, "AT3G16440")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G16440")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Mature Lateral Root Cap
  # The calculation of the specificity gives a score of 74.2% for this cluster.
  # CONCLUSION : Mature Lateral Root Cap -> AT3G16440

  # Subset1.1 Cluster 4
clusters.markers1.1[which(clusters.markers1.1[,6] == 4),][1:5,]
specificity(gala_subset1.1, "AT1G02730")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G02730")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Dividing Lateral Root Cap
  # The calculation of the specificity gives a score of 72.3% for this cluster.
  # CONCLUSION : Dividing Lateral Root Cap -> AT1G02730

  # Subset1.1 Cluster 5
clusters.markers1.1[which(clusters.markers1.1[,6] == 5),][1:5,]
specificity(gala_subset1.1, "AT5G59870")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G59870")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Dividing Lateral Root Cap
  # The calculation of the specificity gives a score of 55.4% for this cluster.
  # CONCLUSION : Dividing Lateral Root Cap -> AT5G59870 Histone.HA2AW

  # Subset1.1 Cluster 6
clusters.markers1.1[which(clusters.markers1.1[,6] == 6),][1:5,]
specificity(gala_subset1.1, "AT2G04025")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G04025")
DimPlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", group.by = "hour")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Columella
  # The calculation of the specificity gives a score of 91.7% for this cluster.
  # CONCLUSION : Columella -> AT2G04025

  # Subset1.1 Cluster 7
clusters.markers1.1[which(clusters.markers1.1[,6] == 7),][1:5,]
specificity(gala_subset1.1, "AT3G43430")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G43430")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: 
  # The calculation of the specificity gives a score of 91.7% for this cluster.
  # CONCLUSION : Ambiguous -> AT3G43430

  # Subset1.1 Cluster 8
clusters.markers1.1[which(clusters.markers1.1[,6] == 8),][1:5,]
specificity(gala_subset1.1, "AT4G25200")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G25200")
  # TAIR indicates that this gene is very specific for the cell type: Stressed Lateral Root Cap
  # The calculation of the specificity gives a score of 98.9% for this cluster.
  # CONCLUSION : Stressed Lateral Root Cap -> AT4G25200

  # Subset1.1 Cluster 9
clusters.markers1.1[which(clusters.markers1.1[,6] == 9),][1:5,]
specificity(gala_subset1.1, "AT1G06135")
FeaturePlot_3D(object = gala_subset1.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G06135")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Lateral Root Cap
  # The calculation of the specificity gives a score of 95.4% for this cluster.
  # CONCLUSION : Lateral Root Cap -> AT1G06135



### Subclustering of the entire dataset (cluster 0 and 2)

cells.subset2 = c(names(gala_all@active.ident[which(gala_all@active.ident ==0)]),
                  names(gala_all@active.ident[which(gala_all@active.ident ==2)]))
gala_cells.subset2 = gala_cells[cells.subset2,]
gala_matrix.subset2 = gala_matrix[,cells.subset2]

gala_subset2 = CreateSeuratObject(counts = gala_matrix.subset2,
                                  project = "Gala Annotation",
                                  assay = "RNA",
                                  names.field = 1,
                                  names.delim = NULL,
                                  meta.data = gala_cells.subset2,
                                  min.cells = 0,
                                  min.features = 0,
                                  row.names = NULL)

# Checking outliers
gala_subset2[["percent.mt"]] <- PercentageFeatureSet(gala_subset2, pattern = "^ATM")
VlnPlot(gala_subset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset2 = NormalizeData(gala_subset2,
                             normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                             scale.factor = 10000,
                             margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                             block.size = NULL,
                             verbose = TRUE)

gala_subset2 = FindVariableFeatures(gala_subset2,
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

gala_subset2 = ScaleData(gala_subset2,
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

gala_subset2 = RunPCA(gala_subset2,
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

gala_subset2 = RunUMAP(gala_subset2,
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

gala_subset2 = FindNeighbors(gala_subset2,
                             reduction = "umap", 
                             dims=1:4,
                             assay = NULL,
                             features =NULL,
                             k.param = 20,
                             return.neighbor = F)

gala_subset2 = FindClusters(gala_subset2,
                            graph.name = NULL,
                            modularity.fxn = 1,
                            initial.membership = NULL,
                            node.sizes = NULL,
                            resolution = 0.05,
                            method = "matrix",
                            algorithm = 3,
                            n.start = 20,
                            n.iter = 20,
                            random.seed = 0,
                            group.singletons = TRUE,
                            temp.file.location = NULL,
                            edge.file.name = NULL,
                            verbose = TRUE)

DimPlot_3D(object = gala_subset2, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

  # We get 9 clusters. A more in-depth study is carried out by making sub-clusters.


### Subclustering of Subset2 (clusters 6 and 8)

cells.subset2.1 = c(names(gala_subset2@active.ident[which(gala_subset2@active.ident == 6)]), 
                    names(gala_subset2@active.ident[which(gala_subset2@active.ident == 8)]))
gala_cells.subset2.1 = gala_cells[cells.subset2.1,]
gala_matrix.subset2.1 = gala_matrix[,cells.subset2.1]


gala_subset2.1 = CreateSeuratObject(counts = gala_matrix.subset2.1,
                                    project = "Gala Annotation",
                                    assay = "RNA",
                                    names.field = 1,
                                    names.delim = NULL,
                                    meta.data = gala_cells.subset2.1,
                                    min.cells = 0,
                                    min.features = 0,
                                    row.names = NULL)

# Checking outliers
gala_subset2.1[["percent.mt"]] <- PercentageFeatureSet(gala_subset2.1, pattern = "^ATM")
VlnPlot(gala_subset2.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset2.1 = NormalizeData(gala_subset2.1,
                               normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                               scale.factor = 10000,
                               margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                               block.size = NULL,
                               verbose = TRUE)

gala_subset2.1 = FindVariableFeatures(gala_subset2.1,
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


gala_subset2.1 = ScaleData(gala_subset2.1,
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

gala_subset2.1 = RunPCA(gala_subset2.1,
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

gala_subset2.1 = RunUMAP(gala_subset2.1,
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

gala_subset2.1 = FindNeighbors(gala_subset2.1,
                               reduction = "umap", 
                               dims=1:4,
                               assay = NULL,
                               features =NULL,
                               k.param = 20,
                               return.neighbor = F)

gala_subset2.1 = FindClusters(gala_subset2.1,
                              graph.name = NULL,
                              modularity.fxn = 1,
                              initial.membership = NULL,
                              node.sizes = NULL,
                              resolution = 0.1,
                              method = "matrix",
                              algorithm = 3,
                              n.start = 20,
                              n.iter = 20,
                              random.seed = 0,
                              group.singletons = TRUE,
                              temp.file.location = NULL,
                              edge.file.name = NULL,
                              verbose = TRUE)

DimPlot_3D(object = gala_subset2.1, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# Search for marker gene expression clusters

clusters.markers2.1 = FindAllMarkers(gala_subset2.1, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset2.1 Cluster 0
clusters.markers2.1[which(clusters.markers2.1[,6] == 0),][1:5,]
specificity(gala_subset2.1, "AT5G15970")
FeaturePlot_3D(object = gala_subset2.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G15970")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Intermediate Cortex
  # The calculation of the specificity gives a score of 62.3% for this cluster.
  # CONCLUSION : Intermediate Cortex -> AT5G15970

  # Subset2.1 Cluster 1
clusters.markers2.1[which(clusters.markers2.1[,6] == 1),][1:5,]
specificity(gala_subset2.1, "AT5G15230")
FeaturePlot_3D(object = gala_subset2.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G15230")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Young Cortex Endodermis
  # The calculation of the specificity gives a score of 98.1% for this cluster. (Initials AT1G73620 are 100% specific)
  # CONCLUSION : Young Cortex Endodermis -> AT5G15230

  # Subset2.1 Cluster 2
clusters.markers2.1[which(clusters.markers2.1[,6] == 2),][1:5,]
specificity(gala_subset2.1, "AT5G13080")
FeaturePlot_3D(object = gala_subset2.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G13080")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Differenciated Cortex
  # The calculation of the specificity gives a score of 86.5% for this cluster.
  # CONCLUSION : Differenciated Cortex-> AT5G13080

  # Subset2.1 Cluster 3
clusters.markers2.1[which(clusters.markers2.1[,6] == 3),][1:5,]
specificity(gala_subset2.1, "AT4G03280")
FeaturePlot_3D(object = gala_subset2.1, dims = c(1,2,3), reduction.method = "umap", feature = "")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Differenciated Cortex
  # The calculation of the specificity gives a score of 95.2% for this cluster.
  # CONCLUSION : Differenciated Cortex-> AT4G03280



### Subclustering of Subset2 (cluster 3)

cells.subset2.2 = names(gala_subset2@active.ident[which(gala_subset2@active.ident == 3)])
gala_cells.subset2.2 = gala_cells[cells.subset2.2,]
gala_matrix.subset2.2 = gala_matrix[,cells.subset2.2]


gala_subset2.2 = CreateSeuratObject(counts = gala_matrix.subset2.2,
                                    project = "Gala Annotation",
                                    assay = "RNA",
                                    names.field = 1,
                                    names.delim = NULL,
                                    meta.data = gala_cells.subset2.2,
                                    min.cells = 0,
                                    min.features = 0,
                                    row.names = NULL)

# Checking outliers
gala_subset2.2[["percent.mt"]] <- PercentageFeatureSet(gala_subset2.2, pattern = "^ATM")
VlnPlot(gala_subset2.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset2.2 = NormalizeData(gala_subset2.2,
                               normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                               scale.factor = 10000,
                               margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                               block.size = NULL,
                               verbose = TRUE)

gala_subset2.2 = FindVariableFeatures(gala_subset2.2,
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


gala_subset2.2 = ScaleData(gala_subset2.2,
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

gala_subset2.2 = RunPCA(gala_subset2.2,
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

gala_subset2.2 = RunUMAP(gala_subset2.2,
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

gala_subset2.2 = FindNeighbors(gala_subset2.2,
                               reduction = "umap", 
                               dims=1:4,
                               assay = NULL,
                               features =NULL,
                               k.param = 20,
                               return.neighbor = F)

gala_subset2.2 = FindClusters(gala_subset2.2,
                              graph.name = NULL,
                              modularity.fxn = 1,
                              initial.membership = NULL,
                              node.sizes = NULL,
                              resolution = 0.2,
                              method = "matrix",
                              algorithm = 3,
                              n.start = 20,
                              n.iter = 20,
                              random.seed = 0,
                              group.singletons = TRUE,
                              temp.file.location = NULL,
                              edge.file.name = NULL,
                              verbose = TRUE)

DimPlot_3D(object = gala_subset2.2, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# Search for marker gene expression clusters

clusters.markers2.2 = FindAllMarkers(gala_subset2.2, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset2.2 Cluster 0
clusters.markers2.2[which(clusters.markers2.2[,6] == 0),][1:5,]
specificity(gala_subset2.2, "AT2G37750")
FeaturePlot_3D(object = gala_subset2.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT2G37750")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Differenciated Endodermis
  # The calculation of the specificity gives a score of 93.5% for this cluster.
  # CONCLUSION : Differenciated Endodermis -> AT2G37750

  # Subset2.2 Cluster 1
clusters.markers2.2[which(clusters.markers2.2[,6] == 1),][1:5,]
specificity(gala_subset2.2, "AT1G54580")
FeaturePlot_3D(object = gala_subset2.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G15230")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Young Endodermis
  # The calculation of the specificity gives a score of 90.2% for this cluster.
  # CONCLUSION : Young Endodermis -> AT5G15230

  # Subset2.2 Cluster 2
clusters.markers2.2[which(clusters.markers2.2[,6] == 2),][1:5,]
specificity(gala_subset2.2, "AT3G22600")
FeaturePlot_3D(object = gala_subset2.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G22600")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: 
  # The calculation of the specificity gives a score of 59.7% for this cluster.
  # CONCLUSION : Intermediate Endodermis -> AT3G22600


### Subclustering of Subset2 (clusters 0, 1, 2, 4, 5 and 7)

cells.subset2.3 = c(names(gala_subset2@active.ident[which(gala_subset2@active.ident == 0)]), 
                    names(gala_subset2@active.ident[which(gala_subset2@active.ident == 1)]),
                    names(gala_subset2@active.ident[which(gala_subset2@active.ident == 2)]),
                    names(gala_subset2@active.ident[which(gala_subset2@active.ident == 4)]),
                    names(gala_subset2@active.ident[which(gala_subset2@active.ident == 5)]),
                    names(gala_subset2@active.ident[which(gala_subset2@active.ident == 7)]))
gala_cells.subset2.3 = gala_cells[cells.subset2.3,]
gala_matrix.subset2.3 = gala_matrix[,cells.subset2.3]


gala_subset2.3 = CreateSeuratObject(counts = gala_matrix.subset2.3,
                                    project = "Gala Annotation",
                                    assay = "RNA",
                                    names.field = 1,
                                    names.delim = NULL,
                                    meta.data = gala_cells.subset2.3,
                                    min.cells = 0,
                                    min.features = 0,
                                    row.names = NULL)

# Checking outliers
gala_subset2.3[["percent.mt"]] <- PercentageFeatureSet(gala_subset2.3, pattern = "^ATM")
VlnPlot(gala_subset2.3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset2.3 = NormalizeData(gala_subset2.3,
                               normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                               scale.factor = 10000,
                               margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                               block.size = NULL,
                               verbose = TRUE)

gala_subset2.3 = FindVariableFeatures(gala_subset2.3,
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

gala_subset2.3 = ScaleData(gala_subset2.3,
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

gala_subset2.3 = RunPCA(gala_subset2.3,
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

gala_subset2.3 = RunUMAP(gala_subset2.3,
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

gala_subset2.3 = FindNeighbors(gala_subset2.3,
                               reduction = "umap", 
                               dims=1:4,
                               assay = NULL,
                               features =NULL,
                               k.param = 20,
                               return.neighbor = F)

gala_subset2.3 = FindClusters(gala_subset2.3,
                              graph.name = NULL,
                              modularity.fxn = 1,
                              initial.membership = NULL,
                              node.sizes = NULL,
                              resolution = 0.2,
                              method = "matrix",
                              algorithm = 3,
                              n.start = 20,
                              n.iter = 20,
                              random.seed = 0,
                              group.singletons = TRUE,
                              temp.file.location = NULL,
                              edge.file.name = NULL,
                              verbose = TRUE)

DimPlot_3D(object = gala_subset2.3, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

  # We get 13 clusters. we are therefore going to study sub-clusterizations.


### Subclustering of Subset2.3 (cluster 6)

cells.subset2.3.1 = names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 6)])
gala_cells.subset2.3.1 = gala_cells[cells.subset2.3.1,]
gala_matrix.subset2.3.1 = gala_matrix[,cells.subset2.3.1]

gala_subset2.3.1 = CreateSeuratObject(counts = gala_matrix.subset2.3.1,
                                    project = "Gala Annotation",
                                    assay = "RNA",
                                    names.field = 1,
                                    names.delim = NULL,
                                    meta.data = gala_cells.subset2.3.1,
                                    min.cells = 0,
                                    min.features = 0,
                                    row.names = NULL)

# Checking outliers
gala_subset2.3.1[["percent.mt"]] <- PercentageFeatureSet(gala_subset2.3.1, pattern = "^ATM")
VlnPlot(gala_subset2.3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset2.3.1 = NormalizeData(gala_subset2.3.1,
                               normalization.method = "LogNormalize",    # Can be : "LogNormalize", "CLR" (centered log ratio), "RC" (relative counts)
                               scale.factor = 10000,
                               margin = 1,                               # If performing CLR normalization, normalize across features (1) or cells (2)
                               block.size = NULL,
                               verbose = TRUE)

gala_subset2.3.1 = FindVariableFeatures(gala_subset2.3.1,
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

gala_subset2.3.1 = ScaleData(gala_subset2.3.1,
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

gala_subset2.3.1 = RunPCA(gala_subset2.3.1,
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

gala_subset2.3.1 = RunUMAP(gala_subset2.3.1,
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

gala_subset2.3.1 = FindNeighbors(gala_subset2.3.1,
                               reduction = "umap", 
                               dims=1:4,
                               assay = NULL,
                               features =NULL,
                               k.param = 20,
                               return.neighbor = F)

gala_subset2.3.1 = FindClusters(gala_subset2.3.1,
                              graph.name = NULL,
                              modularity.fxn = 1,
                              initial.membership = NULL,
                              node.sizes = NULL,
                              resolution = 0.2,
                              method = "matrix",
                              algorithm = 3,
                              n.start = 20,
                              n.iter = 20,
                              random.seed = 0,
                              group.singletons = TRUE,
                              temp.file.location = NULL,
                              edge.file.name = NULL,
                              verbose = TRUE)

DimPlot_3D(object = gala_subset2.3.1, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# Search for marker gene expression clusters

clusters.markers2.3.1 = FindAllMarkers(gala_subset2.3.1, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset2.3.1 Cluster 0
clusters.markers2.3.1[which(clusters.markers2.3.1[,6] == 0),][1:5,]
specificity(gala_subset2.3.1, "AT1G62480")
FeaturePlot_3D(object = gala_subset2.3.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G62480")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Phloem
  # The calculation of the specificity gives a score of 97.6% for this cluster.
  # CONCLUSION : Phloem -> AT1G62480

  # Subset2.3.1 Cluster 1
clusters.markers2.3.1[which(clusters.markers2.3.1[,6] == 1),][1:5,]
specificity(gala_subset2.3.1, "AT4G27520")
FeaturePlot_3D(object = gala_subset2.3.1, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G27520")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Companion Cell
  # The calculation of the specificity gives a score of 89.8% for this cluster.
  # CONCLUSION : Companion Cell -> AT4G27520



### Subclustering of Subset2.3 (clusters 1, 2, 3, 4, 5, 7, 11, 12 )

cells.subset2.3.2 = c(names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 1)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 2)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 3)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 4)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 5)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 7)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 11)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 12)]))
gala_cells.subset2.3.2 = gala_cells[cells.subset2.3.2,]
gala_matrix.subset2.3.2 = gala_matrix[,cells.subset2.3.2]


gala_subset2.3.2 = CreateSeuratObject(counts = gala_matrix.subset2.3.2,
                                      project = "Gala Annotation",
                                      assay = "RNA",
                                      names.field = 1,
                                      names.delim = NULL,
                                      meta.data = gala_cells.subset2.3.2,
                                      min.cells = 0,
                                      min.features = 0,
                                      row.names = NULL)

# Checking outliers
gala_subset2.3.2[["percent.mt"]] <- PercentageFeatureSet(gala_subset2.3.2, pattern = "^ATM")
VlnPlot(gala_subset2.3.2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # No filtering was finally done.

gala_subset2.3.2 = NormalizeData(gala_subset2.3.2,
                                 normalization.method = "LogNormalize",    
                                 scale.factor = 10000,
                                 margin = 1,                               
                                 block.size = NULL,
                                 verbose = TRUE)

gala_subset2.3.2 = FindVariableFeatures(gala_subset2.3.2,
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

gala_subset2.3.2 = ScaleData(gala_subset2.3.2,
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

gala_subset2.3.2 = RunPCA(gala_subset2.3.2,
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

gala_subset2.3.2 = RunUMAP(gala_subset2.3.2,
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

gala_subset2.3.2 = FindNeighbors(gala_subset2.3.2,
                                 reduction = "umap", 
                                 dims=1:4,
                                 assay = NULL,
                                 features =NULL,
                                 k.param = 20,
                                 return.neighbor = F)

gala_subset2.3.2 = FindClusters(gala_subset2.3.2,
                                graph.name = NULL,
                                modularity.fxn = 1,
                                initial.membership = NULL,
                                node.sizes = NULL,
                                resolution = 0.1,
                                method = "matrix",
                                algorithm = 3,
                                n.start = 20,
                                n.iter = 20,
                                random.seed = 0,
                                group.singletons = TRUE,
                                temp.file.location = NULL,
                                edge.file.name = NULL,
                                verbose = TRUE)

DimPlot_3D(object = gala_subset2.3.2, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")

# Search for marker gene expression clusters

clusters.markers2.3.2 = FindAllMarkers(gala_subset2.3.2, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset2.3.2 Cluster 0
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 0),][1:5,]
specificity(gala_subset2.3.2, "AT1G05340")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G05340")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 31.1% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT1G05340

  # Subset2.3.2 Cluster 1
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 1),][1:5,]
specificity(gala_subset2.3.2, "AT4G26320")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G26320")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Xylem Pole Pericycle
  # The calculation of the specificity gives a score of 54.1% for this cluster.
  # CONCLUSION : Xylem Pole Pericycle -> AT4G26320

# Subset2.3.2 Cluster 2
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 2),][1:5,]
specificity(gala_subset2.3.2, "AT4G36430")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G36430")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 74.5% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT4G36430

  # Subset2.3.2 Cluster 3
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 3),][1:5,]
specificity(gala_subset2.3.2, "AT1G48750")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT1G48750")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Mature Pericycle
  # The calculation of the specificity gives a score of 70.0% for this cluster.
  # CONCLUSION : Mature Pericycle -> AT1G48750

  # Subset2.3.2 Cluster 4
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 4),][1:5,]
specificity(gala_subset2.3.2, "AT3G23830")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G23830")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: LRP + dividing + QC
  # The calculation of the specificity gives a score of 74.5% for this cluster.
  # CONCLUSION : LRP -> AT3G23830

  # Subset2.3.2 Cluster 5
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 5),][1:5,]
specificity(gala_subset2.3.2, "AT3G59370")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G59370")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 75.9% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT3G59370

  # Subset2.3.2 Cluster 6
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 6),][1:5,]
specificity(gala_subset2.3.2, "AT3G46230")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G46230")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Phloem Pole Pericycle
  # The calculation of the specificity gives a score of 92.5% for this cluster.
  # CONCLUSION : Phloem Pole Pericycle -> AT3G46230

  # Subset2.3.2 Cluster 7
clusters.markers2.3.2[which(clusters.markers2.3.2[,6] == 7),][1:5,]
specificity(gala_subset2.3.2, "AT5G60530")
FeaturePlot_3D(object = gala_subset2.3.2, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G60530")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Ambiguous
  # The calculation of the specificity gives a score of 94.4% for this cluster.
  # CONCLUSION : Ambiguous -> AT5G60530


### Subclustering of Subset2.3 (clusters 0, 8, 9, 10 )

cells.subset2.3.3 = c(names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 0)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 8)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 9)]),
                      names(gala_subset2.3@active.ident[which(gala_subset2.3@active.ident == 10)]))
gala_cells.subset2.3.3 = gala_cells[cells.subset2.3.3,]
gala_matrix.subset2.3.3 = gala_matrix[,cells.subset2.3.3]


gala_subset2.3.3 = CreateSeuratObject(counts = gala_matrix.subset2.3.3,
                                      project = "Gala Annotation",
                                      assay = "RNA",
                                      names.field = 1,
                                      names.delim = NULL,
                                      meta.data = gala_cells.subset2.3.3,
                                      min.cells = 0,
                                      min.features = 0,
                                      row.names = NULL)

# Checking outliers
gala_subset2.3.3[["percent.mt"]] <- PercentageFeatureSet(gala_subset2.3.3, pattern = "^ATM")
VlnPlot(gala_subset2.3.3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# No filtering was finally done.

gala_subset2.3.3 = NormalizeData(gala_subset2.3.3,
                                 normalization.method = "LogNormalize",    
                                 scale.factor = 10000,
                                 margin = 1,                               
                                 block.size = NULL,
                                 verbose = TRUE)

gala_subset2.3.3 = FindVariableFeatures(gala_subset2.3.3,
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

gala_subset2.3.3 = ScaleData(gala_subset2.3.3,
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

gala_subset2.3.3 = RunPCA(gala_subset2.3.3,
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

gala_subset2.3.3 = RunUMAP(gala_subset2.3.3,
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

gala_subset2.3.3 = FindNeighbors(gala_subset2.3.3,
                                 reduction = "umap", 
                                 dims=1:4,
                                 assay = NULL,
                                 features =NULL,
                                 k.param = 20,
                                 return.neighbor = F)

gala_subset2.3.3 = FindClusters(gala_subset2.3.3,
                                graph.name = NULL,
                                modularity.fxn = 1,
                                initial.membership = NULL,
                                node.sizes = NULL,
                                resolution = 0.1,
                                method = "matrix",
                                algorithm = 3,
                                n.start = 20,
                                n.iter = 20,
                                random.seed = 0,
                                group.singletons = TRUE,
                                temp.file.location = NULL,
                                edge.file.name = NULL,
                                verbose = TRUE)

DimPlot_3D(object = gala_subset2.3.3, dims=c(1:3), reduction.method = "umap", group.by = "seurat_clusters")


# Search for marker gene expression clusters

clusters.markers2.3.3 = FindAllMarkers(gala_subset2.3.3, only.pos =T, min.pct = 0.05, logfc.threshold = 0.25)

  # Subset2.3.3 Cluster 0
clusters.markers2.3.3[which(clusters.markers2.3.3[,6] == 0),][1:5,]
specificity(gala_subset2.3.3, "AT3G16330")
FeaturePlot_3D(object = gala_subset2.3.3, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G16330")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Ambiguous Stele Cell
  # The calculation of the specificity gives a score of 47.6% for this cluster.
  # CONCLUSION : Ambiguous Stele Cell -> AT3G16330

  # Subset2.3.3 Cluster 1
clusters.markers2.3.3[which(clusters.markers2.3.3[,6] == 1),][1:5,]
specificity(gala_subset2.3.3, "AT4G35100")
FeaturePlot_3D(object = gala_subset2.3.3, dims = c(1,2,3), reduction.method = "umap", feature = "AT4G35100")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Metaxylem
  # The calculation of the specificity gives a score of 70.1% for this cluster.
  # CONCLUSION : Ambiguous Stele Cells -> AT4G35100


  # Subset2.3.3 Cluster 2
clusters.markers2.3.3[which(clusters.markers2.3.3[,6] == 2),][1:5,]
specificity(gala_subset2.3.3, "AT3G18280")
FeaturePlot_3D(object = gala_subset2.3.3, dims = c(1,2,3), reduction.method = "umap", feature = "AT3G18280")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Metaxylem
  # The calculation of the specificity gives a score of 40.2% for this cluster.
  # CONCLUSION :  Metaxylem-> AT3G18280

  # Subset2.3.3 Cluster 3
clusters.markers2.3.3[which(clusters.markers2.3.3[,6] == 3),][1:5,]
specificity(gala_subset2.3.3, "AT5G12870")
FeaturePlot_3D(object = gala_subset2.3.3, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G12870")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Protoxylem
  # The calculation of the specificity gives a score of 99.6% for this cluster.
  # CONCLUSION : Protoxylem -> AT5G12870

  # Subset2.3.3 Cluster 4
clusters.markers2.3.3[which(clusters.markers2.3.3[,6] == 4),][1:5,]
specificity(gala_subset2.3.3, "AT5G59870")
FeaturePlot_3D(object = gala_subset2.3.3, dims = c(1,2,3), reduction.method = "umap", feature = "AT5G59870")
  # De Rybel's atlas indicates that this gene is very specific for the cell type: Ambiguous Stele Cell
  # The calculation of the specificity gives a score of 86.3% for this cluster.
  # CONCLUSION : Ambiguous Stele Cell -> AT5G59870



### Cell Annotation

## Subset 1 
cells.cluster.subset1 = gala_subset1@active.ident
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 1)]),4] = "Atrichoblast"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 1)]),5] = "AT2G33790"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 2)]),4] = "Epidermis_to_Trichoblast"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 2)]),5] = "AT1G27140"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 3)]),4] = "Epidermis_Initials"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 3)]),5] = "AT4G25630"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 5)]),4] = "Trichoblast"
gala_cells[names(cells.cluster.subset1[which(cells.cluster.subset1 == 5)]),5] = "AT3G62680"

## Subset 1.1
cells.cluster.subset1.1 = gala_subset1.1@active.ident
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 0)]),4] = "Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 0)]),5] = "AT1G14870"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 1)]),4] = "Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 1)]),5] = "AT1G28290"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 2)]),4] = "Young_Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 2)]),5] = "AT5G02960"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 3)]),4] = "Mature_Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 3)]),5] = "AT3G16440"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 4)]),4] = "Dividing_Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 4)]),5] = "AT1G02730"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 5)]),4] = "Dividing_Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 5)]),5] = "AT5G59870"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 6)]),4] = "Columella"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 6)]),5] = "AT2G04025"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 7)]),4] = "Ambiguous"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 7)]),5] = "AT3G43430"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 8)]),4] = "Stressed_Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 8)]),5] = "AT4G25200"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 9)]),4] = "Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset1.1[which(cells.cluster.subset1.1 == 9)]),5] = "AT1G06135"


## Subset 2.1
cells.cluster.subset2.1 = gala_subset2.1@active.ident
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 0)]),4] = "Intermediate_Cortex"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 0)]),5] = "AT5G15970"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 1)]),4] = "Young_Cortex"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 1)]),5] = "AT5G15230"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 2)]),4] = "Differenciated_Cortex"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 2)]),5] = "AT5G13080"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 3)]),4] = "Differenciated_Cortex"
gala_cells[names(cells.cluster.subset2.1[which(cells.cluster.subset2.1 == 3)]),5] = "AT4G03280"


## Subset 2.2
cells.cluster.subset2.2 = gala_subset2.2@active.ident
gala_cells[names(cells.cluster.subset2.2[which(cells.cluster.subset2.2 == 0)]),4] = "Differenciated_Endodermis"
gala_cells[names(cells.cluster.subset2.2[which(cells.cluster.subset2.2 == 0)]),5] = "AT2G37750"
gala_cells[names(cells.cluster.subset2.2[which(cells.cluster.subset2.2 == 1)]),4] = "Young_Endodermis"
gala_cells[names(cells.cluster.subset2.2[which(cells.cluster.subset2.2 == 1)]),5] = "AT3G15230"
gala_cells[names(cells.cluster.subset2.2[which(cells.cluster.subset2.2 == 2)]),4] = "Intermediate_Endodermis"
gala_cells[names(cells.cluster.subset2.2[which(cells.cluster.subset2.2 == 2)]),5] = "AT3G22600"


## Subset 2.3.1
cells.cluster.subset2.3.1 = gala_subset2.3.1@active.ident
gala_cells[names(cells.cluster.subset2.3.1[which(cells.cluster.subset2.3.1 == 0)]),4] = "Phloem"
gala_cells[names(cells.cluster.subset2.3.1[which(cells.cluster.subset2.3.1 == 0)]),5] = "AT1G62480"
gala_cells[names(cells.cluster.subset2.3.1[which(cells.cluster.subset2.3.1 == 1)]),4] = "Companion_Cell"
gala_cells[names(cells.cluster.subset2.3.1[which(cells.cluster.subset2.3.1 == 1)]),5] = "AT4G27520"


## Subset 2.3.2
cells.cluster.subset2.3.2 = gala_subset2.3.2@active.ident
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 0)]),4] = "Phloem_Pole_Pericycle"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 0)]),5] = "AT1G05340"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 1)]),4] = "Xylem_Pole_Pericycle"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 1)]),5] = "AT4G26320"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 2)]),4] = "Phloem_Pole_Pericycle"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 2)]),5] = "AT4G36430"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 3)]),4] = "Mature_Pericycle"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 3)]),5] = "AT1G48730"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 4)]),4] = "Lateral_Root_Primordia"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 4)]),5] = "AT3G23830"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 5)]),4] = "Phloem_Pole_Pericycle"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 5)]),5] = "AT3G59370"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 6)]),4] = "Phloem_Pole_Pericycle"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 6)]),5] = "AT3G46230"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 7)]),4] = "Ambiguous"
gala_cells[names(cells.cluster.subset2.3.2[which(cells.cluster.subset2.3.2 == 7)]),5] = "AT5G60530"


## Subset 2.3.3
cells.cluster.subset2.3.3 = gala_subset2.3.3@active.ident
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 0)]),4] = "Ambiguous_Stele_Cell"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 0)]),5] = "AT3G16330"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 1)]),4] = "Ambiguous_Stele_Cell"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 1)]),5] = "AT4G35100"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 2)]),4] = "Metaxylem"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 2)]),5] = "AT3G18280"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 3)]),4] = "Protoxylem"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 3)]),5] = "AT5G12870"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 4)]),4] = "Young_Lateral_Root_Cap"
gala_cells[names(cells.cluster.subset2.3.3[which(cells.cluster.subset2.3.3 == 4)]),5] = "AT5G59870"


### Storing new annotations
write.table(gala_cells, "GSE158761_barcodes_reannotated.tsv")
  # To view the new cell type assignment, use the group.by = "cell_type_ips2" parameter.
