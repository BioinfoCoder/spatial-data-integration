library(RANN)
library(scater)
library(scuttle)
install.packages("harmony")

library(harmony)
library(ggspavis)
devtools::install_github("lmweber/ggspavis")
BiocManager:: install ("VisiumIO") 
library(VisiumIO)
library(patchwork)
install.packages("cli")
install.packages("RPostgreSQL")
library(RPostgreSQL)
if (!require("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("estellad/OSTA.data")

library(OSTA.data)
# Install devtools, if necessary
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")


devtools::install_github("edward130603/BayesSpace")
library(BayesSpace)
library(Seurat)
library(SpatialExperiment)
BiocManager:: install ("SpatialExperimentIO") 

library(SpatialExperimentIO)
# set seed for random number generation
# in order to make results reproducible
set.seed(194849)
# Visium
# Visium
id <- "Visium_HumanBreast_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
obj <- TENxVisium(
  spacerangerOut=td, 
  images="lowres", 
  format="h5")
(vis <- import(obj))

# Xenium
id <- "Xenium_HumanBreast1_Janesick"
pa <- OSTA.data_load(id)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
(xen <- readXeniumSXE(td))

# also retrieve cell subpopulation labels
df <- read.csv(file.path(td, "annotation.csv"))
xen$anno <- df$Annotation[match(xen$cell_id, df$Barcode)]

# simplify spatial coordinate names
spatialCoordsNames(vis) <- spatialCoordsNames(xen) <- c("x", "y")
# use gene symbols as feature names
rownames(vis) <- make.unique(rowData(vis)$Symbol)
rownames(xen) <- make.unique(rowData(xen)$Symbol)

# affine matrix for aligning Xenium onto Visium
mtx <- matrix(nrow=2, byrow=TRUE, c(
  8.82797498e-02, -1.91831377e+00, 1.63476055e+04,
  1.84141210e+00,  5.96797885e-02, 4.12499099e+03),
  dimnames=list(c("x", "y"), c("x", "y", "z")))

# stash original coordinates
old <- spatialCoords(xen)
colData(xen)[c(".x", ".y")] <- old
# apply affine transformation
new <- old %*% t(mtx[, -3]) # scale/rotate &
new <- sweep(new, 2, mtx[, 3], `+`) # offset
spatialCoords(xen) <- new

df <- data.frame(spatialCoords(vis))
fd <- data.frame(spatialCoords(xen))
ggplot() + coord_equal() + theme_void() +
  geom_point(aes(x, y), df, col="grey", stroke=0, size=1) +
  geom_point(aes(x, y), fd, col="blue", stroke=0, size=0.1) 

# do a fixed-radius search to get cell 
# centroids that fall on a given spot
nns <- nn2(
  searchtype="radius", radius=55/2/0.2125, k=200,
  data=spatialCoords(xen), query=spatialCoords(vis))

# get cell indices and number of cells per spot
vis$n_cells <- rowSums((idx <- nns$nn.idx) > 0)
plotSpots(vis, annotate="n_cells")

# aggregate Xenium data into pseudo-spots
xem <- aggregateAcrossCells(xen[, c(t(idx))], rep.int(seq(ncol(vis)), vis$n_cells))
# propagate Visium data's spatial coordinates, excluding empty pseudo-spots
spatialCoords(xem) <- spatialCoords(vis)[vis$n_cells > 0, ] 
colnames(xem) <- colnames(vis)[vis$n_cells > 0]
xem$in_tissue <- 1
imgData(xem) <- imgData(vis)

sub <- list(gs=intersect(rownames(vis), rownames(xen)))
vis <- addPerCellQCMetrics(vis, subsets=sub)
xem <- addPerCellQCMetrics(xem, subsets=sub)

plotVisium(vis, 
           annotate="subsets_gs_sum",
           zoom=TRUE, facets=NULL) +
  ggtitle("Visium") +
  plotVisium(xem, 
             annotate="subsets_gs_sum",
             zoom=TRUE, facets=NULL) +
  ggtitle("Xenium") +
  plot_layout(nrow=1) 

df <- data.frame(
  n_cells=xem$ncells,
  Xenium=xem$subsets_gs_sum,
  Visium=vis[, colnames(xem)]$subsets_gs_sum)
ggplot(df, aes(Xenium, Visium, col=n_cells)) + 
  scale_color_gradientn(colors=rev(hcl.colors(9, "Mako"))) +
  geom_point(alpha=0.5) + theme_bw() + theme(aspect.ratio=1) 
