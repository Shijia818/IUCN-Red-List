
library(dplyr)
library(raster)
library(doParallel)
library(ape)
library(stringr)
library(apTreeshape)
library(picante)
library(data.table)
library(foreach)
library(rgdal)
library(ggtree)

#################
setwd("D:/Shijia/Threatened_species")

##### future #####
biomod.spdis <- get(load("./distribution_data/current/biomod_spdis.RData"))
####RLI_index <- read.csv("./Zonation/Threatened_Species/Species/RLI_index.csv")
####pos <- which(RLI_index$stable_rcp85!= "LC" & RLI_index$stable_rcp85!= "NT")
####IUCN.future<- RLI_index[pos,]
####biomod.spdis.future <- biomod.spdis[,which(colnames(biomod.spdis) %in% IUCN.future$code)]

RLI_weight <- read.csv("./Zonation/Threatened_Species/Species/Weight/IUCN.total/stable/RCP85.csv")
biomod.spdis.future <- biomod.spdis[,which(colnames(biomod.spdis) %in% RLI_weight$Code)]


############
#########
temp <- as.data.frame(rownames(biomod.spdis.future))
biomod.spdis.future <-cbind(temp,biomod.spdis.future)
colnames(biomod.spdis.future)[1]<-"Gridcode"

#################
basic_shp <-  shapefile('./grid20km/grid20kmChina.shp')
basic_shp$GRIDCODE <- as.character(basic_shp$GRIDCODE)

###Prepare point to extract values from shapefiles
basic_raster <- raster('./raster/China_r.tif')
coord <- coordinates(basic_raster) %>% data.frame()
coord <- coord[!is.na(values(basic_raster)),]
coordinates(coord) <- ~x+y 
crs(coord) <- crs(basic_shp) 

#####
basic_shp2 <- basic_shp["GRIDCODE"]
basic_shp2@data <-left_join(basic_shp2@data, biomod.spdis.future, by = c("GRIDCODE" = 'Gridcode'))
basic_shp2$"GRIDCODE" <- NULL
layers <- names(basic_shp2)
cl <- parallel::makeCluster(30)
registerDoParallel(cl)
foreach(ii = 1:ncol(basic_shp2),
        .packages = c('raster')) %dopar% {
          basic_raster2 <- basic_raster
          values(basic_raster2)[!is.na(values(basic_raster))] <- over(coord, basic_shp2[ii])[, 1]
          raster::writeRaster(basic_raster2,
                              file.path("./Zonation/Threatened_Species/Species/Layers/Total/stable/RCP85/", layers[ii]),
                              format = 'GTiff',
                              overwrite = TRUE)
          NULL
        }
stopCluster(cl)
setwd("D:/Shijia/Threatened_species/Zonation/Threatened_Species/Species/Layers/Total/stable/RCP85")
allfile = dir() 
txtfile <- grep("*.tif.aux.xml", allfile)
file.remove(allfile[txtfile])





###########################################################
########### Phylogeny-based analysis ###################
setwd("D:/SHIJIA/Threatened_species")
source("./Analyses/Zonation/Threatened_Species/Phylogeny/phylo_branch_matrix.R")
names_match <-read.csv("./names_match_10027.csv")
biomod.spdis <- get(load("./distribution_data/current/biomod_spdis.RData"))
colnames(biomod.spdis) <- names_match$name
RLI_index <- read.csv("./Analyses/Zonation/Threatened_Species/Species/RLI_index.csv")
tree<-get(load("./Analyses/Zonation/Threatened_Species/Phylogeny/Random100.RData"))

#################################################
pos <- which(RLI_index$Qin!="LC" & RLI_index$Qin!="NT" & RLI_index$Qin!="DD")
#pos <-  which(RLI_index$stable_rcp26!= "LC" & RLI_index$stable_rcp26!= "NT")
#pos <- which((RLI_index$stable_rcp26!= "LC" & RLI_index$stable_rcp26!= "NT")|(RLI_index$Qin!="LC" & RLI_index$Qin!="NT" & RLI_index$Qin!="DD"))
RLI_index_1 <- RLI_index[pos, ]
biomod.spdis.1 <- biomod.spdis[,pos]


#pos.tree <- colnames(biomod.spdis.1)[which(colnames(biomod.spdis.1) %in% tree[[1]]$tip.label==T)]
#biomod.spdis.final <- biomod.spdis.1[,pos.tree]

#res <- which(tree[[1]]$tip.label %in% colnames(biomod.spdis.final) == F)
#tree.final<-drop.tip(tree[[1]],tree[[1]]$tip.label[res])

######Integration of occurrence matrix and phylogenetic data ########
basic_shp <-  shapefile('./grid20km/grid20kmChina.shp')
basic_shp$GRIDCODE <- as.character(basic_shp$GRIDCODE)

###Prepare point to extract values from shapefiles
basic_raster <- raster('./raster/China_r.tif')
coord <- coordinates(basic_raster) %>% data.frame()
coord <- coord[!is.na(values(basic_raster)),]
coordinates(coord) <- ~x+y 
crs(coord) <- crs(basic_shp) 

#########
dir.create('./Zonation/Threatened_Species/Phylogeny/tree_100')

#### i=1 ###

for(i in 1:length(tree)) {
  message(i)
  pos.tree <- colnames(biomod.spdis.1)[which(colnames(biomod.spdis.1) %in% tree[[i]]$tip.label==T)]
  biomod.spdis.final <- biomod.spdis.1[,pos.tree]
  res <- which(tree[[i]]$tip.label %in% colnames(biomod.spdis.final) == F)
  tree.final<-drop.tip(tree[[i]],tree[[i]]$tip.label[res])
  dirs <- file.path('./Zonation/Threatened_Species/Phylogeny/tree_100/', paste0('Phylotree_', i))
  dir.create(dirs)
  m <-phylo_branch_matrix(site_sp_matrix = biomod.spdis.final, tree_file = tree.final)
  data.table::fwrite(m, file.path(dirs, paste0("phylo_", i, '.gz')))
}

## Transform each matrix to raster 
mlist <- list.dirs('./Zonation/Threatened_Species/Phylogeny/tree_100') %>% grep('Phylotree_',., value = T) %>% list.files(pattern = '.gz$', full.names = TRUE)

for (i in 1:length(mlist)){
  message(i)
  dirs <- file.path('./Zonation/Threatened_Species/Phylogeny/tree_100/', paste0('Phylotree_', i))
  m <- data.table::fread(mlist[i]) %>% data.frame()
  m$ncell <- as.character(m$ncell)
  basic_shp2 <- basic_shp["GRIDCODE"]
  basic_shp2@data <-left_join(basic_shp2@data, m, by = c("GRIDCODE" = 'ncell'))
  basic_shp2$"GRIDCODE" <- NULL
  
  layers <- names(basic_shp2)
  layers <- gsub('X', '', layers)
  cl <- parallel::makeCluster(30)
  registerDoParallel(cl)
  foreach(ii = 1:ncol(basic_shp2),
          .packages = c('raster')) %dopar% {
            basic_raster2 <- basic_raster
            values(basic_raster2)[!is.na(values(basic_raster))] <-
              over(coord, basic_shp2[ii])[, 1]
            raster::writeRaster(basic_raster2,filename=paste("D:/Shijia/Threatened_species/Zonation/Threatened_Species/Phylogeny/Layers/tree_100/","tree_",i,"/",layers[ii],".tif",sep=""))
          }
  
  stopCluster(cl)
}

#### delete useless files
allfile = dir() 
txtfile <- grep("*.tif.aux.xml", allfile)
file.remove(allfile[txtfile])


#### assign each branch a weight ####

############ descendant of node############
IUCN_Qin <- read.csv("./Analyses/Zonation/Threatened_Species/Species/Weight/IUCN.Qin.csv")

for(i in 1:length(tree)) {
message(i)
pos.tree <- colnames(biomod.spdis.1)[which(colnames(biomod.spdis.1) %in% tree[[i]]$tip.label==T)]
biomod.spdis.final <- biomod.spdis.1[,pos.tree]
res <- which(tree[[i]]$tip.label %in% colnames(biomod.spdis.final) == F)
tree.final<-drop.tip(tree[[i]],tree[[i]]$tip.label[res])
site_sp_matrix <- biomod.spdis.final
tree_file <- tree.final
sp_names <- colnames(site_sp_matrix)
site_species_mat <- site_sp_matrix
branches <- cbind(tree_file$edge, data.frame(edge_length=tree_file$edge.length))
ext_branches <-branches[which(branches[, 2] <= length(tree_file$tip.label)), ] ### extract descendant ###
labs <- data.frame(tip.label = tree_file$tip.label)
ext_branches.n <- merge(ext_branches, labs, by.x = '2', by.y = 0) ## by.y =0 represent row #
int_branches <-branches[which(!branches[, 2] %in% 1:length(tree_file$tip.label)),]
dimnames <-list(ext_branches.n$"tip.label", int_branches$"2")
tip_branch_mat <-matrix(0, nrow(ext_branches.n), nrow(int_branches), dimnames = dimnames)
idx_first_branch <- min(int_branches$"2")
for (brch in int_branches$"2") {
  tips.nums <- picante::internal2tips(tree_file, brch)
  tip_branch_mat[tips.nums, brch - idx_first_branch + 1] <- 1
}
#######
newnames <- fortify(tree_file) %>% dplyr::select(node, label) %>% arrange(label) %>% na.omit() %>% data.frame()
tip_branch_mat_1 <- data.frame(tip_branch_mat)
newnames1 <- newnames[match(rownames(tip_branch_mat_1), newnames$label),]
data.ext <- matrix(0,nrow(newnames1),nrow(newnames1))
rownames(data.ext) <- newnames1$label
colnames(data.ext) <- newnames1$node
data.ext <- data.frame(data.ext)
for(j in 1:ncol(data.ext)){
  data.ext[j,j] <- 1
}
all_dist <- cbind(data.ext,tip_branch_mat_1)
#save(all_dist,file="./Zonation/Phylogeny/weight/Phylotree_1/all_dist.RData")

########
#all_dist <- get(load("./Zonation/Phylogeny/weight/Phylotree_1/all_dist.RData"))
IUCN_Qin_1 <- IUCN_Qin[match(rownames(all_dist), IUCN_Qin$Species),]

#IUCN_total <- read.csv("./Zonation/Threatened_Species/Species/Weight/IUCN.total/stable/RCP26.csv")
#IUCN_total_1 <- IUCN_total[match(rownames(all_dist), IUCN_total$Species),]

###
res <- apply(all_dist, 2, sum)
weight_total <- tibble(code = colnames(all_dist), sum = res, Weight = 0)

for(k in 1:nrow(weight_total)){
  index <- all_dist[,k]
  index.weight <- sum(index * IUCN_Qin_1$Weight)  
  weight_total[k,3] <- index.weight/weight_total[k,2]
  write.csv(weight_total, file = paste("./Analyses/Zonation/Threatened_Species/Phylogeny/Weight/","tree_",i,".csv",sep=""))
}
   }






