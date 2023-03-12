library(ggplot2)
library(cowplot)
library(VennDiagram)
library(terra)
library(dplyr)
library(raster)
library(doParallel)
library(adephylo)
library(ape)
library(stringr)
library(apTreeshape)
library(picante)
library(data.table)
library(foreach)
library(rgdal)
library(ggtree)
library(ggsci)
library(RColorBrewer)

#######
setwd("D:/SHIJIA/Threatened_species")

#######
RLI_index <- read.csv("./Analyses/endanger_analysis/RLI_index.csv") ### conservation status based on Qin and Peng ###
dat.Qin <- RLI_index[which(RLI_index$Qin=="EX" | RLI_index$Qin=="CR" | RLI_index$Qin=="EN" | RLI_index$Qin=="VU"),]
sp.Qin <- dat.Qin$Species

dat.full.rcp26 <- RLI_index[which(RLI_index$full_rcp26=="EX" | RLI_index$full_rcp26=="CR" | RLI_index$full_rcp26=="EN" | RLI_index$full_rcp26=="VU"),]
sp.full.rcp26 <- dat.full.rcp26$Species

dat.full.rcp60 <- RLI_index[which(RLI_index$full_rcp60=="EX" | RLI_index$full_rcp60=="CR" | RLI_index$full_rcp60=="EN" | RLI_index$full_rcp60=="VU"),]
sp.full.rcp60 <- dat.full.rcp60$Species

dat.full.rcp85 <- RLI_index[which(RLI_index$full_rcp85=="EX" | RLI_index$full_rcp85=="CR" | RLI_index$full_rcp85=="EN" | RLI_index$full_rcp85=="VU"),]
sp.full.rcp85 <- dat.full.rcp85$Species

dat.buffer.rcp26 <- RLI_index[which(RLI_index$buffer_rcp26=="EX" | RLI_index$buffer_rcp26=="CR" | RLI_index$buffer_rcp26=="EN" | RLI_index$buffer_rcp26=="VU"),]
sp.buffer.rcp26 <- dat.buffer.rcp26$Species

dat.buffer.rcp60 <- RLI_index[which(RLI_index$buffer_rcp60=="EX" | RLI_index$buffer_rcp60=="CR" | RLI_index$buffer_rcp60=="EN" | RLI_index$buffer_rcp60=="VU"),]
sp.buffer.rcp60 <- dat.buffer.rcp60$Species

dat.buffer.rcp85 <- RLI_index[which(RLI_index$buffer_rcp85=="EX" | RLI_index$buffer_rcp85=="CR" | RLI_index$buffer_rcp85=="EN" | RLI_index$buffer_rcp85=="VU"),]
sp.buffer.rcp85 <- dat.buffer.rcp85$Species

dat.stable.rcp26 <- RLI_index[which(RLI_index$stable_rcp26=="EX" | RLI_index$stable_rcp26=="CR" | RLI_index$stable_rcp26=="EN" | RLI_index$stable_rcp26=="VU"),]
sp.stable.rcp26 <- dat.stable.rcp26$Species

dat.stable.rcp60 <- RLI_index[which(RLI_index$stable_rcp60=="EX" | RLI_index$stable_rcp60=="CR" | RLI_index$stable_rcp60=="EN" | RLI_index$stable_rcp60=="VU"),]
sp.stable.rcp60 <- dat.stable.rcp60$Species

dat.stable.rcp85 <- RLI_index[which(RLI_index$stable_rcp85=="EX" | RLI_index$stable_rcp85=="CR" | RLI_index$stable_rcp85=="EN" | RLI_index$stable_rcp85=="VU"),]
sp.stable.rcp85 <- dat.stable.rcp85$Species


#####################################
##### spatial patterns in RLI ##########

id <- c(0,0,1,2,3,4,5)
status <- c("LC","DD","NT","VU","EN","CR","EX")
data.id<-data.frame(id,status)

biomod.spdis <- get(load("./distribution_data/current/biomod_spdis.RData"))
RLI_index <- read.csv("./Analyses/endanger_analysis/RLI_index.csv")
data.RLI <- as.data.frame(matrix(NA,23718,10))
rownames(data.RLI) <- rownames(biomod.spdis)
colnames(data.RLI) <- c("Qin","full_rcp26","full_rcp60","full_rcp85","buffer_rcp26","buffer_rcp60","buffer_rcp85","stable_rcp26","stable_rcp60","stable_rcp85")
for(i in 1:23718){
  res <- which(biomod.spdis[i,]==1)
  spnames <- colnames(biomod.spdis)[res]
  RLI_index_1 <- RLI_index[which(RLI_index$code %in% spnames),]
  
  data.id.Qin <-data.id[match(RLI_index_1$Qin, data.id$status),]
  N <- length(which(data.id.Qin$status!="DD"))
  data.RLI[i,1] <- 1-(sum(data.id.Qin$id)/(5*N))
  
  data.id.full26 <-data.id[match(RLI_index_1$full_rcp26, data.id$status),]
  data.RLI[i,2] <- 1-(sum(data.id.full26$id)/(5*length(spnames)))
  
  data.id.full60 <-data.id[match(RLI_index_1$full_rcp60, data.id$status),]
  data.RLI[i,3] <- 1-(sum(data.id.full60$id)/(5*length(spnames)))
  
  data.id.full85 <-data.id[match(RLI_index_1$full_rcp85, data.id$status),]
  data.RLI[i,4] <- 1-(sum(data.id.full85$id)/(5*length(spnames)))
  
  data.id.buffer26 <-data.id[match(RLI_index_1$buffer_rcp26, data.id$status),]
  data.RLI[i,5] <- 1-(sum(data.id.buffer26$id)/(5*length(spnames)))
  
  data.id.buffer60 <-data.id[match(RLI_index_1$buffer_rcp60, data.id$status),]
  data.RLI[i,6] <- 1-(sum(data.id.buffer60$id)/(5*length(spnames)))
  
  data.id.buffer85 <-data.id[match(RLI_index_1$buffer_rcp85, data.id$status),]
  data.RLI[i,7] <- 1-(sum(data.id.buffer85$id)/(5*length(spnames)))
  
  data.id.stable26 <-data.id[match(RLI_index_1$stable_rcp26, data.id$status),]
  data.RLI[i,8] <- 1-(sum(data.id.stable26$id)/(5*length(spnames)))
  
  data.id.stable60 <-data.id[match(RLI_index_1$stable_rcp60, data.id$status),]
  data.RLI[i,9] <- 1-(sum(data.id.stable60$id)/(5*length(spnames)))
  
  data.id.stable85 <-data.id[match(RLI_index_1$stable_rcp85, data.id$status),]
  data.RLI[i,10] <- 1-(sum(data.id.stable85$id)/(5*length(spnames)))
}


data.gap <- read.csv("./Analyses/New_combined/Gap/data.gap.csv")
RLI_index <- read.csv("./Analyses/New_combined/RLI_index.csv")

#### no dispersal ####

data <- RLI_index[which(RLI_index$Qin != "LC" & RLI_index$Qin != "DD" & RLI_index$Qin != "NT"),]
spnames.Qin <- data$code
data_future<- RLI_index[which(RLI_index$stable_rcp85 != "LC" & RLI_index$stable_rcp85 != "DD" & RLI_index$stable_rcp85 != "NT"),]
spnames.future <- data_future$code
res <- which(spnames.future %in% spnames.Qin == F)
spcode <- spnames.future[res]
status <- data_future$stable_rcp60[res]
data.gap <- data.gap[which(data.gap$sp %in% spcode),]
data.final <- tibble(sp = spcode, status = status, gap = data.gap$category)

write.csv(data.final,file="./Analyses/New_combined/Gap/RCP85/stable85.csv")





