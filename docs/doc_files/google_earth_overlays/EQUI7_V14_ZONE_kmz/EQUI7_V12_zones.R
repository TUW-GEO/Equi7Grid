## Fig. 6 in http://www.sciencedirect.com/science/article/pii/S0098300414001629

library(rgdal)
library(sp)
library(maptools)
setwd("D:\\SoilGrids1km\\Equi7_Grid_V12_Public_Package\\Grids")
lst <- list.files(pattern="*_PROJ_ZONE.shp$", full.names = TRUE, recursive = TRUE)
equi7 <- list(NULL)
for(i in 1:length(lst)){
  equi7[[i]] <- readOGR(lst[i], layer=strsplit(basename(lst[i]),"\\.")[[1]][1])
}

for(i in 1:length(lst)){
  names(equi7)[i] <- strsplit(strsplit(basename(lst[i]), "\\.")[[1]][1], "_")[[1]][3]
}
save(equi7, file="equi7.rda", compress="xz") 

lst.ll <- list.files(pattern="*_GEOG_ZONE.shp$", full.names = TRUE, recursive = TRUE)
equi7.ll <- list(NULL)
for(i in 1:length(lst.ll)){
  equi7.ll[[i]] <- readOGR(lst.ll[i], layer=strsplit(basename(lst.ll[i]),"\\.")[[1]][1])
}

library(plotKML)
kml_open("equi7.kml")
for(i in 1:length(lst)){
  kml_layer(equi7.ll[[i]], subfolder.name=strsplit(basename(lst.ll[i]),"\\.")[[1]][1], colour=ZONE, colour_scale=rep("#FFFF00", 2), alpha=.4)
}
kml_close("equi7.kml")
