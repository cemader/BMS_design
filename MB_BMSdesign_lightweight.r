
library(MBHdesign)
library(spsurvey)
library(maptools)
library(rgeos)
library(raster)
library(car)
library(sf)
library(raster)
library(tidyverse)
library(sp)
library(rgdal)
library(lakemorpho)
library(dplyr)
library(car)
library(doBy)
library(gdata)
library(ggplot2)
library(nngeo)
library(reshape2)
library(viridis)
library(scales)


setwd("C:/Users/maderc/Documents/BMS/Manitoba_selection")
#load("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
#setwd("C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Saskatchewan Breeding Bird Atlas/Boreal Site Selection/")



###################################################################################################################################
##########################################          IMPORT THE SAMPLING FRAME            ##########################################
###################################################################################################################################


# Albers equal-area projection parameters based on above
#aea.proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
laea.proj <- "+proj=laea +lat_0=0 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" # use this because it is lccc10 projection
NAD.proj <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs " #Manitoba hexagons are in this projection

hexagons<-readOGR("Mb_BOREAL_Hex_CM_UTM.shp", p4s= NAD.proj)
LCCC10 <- raster("C:/Users/maderc/Documents/BMS/Manitoba_selection/NA_LandCover_2010_25haMMU/NA_LandCover_2010_25haMMU.tif")
hexagons <- spTransform(hexagons, crs(LCCC10)) #transform to be in same projection as raster


#plot(hexagons, col=hexagons$Min_Cost) takes needless time
hexagons$cost_prob<-1/sqrt(hexagons$Min_Cost)


writeOGR(hexagons, "C:/Users/maderc/Documents/BMS/Manitoba_selection", "hexagons", driver = "ESRI Shapefile", overwrite = TRUE)

## now use the final_hexagons to clip lccc10 raster
# read in lccc10 file

###################------------    the following section is to create habitat selected MB.rast dataset.
##################------------     Can be removed after first successful creation of MB.rast after masking and removal of unwanted habitat classes
#################-------------     Temporarily replaced with the following load of dataset, because mask() takes around 45 minutes

MB.rast <- raster("C:/Users/maderc/Documents/BMS/Manitoba_selection/MB.rast.tif")
#plot(MB.rast)

#MB.rast <- crop(LCCC10,as(hexagons,"Spatial")) # LCCC raster for MB only
#MB.rast <- mask(MB.rast,hexagons,"Spatial", overwrite = TRUE)

# remove habitat classes not wanted in the model
#MB.rast <- calc(MB.rast, fun=function(x){x[x==19] <- NA; return(x)}) # remove snow and ice
#MB.rast <- calc(MB.rast, fun=function(x){x[x==18] <- NA; return(x)}) # water
#MB.rast <- calc(MB.rast, fun=function(x){x[x==17] <- NA; return(x)}) # urban
#MB.rast <- calc(MB.rast, fun=function(x){x[x==15] <- NA; return(x)}) # cropland
#MB.rast <- calc(MB.rast, fun=function(x){x[x==0] <- NA; return(x)}) # pixels not included in the cropping
#writeRaster(MB.rast,filename="MB.rast.tif",format="GTiff",overwrite=TRUE)

rm(LCCC10) # clear some memory

##############################################################################################################
############################ Now calculating Habitat Probabilities for ecoregion #############################


eco.list <- c("Aspen Parkland","Boreal Transition","Coastal Hudson Bay Lowland","Churchill River Upland",
              "Hudson Bay Lowland", "Hayes River Upland","Interlake Plain","Kazan River Upland", "Lake Manitoba Plain", 
              "Lac Seul Upland", "Lake of the Woods", "Mid-Boreal Lowland","Mid-Boreal Uplands","Maguse River Upland",
              "Selwyn Lake Upland" )




########################  ---- up to this point updated for MB Dec 4  --CM

##############################################################################################################
############################# Add hab prob and cost prob to NL hexagon layer ##################################

#eco.hab.hex <- raster::extract(HabProbraster,hexagons,fun=sum,na.rm=TRUE,df=TRUE,sp=TRUE) # Takes a long time. new hexagon file with hab probs
#writeOGR(eco.hab.hex, "C:/Users/maderc/Documents/BMS/Manitoba_selection", "ecohabhex", driver = "ESRI Shapefile") #Wrote in new shapefile, temporarily replaced this step with load

eco.hab.hex<-readOGR("ecohabhex.shp", p4s= laea.proj)

hexagons$hab_prob<-eco.hab.hex@data$layer



###### CALCULATE THE INCLUSION PROBABILITY
###### NOTE 1: The variable Balanced is the habitat based inclusion probabilities for Balanced draw (weighted in Mahon parlance) 
###### NOTE 2: The variable invsqrCost is the inverse square root of the minimum cost (across differeint methods) to access the hexagon
###### NOTE 3: HOLDING HABITAT CONSTANT SO IT IS EASY TO TEST WHETHER COST IS HAVING ANY EFFECT

hexagons$p<-(hexagons$hab_prob*hexagons$cost_prob)/sum(hexagons$hab_prob*hexagons$cost_prob)

### dividing by sum puts on probability scale (i.e. all hexagons will sum to 1)

hexagons<-subset(hexagons, p!=0) #remove zeroes
plot(hexagons, col= hexagons$p)

#row.names(hexagons) <- seq(1,nrow(hexagons))  ### this will help for binding data later

rm(HabProbraster)
rm(eco.hab.hex)
rm(MB.rast)

##### Target sample sizes
#samplesize<-c(57,20,34,67,12,11,1,11)  # target sample sizes for priority squares
samplesize<-c(3,1,70,24,95,44,36,16,50,1,14,56,15,1,77)  # CODE TEST PLACEHOLDER #### target sample sizes for priority squares BASED on N=400  
Ecoregion<-c("Aspen Parkland","Boreal Transition","Coastal Hudson Bay Lowland","Churchill River Upland",
             "Hudson Bay Lowland", "Hayes River Upland","Interlake Plain","Kazan River Upland", "Lake Manitoba Plain", 
             "Lac Seul Upland", "Lake of the Woods", "Mid-Boreal Lowland","Mid-Boreal Uplands","Maguse River Upland",
             "Selwyn Lake Upland")
sample.size<-data.frame(samplesize, Ecoregion)
sample.size

### May need to add more oversamples for Boreal and Taiga Shield later



#sub<-subset(hexagons, LegacySite==1)
#table(sub$Ecoregion)
#sample.size$Legacy<- c(50,11,26,42,13,7,1,15)
#colnames(sample.size)<-c("Target","Ecoregion","Legacy")
#sample.size$samplesize<-sample.size$Target - sample.size$Legacy
#sample.size$samplesize<-ifelse(sample.size$samplesize< 1, 1,sample.size$samplesize)



#############################################################################################################################
###                                      
###                                      
###                                      		Select the ecoregion to draw from
###                                      
###                                      
#############################################################################################################################

### Total sample size target across ecoregions
N<- sum(sample.size$samplesize)


# inclusion probabilties 
p <- as.matrix(hexagons$p)
p <- N * p / sum( p) #standardise to get n samples  




#############################################################################################################################
###                                      
###                                      
###                                      TEST THE COMPETING DESIGNS by drawing repeated random samples
###                                      
###                                      
#############################################################################################################################



writeSpatialShape(hexagons, "Sample_Frame")



rm(hexagons)


set.seed(3911548)  # Set seed 
# create mdcaty.  keeps variable it is based on in output in addition to mdcaty.  useful sometimes
# also want sum of inclusion probabilities to sum to expected sample size.
#attframe <- read.dbf("Sample_Frame")
attframe <- read_sf("Sample_Frame.shp") 

summary(attframe$p)


attframe$mdcaty <- N * attframe$p/sum(attframe$p)  

sum(attframe$mdcaty)



###  SPECIFY THE STRATIFIED SAMPLING DESIGN
oversample.size = 1 ### 2x oversample


################### RESTART HERE
# DESIGN 3: FULL BMS DESIGN -----------------------------------------------------

attframe$mdcaty <- N * attframe$p/sum(attframe$p)


#### STRATIFIED DESIGN WITH THE USE OF INCLUSION PROBABILITIES
# NO OVERSAMPLE because just testing the design with Panel one
Stratdsgn <- list("Aspen Parkland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Aspen Parkland")$samplesize), over=0, seltype="Continuous"),
                  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0, seltype="Continuous"),
                  "Coastal Hudson Bay Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Coastal Hudson Bay Lowland")$samplesize), over=0,seltype="Continuous"),
                  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0,seltype="Continuous"),
                  "Hudson Bay Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Hudson Bay Lowland")$samplesize), over=0,seltype="Continuous"),
                  "Hayes River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Hayes River Upland")$samplesize), over=0,seltype="Continuous"),
                  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=0,seltype="Continuous"),
                  "Kazan River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Kazan River Upland")$samplesize), over=0,seltype="Continuous"),
                  "Lake Manitoba Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Lake Manitoba Plain")$samplesize), over=0,seltype="Continuous"),
                  "Lac Seul Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Lac Seul Upland")$samplesize), over=0,seltype="Continuous"),
                  "Lake of the Woods"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Lake of the Woods")$samplesize), over=0,seltype="Continuous"),
                  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Continuous"),
                  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Continuous"),
                  "Maguse River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Maguse River Upland")$samplesize), over=0,seltype="Continuous"),
                  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=0,seltype="Continuous")
)
Stratdsgn

#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=0,seltype="Continuous"),

#Stratdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize*oversample.size,0), seltype="Continuous"),
#  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize*oversample.size,0), seltype="Continuous"),
#  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Interlake Plain")$samplesize*oversample.size,0),seltype="Continuous"),
#  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=round(subset(sample.size, Ecoregion=="Boreal Transition")$samplesize*oversample.size,0),seltype="Continuous")
#)
#Stratdsgn


n.simulations=3 #### run 100 random draws
full.design.samp<-list() ###


#### this differs because it saves every iteration, other designs ONLY saving the last iteration
for (i in 1:n.simulations){
  full.design.samp[i] <- grts(design=Stratdsgn,
                              DesignID="ET_ID",
                              type.frame="area",
                              src.frame="shapefile",
                              in.shape="Sample_Frame",
                              att.frame="attframe",
                              stratum="REGION_NAM",
                              mdcaty="mdcaty",
                              shapefile=FALSE,
                              out.shape=paste0("Full_BMS_design", i)) 
}

#for (i in 1:n.simulations){
#full.design.samp[i] <- grts(design=Stratdsgn,
#DesignID="ET_ID",
#type.frame="area",
#src.frame="shapefile",
#in.shape="Sample_Frame",
#att.frame=attframe,
#stratum="Ecoregion",
#mdcaty="mdcaty",
#shapefile=TRUE,
#out.shape="Full_BMS_design") 
#}


#### get cost estimates from the draws (above)
cost.est<-rep (NA, length(full.design.samp))
for (i in 1:length(full.design.samp)){
  cost.est[[i]] <- sum((full.design.samp[[i]])$Cost)
}

sp.balance<-rep (NA, length(full.design.samp))
for (i in 1:length(full.design.samp)){
  sp.balance[[i]] <- as.matrix(unlist(spbalance(full.design.samp[[i]], frame, tess_ind = F, sbc_ind = T)$sbc$J_subp))[1,1] 
}



hab = list()
for (i in 1:length(full.design.samp)) {
  # ... make some data
  dat <- data.frame(VALUE_1= sum(full.design.samp[[i]]@data$VALUE_1),VALUE_2= sum(full.design.samp[[i]]@data$VALUE_2),VALUE_5= sum(full.design.samp[[i]]@data$VALUE_5),VALUE_6= sum(full.design.samp[[i]]@data$VALUE_6),VALUE_8= sum(full.design.samp[[i]]@data$VALUE_8),VALUE_10= sum(full.design.samp[[i]]@data$VALUE_10),VALUE_14= sum(full.design.samp[[i]]@data$VALUE_14),VALUE_15= sum(full.design.samp[[i]]@data$VALUE_15),VALUE_16= sum(full.design.samp[[i]]@data$VALUE_16),VALUE_17= sum(full.design.samp[[i]]@data$VALUE_17),VALUE_18= sum(full.design.samp[[i]]@data$VALUE_18))
  dat$i <- i  # maybe you want to keep track of which iteration produced it?
  hab[[i]] <- dat # add it to your list
}

habitat = do.call(rbind, hab)


Design.data.full.design.draw<-data.frame(Design=rep("Full BMS design", length(cost.est)), Cost=cost.est, Spatial.Balance= sp.balance)
Design.data.full.design.draw<-data.frame(Design.data.full.design.draw, habitat)
Design.data.full.design.draw
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")


# DESIGN 4: SPATIAL BALANCE ONLY ------------------------------------------
attframe <- read.dbf("Sample_Frame")
attframe <-read_sf("Sample_Frame.shp")

summary(attframe$p)


attframe$mdcaty <- N * attframe$p/sum(attframe$p)  
sum(attframe$mdcaty)

# NO OVERSAMPLE because just testing the design with Panel one
Spatdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=0, seltype="Equal"),
                 "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=0, seltype="Equal"),
                 "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=0,seltype="Equal"),
                 "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=0,seltype="Equal"),
                 "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=0,seltype="Equal"),
                 "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=0,seltype="Equal"),
                 "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=0,seltype="Equal")
)
Spatdsgn

#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=0,seltype="Equal"),

#Spatdsgn <- list("Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Athabasca Plain")$samplesize*oversample.size,0), seltype="Equal"),
#  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Churchill River Upland")$samplesize*oversample.size,0), seltype="Equal"),
#  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Lowland")$samplesize*oversample.size,0),seltype="Equal"),
#  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize), over=round(subset(sample.size, Ecoregion=="Mid-Boreal Uplands")$samplesize*oversample.size,0),seltype="Equal"),
#  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Selwyn Lake Upland")$samplesize*oversample.size,0),seltype="Equal"),
#  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize), over=round(subset(sample.size, Ecoregion=="Tazin Lake Upland")$samplesize*oversample.size,0),seltype="Equal"),
#  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Interlake Plain")$samplesize), over=round(subset(sample.size, Ecoregion=="Interlake Plain")$samplesize*oversample.size,0),seltype="Equal"),
#  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, Ecoregion=="Boreal Transition")$samplesize), over=round(subset(sample.size, Ecoregion=="Boreal Transition")$samplesize*oversample.size,0),seltype="Equal")
#)
#Spatdsgn


n.simulations=100 #### run 100 random draws
spatial.design.samp<-list() ###

for (i in 1:n.simulations){
  spatial.design.samp[i] <- grts(design=Spatdsgn,
                                 DesignID="ET_ID",
                                 type.frame="area",
                                 src.frame="shapefile",
                                 in.shape="Sample_Frame",
                                 att.frame=attframe,
                                 stratum="Ecoregion",
                                 mdcaty="mdcaty",
                                 shapefile=TRUE,
                                 out.shape="Spatial_design") 
}



#### get cost estimates from the draws (above)
cost.est<-rep (NA, length(spatial.design.samp))
for (i in 1:length(spatial.design.samp)){
  cost.est[[i]] <- sum((spatial.design.samp[[i]])$Cost)
}

sp.balance<-rep (NA, length(spatial.design.samp))
for (i in 1:length(spatial.design.samp)){
  sp.balance[[i]] <- as.matrix(unlist(spbalance(spatial.design.samp[[i]], frame, tess_ind = F, sbc_ind = T)$sbc$J_subp))[1,1] 
}



hab = list()
for (i in 1:length(spatial.design.samp)) {
  # ... make some data
  dat <- data.frame(VALUE_1= sum(spatial.design.samp[[i]]@data$VALUE_1),VALUE_2= sum(spatial.design.samp[[i]]@data$VALUE_2),VALUE_5= sum(spatial.design.samp[[i]]@data$VALUE_5),VALUE_6= sum(spatial.design.samp[[i]]@data$VALUE_6),VALUE_8= sum(spatial.design.samp[[i]]@data$VALUE_8),VALUE_10= sum(spatial.design.samp[[i]]@data$VALUE_10),VALUE_14= sum(spatial.design.samp[[i]]@data$VALUE_14),VALUE_15= sum(spatial.design.samp[[i]]@data$VALUE_15),VALUE_16= sum(spatial.design.samp[[i]]@data$VALUE_16),VALUE_17= sum(spatial.design.samp[[i]]@data$VALUE_17),VALUE_18= sum(spatial.design.samp[[i]]@data$VALUE_18))
  dat$i <- i  # maybe you want to keep track of which iteration produced it?
  hab[[i]] <- dat # add it to your list
}

habitat = do.call(rbind, hab)



Design.data.spatial.draw<-data.frame(Design=rep("Spatial only", length(cost.est)), Cost=cost.est, Spatial.Balance= sp.balance)
Design.data.spatial.draw<-data.frame(Design.data.spatial.draw, habitat)
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
Design.data.spatial.draw


Design.data<-rbind(Design.data.cost.draw,Design.data.habitat.draw,Design.data.full.design.draw,Design.data.spatial.draw)
Design.data
save.image("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
write.csv(Design.data, file = "Design_comparison_largeN.csv") 

#library(car)
boxplot(Cost~Design, data=Design.data, ylab="Program Cost")
Design.data$factor<-factor(Design.data$Design, levels=c("Cost only", "Full BMS design", "Spatial only", "Habitat only"))

tiff(filename = "C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS Manuscript/SK_Cost_figure.tif",width = 174, height = 174, units = "mm",  res=600, compression = c("lzw"),restoreConsole = TRUE)
boxplot(Cost~factor, data=Design.data, ylab="Program Cost over Ten Years ($)", yaxt="n")
axis(2, at=axTicks(2), labels=sprintf("$%s", axTicks(2)))
dev.off()

tiff(filename = "C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS Manuscript/SK_spatial_balance_figure.tif",width = 174, height = 174, units = "mm",  res=600, compression = c("lzw"),restoreConsole = TRUE)
boxplot(Spatial.Balance~factor, data=Design.data, ylab="Spatial Balance")
dev.off()

tiff(filename = "C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Boreal Monitoring Strategy/BMS Manuscript/SK_Habitat_figure.tif",width = 174, height = 174, units = "mm",  res=600, compression = c("lzw"),restoreConsole = TRUE)
Boxplot(SSD~Design, data=Designs, ylab="Sum of Squared Differences", id.method="none")
dev.off()
