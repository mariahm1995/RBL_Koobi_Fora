####Carn Ecometric Koobi Fora projections
#Updated 11/7/2024

#############RBL Analyses
#Worldwide & N. America
#Precip (log) & temp
#Sensitivity analysis
#Seven paleo sites
#Updated 12/8/2023


#install.packages(c("sp", "raster", "sf", "ggplot2"))
library(sp)
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
# library(rgdal)

# Input data ----

#Polygon shape files for continental shape
continents <- shapefile("inputs/continents.shp")
continents.sf <- st_as_sf(continents)

#All global environmental data at 50km equidistant points
points <- read.csv("inputs/S_L_allpointsdata2.csv")
head(points)
dim(points)

##convert sampling points to spatial points and sample species at each sampling locality
sp <- SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(continents))
sp <- raster::intersect(sp, continents)
head(sp)
length(sp)

#Save continent to points
points$cont<- sp@data$CONTINENT

##Load RBL trait data & format file
traits <- read.csv("inputs/RBL_global_Dec2023.csv", header = T, na= ".", stringsAsFactors = FALSE)
traits$TaxonName <- gsub("_", " ", traits$TaxonName)
rownames(traits) <- as.character(traits$TaxonName)
# colnames(traits) <- c("binomial", "RBL", "IUCN", "Vert", "Non-vert", "Diet", 
#                       "Diet_vert", "Body_size", "Diet breadth", "Home Range", 
#                       "Geographic Range", "Publication", "Sample size", "Notes")
colnames(traits) <- c( "binomial", "RBL", "IUCN", "Vert", "Non-vert",
                       "Diet", "Diet_vert", "Body_size", "Diet_breadth", "Home_range",     
                       "Geog_range", "Museum_Pub", "N", "Notes", "X",
                       "Invert", "Scav", "Fruit", "Nect", "Seed",
                       "Plant", "Vert2", "Diet_old", "Diet_vert2")
traits <- traits[, !names(traits) %in% "X"]

# Trait value graphs ----

#Diet_diversity <- diet vert
names(traits)
table(traits$Diet)

#Rank order plot
traits$binomial <- factor(traits$binomial, levels = traits$binomial[order(traits$RBL)])
rankorderspecies <- ggplot(traits, aes(x = binomial, y = RBL)) +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T), linetype = "dashed", color = "black") +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T) - sd(traits$RBL), linetype = "solid", color = "black") +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T) + sd(traits$RBL), linetype = "solid", color = "black") +
  geom_segment(aes(xend = binomial, yend = 1.3)) +
  geom_point(aes(x = binomial, y = RBL), color = "black", pch = 16, cex = 2) +
  labs(x = "Species", y = "RBL") +
  scale_y_continuous(expand = c(0, 0), limits = c(0.3, 1.1), breaks = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)) +
  theme(axis.line = element_line(color = "black"), 
        panel.background = element_rect(color = "black", fill = NA), 
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", angle = 90, size = 8, hjust = 1, vjust = 0.25),
        axis.title = element_text(color = "black", size = 14))
rankorderspecies

# Graph of 
mid <- mean(traits$Vert)

traits$binomial <- factor(traits$binomial, levels = traits$binomial[order(traits$RBL)])
ggplot(traits, aes(x = binomial, y = RBL, col = Vert)) +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T), linetype = "dashed", color = "black") +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T) - sd(traits$RBL), linetype = "solid", color = "black") +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T) + sd(traits$RBL), linetype = "solid", color = "black") +
  geom_segment(aes(x = binomial, xend = binomial, y = 0, yend = RBL), lwd=0.3) +
  geom_point(aes(x = binomial, y = RBL), pch = 16, cex = 2.5) +
  scale_color_gradient2(midpoint = mid, low = "turquoise3", mid = "bisque2", high = "darkorchid4", space = "Lab", name = "Percent \nvertebrate \nprey") +
  labs(x = "Species", y = "RBL") +
  scale_x_discrete(breaks = levels(traits$binomial)[c(T, rep(F, 3))]) +
  theme(axis.line = element_line(color = "black"), 
        panel.background = element_rect(color = "black", fill = NA), 
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", angle = 90, size = 15, hjust = 1, vjust = 0.25, face="italic"),
        axis.title = element_text(color = "black", size = 24)) +
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12)) +
  theme(legend.key.size = unit(2, 'lines'))

#save as 1500 x 700
scale_color_gradient(low="turquoise", high="darkmagenta", name = "Vertebrate Prey") 

# Geography ----
#geography <- shapefile("TERRESTRIAL_MAMMALS.shp")
#save(geography, file="geography.Rdata")
load("inputs/geography.Rdata")
head(geography)
geography_carn <- subset(geography, order_ == "CARNIVORA")
head(geography)

#Sample species at each sampling location.
# o <- over(sp, geography_carn, returnList = T)
# save(o, file="o.Rdata")
load("inputs/o.Rdata")
head(o)
length(o)

points.iucn <- points

#Log transformation of precipitation to represent habitat patterns
points.iucn$LogPrecip <- log(points$BIO12+1)

#Summarize traits for community-level distribution: calculating mean & standard deviation
points.iucn$mean_RBL <- unlist(lapply(o, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
hist(points.iucn$mean_RBL)
points.iucn$sd_RBL <- unlist(lapply(o, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))
hist(points.iucn$sd_RBL)

#Overall mammal species richness at each point
points.iucn$RBL_richness <- unlist(lapply(o, function(x) length(traits[x$binomial, "RBL"])))
hist(points.iucn$RBL_richness)

#Calculate richness at each point with traits data
points.iucn$count <- unlist(lapply(o, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))
hist(points.iucn$count)

#Min 3 species per community
points.rest3 <- points.iucn[points.iucn$count > 2, ]
head(points.rest3)
dim(points.rest3)

# Global analysis ----

#Bin distributions 25x25
mRBL3 <- range(points.rest3$mean_RBL, na.rm=T)
sdRBL3 <- range(points.rest3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrks3 <- seq(mRBL3[1]-0.001, mRBL3[2]+0.001, diff(mRBL3)/25)
sdbrks3 <- seq(sdRBL3[1]-0.001, sdRBL3[2]+0.001, diff(sdRBL3)/25)

#Assign bin codes
mbc3<-.bincode(points.rest3$mean_RBL, breaks = mbrks3)
sdbc3 <-.bincode(points.rest3$sd_RBL, breaks = sdbrks3)

#Assign mean & sd bin code to each community
points.rest3$mbc3<-mbc3
points.rest3$sdbc3<-sdbc3

###################Precipitation
mod <- list()
obj <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbc3 == i & sdbc3 == j, na.rm = T) != 0) {
      dat <- points.rest3$LogPrecip[which(mbc3 == i & sdbc3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}


##make raster for body mass and temperature
r <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r <- setValues(r, obj)
rdf <- as.data.frame(r, xy=TRUE)
names(rdf) <- c("Mean", "SD", "Precip")
head(rdf)

breaksm <- c((mbrks3[1]+(mRBL3[2]-mRBL3[1])/25),
             (mbrks3[13]+(mRBL3[2]-mRBL3[1])/25),
             (mbrks3[24]+(mRBL3[2]-mRBL3[1])/25))
breaksm <- trunc(breaksm*10^2)/10^2

breaksd <- c((sdbrks3[1]+(sdbrks3[2]-sdbrks3[1])/25), 
             (sdbrks3[13]+(sdbrks3[2]-sdbrks3[1])/25),
             (sdbrks3[24]+(sdbrks3[2]-sdbrks3[1])/25))
breaksd <- trunc(breaksd*10^2)/10^2

#Ecometric space for precipitation
color <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r) - minValue(r)))
precip_ecospace <-ggplot(data = rdf) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color), name = "Precipitation (mm)", 
                       na.value = "transparent", 
                       breaks = waiver(), labels = waiver()) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.44), labels = breaksd) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.5), labels = breaksm) +
  coord_fixed() +
  theme_bw()
precip_ecospace

######Calculate maximum likelihood for all bins
modmax3 <- array(NA, dim = length(points.rest3[ , 1]))
mod <- list()
for (i in 1:length(points.rest3[ , 1])) {
  if(!(is.na(points.rest3$mbc3[i]) | is.na(points.rest3$sdbc3[i]))) {
    dat <- points.rest3$LogPrecip[which(points.rest3$mbc3 == points.rest3$mbc3[i] & points.rest3$sdbc3 == points.rest3$sdbc3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax3

##histogram of estimated precipitation based on community trait values
hist(modmax3)

#calculate precipitation anomaly
anom_precip <- points.rest3$LogPrecip - modmax3
#Take away the log (nL) so that the value is more easily interpretable
anom_precip_NL <- exp(anom_precip) - 1
plot(anom_precip, anom_precip_NL)

points.rest3$precipanom <- as.numeric(anom_precip)
summary(points.rest3$precipanom, na.rm = T)
points.rest3$precipest <- as.numeric(modmax3)
summary(anom_precip_NL, na.rm = T)

#Correlation of observed ~ expected precipitation
precipcor <- lm(points.rest3$precipest ~ points.rest3$LogPrecip)
precipcor
summary(precipcor)
cor_precipest <- cor.test(points.rest3$precipest, points.rest3$LogPrecip, method = "pearson")
#cor=0.7007268 
cor_precipest

##Actual precip map
sp_points <- SpatialPointsDataFrame(points.rest3[,2:3],points.rest3)
sf_points <- st_as_sf(sp_points)

act_precip<-ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void()  +
  geom_sf(data = sf_points, aes(color = LogPrecip), size = 2, pch = 16) +
  scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Observed \nPrecipitation \n Log (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
  labs(color = "LogPrecip") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
act_precip

#RBL variation map
color_RBL <- c()
rbl_s_map<-ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = sd_RBL), size = 2, pch = 16)  +
  scale_color_gradientn(limits=as.numeric(c(0, 0.30)), colors = rev(color_RBL), name = "RBL \nStandard deviation", na.value = "transparent", breaks = c(0, 0.1, 0.2, 0.3), labels = c(0, 0.1, 0.2, 0.3)) +
  labs(color = "sd_RBL") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") 
rbl_s_map

##Map of predicted precip
est_precip<-ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() + 
  geom_sf(data = sf_points, aes(color = modmax3), size = 2, pch = 16) +
  scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Estimated \nPrecipitation \n Log(mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
  labs(color = "LogPrecip") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
est_precip

##Anomaly map
anomlimits <- range(points.rest3$precipanom, na.rm = T)
anombreaks <- seq(anomlimits[1], anomlimits[2], 3)
anomlab <- round(anombreaks, 2)
#anombreaks <- c(-6.6, -4.4, -2.28, -0.12, 2.04, 4.2) #breaks from range(points.rest3$precipanom, na.rm=T)
#anomlab <- c(-6.6, -4.4, -2.28, -0.12, 2.04, 4.2)
coloranom <- colorRampPalette(c("aquamarine","lightyellow","steelblue"))(round(maxValue(r) - minValue(r)))
anommap <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = points.rest3$precipanom), size = 1.6, pch = 16) +
  scale_color_gradientn(colors = coloranom, name = "Precipitation \nAnomaly Log (mm)", na.value = "transparent", breaks = anombreaks, labels = anomlab) +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
anommap



