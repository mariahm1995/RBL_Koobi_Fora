####Carn Ecometric Koobi Fora projections
# Title:        [Descriptive Title of the Script]
# Description:  [Brief description of what the script does, e.g., analysis of X, 
#                generation of figures, data cleaning, etc.]
# Author:       Maria A. Hurtado-Materon and Leila Siciliano-Martina
# Last Updated: 2025-04-10
# R Version:    4.4.3 (2025-02-28) -- "Trophy Case"
# Script Type:  Data cleaning, analysis, visualization and model fitting
#
# Associated Paper: [Title of Paper]
#                   [Journal/Repository] (optional)
#                   [DOI or link if available]
#
# Data input:   [Brief note or path to input data]
# Output:       [Brief note or path to output files (plots, tables, models, etc.)]
#

# Dependence ----
#install.packages(c("sp", "raster", "sf", "ggplot2"))
library(sp)
library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

# Data input ----
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

#Diet_diversity <- diet vert
names(traits)
table(traits$Diet)

# Trait distribution graphs ----
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

# Ecometrics model precipitation ----

#geography <- shapefile("TERRESTRIAL_MAMMALS.shp")
#save(geography, file="geography.Rdata")
load("inputs/geography.Rdata")
geography_carnivora <- geography[geography$order_ == "CARNIVORA", ]

#Sample species at each sampling location.
# o <- over(sp, geography_carnivora, returnList = T)
# save(o, file="oCarnivora.Rdata")
load("outputs/oCarnivora.Rdata")
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

#Richness of just carnivoran species
# carn_rich <- array(NA, dim = length(o))
# 
# for(i in 1:length(o)){
#   carn_rich[i] <- sum(o[[i]]$order_=="CARNIVORA")
# }
# hist(carn_rich)
# summary(carn_rich)
# points.iucn$carn_rich <- carn_rich
# summary(points.iucn$carn_rich)

#overall carnivoran richness (excluding zero-species communities)
# carn_no_zero <- points.iucn$carn_rich[points.iucn$carn_rich > 0]
# summary(carn_no_zero)
# carn_three_plus <- points.iucn$carn_rich[points.iucn$count > 2]
#Global summary of carnivoran species richness in communities w/ at least 3 carnivoran species
# summary(carn_three_plus)

####Global 
#Min 3 species per community
points.rest3 <- points.iucn[points.iucn$count > 2, ]
head(points.rest3)
dim(points.rest3)

summary(points.rest3$mean_RBL)
sd(points.rest3$mean_RBL)

#Bin distributions 25x25
mRBL3 <- range(points.rest3$mean_RBL, na.rm = T)
sdRBL3 <- range(points.rest3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrks3 <- seq(mRBL3[1] - 0.001, mRBL3[2] + 0.001, diff(mRBL3)/25)
sdbrks3 <- seq(sdRBL3[1] - 0.001, sdRBL3[2] + 0.001, diff(sdRBL3)/25)

#Assign bin codes
mbc3 <- .bincode(points.rest3$mean_RBL, breaks = mbrks3)
sdbc3 <- .bincode(points.rest3$sd_RBL, breaks = sdbrks3)

#Assign mean & sd bin code to each community
points.rest3$mbc3 <- mbc3
points.rest3$sdbc3 <- sdbc3

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
color <- c("#bc6c25", "#fefae0", "#606c38")
precip_ecospace <-ggplot(data = rdf) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color), name = "Precipitation (loge(mm))", 
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
# color_RBL <- c()
# rbl_s_map<-ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void() +
#   geom_sf(data = sf_points, aes(color = sd_RBL), size = 2, pch = 16)  +
#   scale_color_gradientn(limits=as.numeric(c(0, 0.30)), colors = rev(color_RBL), name = "RBL \nStandard deviation", na.value = "transparent", breaks = c(0, 0.1, 0.2, 0.3), labels = c(0, 0.1, 0.2, 0.3)) +
#   labs(color = "sd_RBL") +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") 
# rbl_s_map

##Map of predicted precip
# est_precip <- ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void() + 
#   geom_sf(data = sf_points, aes(color = modmax3), size = 2, pch = 16) +
#   scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Estimated \nPrecipitation \n Log(mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
#   labs(color = "LogPrecip") +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
#   theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
#         legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
# est_precip

##Anomaly map
# anomlimits <- range(points.rest3$precipanom, na.rm = T)
# anombreaks <- seq(anomlimits[1], anomlimits[2], 3)
# anomlab <- round(anombreaks, 2)
# #anombreaks <- c(-6.6, -4.4, -2.28, -0.12, 2.04, 4.2) #breaks from range(points.rest3$precipanom, na.rm=T)
# #anomlab <- c(-6.6, -4.4, -2.28, -0.12, 2.04, 4.2)
# coloranom <- colorRampPalette(c("aquamarine","lightyellow","steelblue"))(round(maxValue(r) - minValue(r)))
# anommap <- ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void() +
#   geom_sf(data = sf_points, aes(color = points.rest3$precipanom), size = 1.6, pch = 16) +
#   scale_color_gradientn(colors = coloranom, name = "Precipitation \nAnomaly Log (mm)", na.value = "transparent", breaks = anombreaks, labels = anomlab) +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
# anommap

# Fossil analysis precipitation ----
fossil <- read.csv("inputs/Koobi_fora_full.csv", header = T)
dim(fossil)
fossil <- fossil[1:66, 1:9]

continents1 <- continents
continents1 <- st_crs(4326)

#Fossil map
fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))
fossmap <- act_precip + geom_sf(data = fossil.sf, col = "black", size = 2, pch = 17)
fossmap

#aggregating fossil data
fossilmean <- aggregate(fossil[ , 9], list(fossil$Time, fossil$Lat, fossil$Long), mean, na.rm = TRUE)
fossilsd <- aggregate(fossil[ , 9], list(fossil$Time), sd, na.rm = TRUE)
colnames(fossilmean) <- c("Time", "Lat", "Long", "Mean")
fossilmean <- merge(fossilmean, fossilsd, by.x = "Time", by.y = "Group.1")
colnames(fossilmean)[5] <- "SD"
fossildata <- fossilmean
head(fossildata)

fossildata$fossilmbc1 <- .bincode(fossildata$Mean, breaks = mbrks3)
fossildata$fossilsdbc1 <- .bincode(fossildata$SD, breaks = sdbrks3)

#calculate fossil precip estimates with 5% limits
fossilmodmax <- list() 
fossilmod <- list()
for (i in 1:length(fossildata[ , 1])) {
  if(!(is.na(fossildata$fossilmbc1[i]) | is.na(fossildata$fossilsdbc1[i]))) {
    if(length(which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])) < 1) {
      fossilmod[[i]] <- NA 
    } else {
      dat <- points.rest3$LogPrecip[which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])]
      fossilmod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      fossilmodmax$precipest[i] <- fossilmod[[i]]$x[which.max(fossilmod[[i]]$y)]
      maxlength <- length(which(fossilmod[[i]]$x < fossilmodmax$precipest[i]))
      bound.05 <- (length(fossilmod[[i]]$x)) * 0.05
      fossilmodmax$minlimit[i] <- fossilmod[[i]]$x[maxlength - bound.05]
      fossilmodmax$maxlimit[i] <- fossilmod[[i]]$x[maxlength + bound.05]
    }
  }
}

fossildata$precipest <- c(fossilmodmax$precipest, NA)
fossildata$minlimit <- c(fossilmodmax$minlimit, NA)
fossildata$maxlimit <- c(fossilmodmax$maxlimit, NA)

fossildata$precipestUN <- c(exp(fossilmodmax$precipest)-(1), NA)
fossildata$minlimitUN <- c(exp(fossilmodmax$minlimit)-(1), NA)
fossildata$maxlimitUN <- c(exp(fossilmodmax$maxlimit)-(1), NA)

#plot fossil sites in ecometric space
df1 <- data.frame(xmin = fossildata$fossilmbc1 - 1,
                  xmax = fossildata$fossilmbc1,
                  ymin = fossildata$fossilsdbc1 - 1,
                  ymax = fossildata$fossilsdbc1)

text <- data_frame(point = c("A","B","C","D","E","F"),
                   x = fossildata$fossilmbc1 -.5,
                   y = fossildata$fossilsdbc1 + .5)
text$x[c(2)] <- text$x[c(2)] - 1
text$x[c(5)] <- text$x[c(5)] +1

fossilEMS <- precip_ecospace + 
  geom_rect(data = df1, aes(xmin = as.numeric(xmin), 
                            xmax = as.numeric(xmax), ymin = as.numeric(ymin),
                            ymax = as.numeric(ymax)), colour = "black", 
            alpha = 0, linewidth = 0.75) + geom_text(text, mapping = aes(x = x,
                                                                         y = y,
                                                                         label = point))
fossilEMS

#select sampling point and observed precipitation nearest to fossil site
fossil.sf <- st_as_sf(fossildata, coords = c(3, 2), crs = crs(continents1))

nearest <- list()
near.point <- list()
for(i in 1:dim(fossildata)[1]) {
  nearest[[i]] <- sf_points[which.min(st_distance(sf_points, fossil.sf[i,])), ]
  near.point <- rbind(near.point, unlist(nearest[[i]], recursive = FALSE)[1])
}

fossildata$near.point <- as.numeric(near.point)
fossildata <- merge(fossildata, sf_points, by.x = "near.point", by.y = "GlobalID")
fossildata

write.csv(fossildata, "outputs/fossildataGlobal.csv", row.names=TRUE)

fossilEMSmod <- fossilEMS + geom_rect(data = fossildata, 
                                      aes(xmin = as.numeric(mbc3) - 1,
                                          xmax = as.numeric(mbc3),
                                          ymin = as.numeric(sdbc3) - 1,
                                          ymax = as.numeric(sdbc3)), 
                                      colour = "#c44536", alpha = 0, size = 1)
fossilEMSmod

# Geographic space of fossil bin map
act_precip_fossil <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = LogPrecip), size = 2, pch = 16) +
  scale_color_gradientn(
    limits = c(-2, 10), 
    colors = rev(color), 
    name = "Observed precipitation Log(mm)", 
    na.value = "transparent", 
    breaks = c(0, 2, 4, 6, 8, 10), 
    labels = c(0, 2, 4, 6, 8, 10)
  ) +
  labs(color = "LogPrecip") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  theme(
    legend.position = "top",  # ← move to top
    legend.direction = "horizontal",  # ← make horizontal
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", colour = "white")
  )

Time_A <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[1]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[1])
Time_BDE <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[2]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[2])
Time_C <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[3]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[3])

Time_A_points <- SpatialPointsDataFrame(Time_A[,2:3], Time_A)
sf_A_modern_points <- st_as_sf(Time_A_points)

Time_BDE_points <- SpatialPointsDataFrame(Time_BDE[,2:3], Time_BDE)
sf_BDE_modern_points <- st_as_sf(Time_BDE_points)

Time_C_points <- SpatialPointsDataFrame(Time_C[,2:3], Time_C)
sf_C_modern_points <- st_as_sf(Time_C_points)


Time_A_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_A_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  # scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Observed \nPrecipitation \n Log (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
  labs(title = "Time A") + theme(legend.position="none") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
  # theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
  #       legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))

Time_BDE_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_BDE_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  labs(title = "Times B D E") + theme(legend.position="none") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")

Time_C_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_C_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  labs(title = "Time C") + theme(legend.position="none")+
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")

combined_plot <- act_precip_fossil + Time_A_fig + Time_BDE_fig+ Time_C_fig

ggsave("figures/combined_plot_prec.png", plot = combined_plot, width = 12, height = 8, dpi = 300)

# Ecometrics model temperature ----

#Log transformation of precipitation to represent habitat patterns
points.rest3$temp <- points.rest3$BIO1/10
hist(points.rest3$temp)

mod <- list()
obj <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbc3 == i & sdbc3 == j, na.rm = T) != 0) {
      dat <- points.rest3$temp[which(mbc3 == i & sdbc3 == j)]
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
names(rdf) <- c("Mean", "SD", "Temp")
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
color <- c("#457b9d", "#fefae0", "#e63946")
temp_ecospace <-ggplot(data = rdf) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = color, name = "Temperature (C)", 
                       na.value = "transparent", 
                       breaks = waiver(), labels = waiver()) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.44), labels = breaksd) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.5), labels = breaksm) +
  coord_fixed() +
  theme_bw()
temp_ecospace

######Calculate maximum likelihood for all bins
modmax3 <- array(NA, dim = length(points.rest3[ , 1]))
mod <- list()
for (i in 1:length(points.rest3[ , 1])) {
  if(!(is.na(points.rest3$mbc3[i]) | is.na(points.rest3$sdbc3[i]))) {
    dat <- points.rest3$temp[which(points.rest3$mbc3 == points.rest3$mbc3[i] & points.rest3$sdbc3 == points.rest3$sdbc3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax3

##histogram of estimated precipitation based on community trait values
hist(modmax3)

#calculate precipitation anomaly
anom_temp <- points.rest3$temp - modmax3
#Take away the log (nL) so that the value is more easily interpretable
# anom_precip_NL <- exp(anom_precip) - 1
# plot(anom_precip, anom_precip_NL)

points.rest3$anom_temp <- as.numeric(anom_temp)
summary(points.rest3$anom_temp, na.rm = T)
points.rest3$tempest <- as.numeric(modmax3)
summary(anom_precip_NL, na.rm = T)

#Correlation of observed ~ expected precipitation
tempcor <- lm(points.rest3$tempest ~ points.rest3$temp)
tempcor
summary(tempcor)
cor_precipest <- cor.test(points.rest3$precipest, points.rest3$LogPrecip, method = "pearson")
#cor=0.7007268 
cor_precipest

##Actual precip map
sp_points <- SpatialPointsDataFrame(points.rest3[,2:3],points.rest3)
sf_points <- st_as_sf(sp_points)

# act_precip <- ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void()  +
#   geom_sf(data = sf_points, aes(color = LogPrecip), size = 2, pch = 16) +
#   scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Observed \nPrecipitation \n Log (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
#   labs(color = "LogPrecip") +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
#   theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
#         legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
# act_precip

#RBL variation map
# color_RBL <- c()
# rbl_s_map<-ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void() +
#   geom_sf(data = sf_points, aes(color = sd_RBL), size = 2, pch = 16)  +
#   scale_color_gradientn(limits=as.numeric(c(0, 0.30)), colors = rev(color_RBL), name = "RBL \nStandard deviation", na.value = "transparent", breaks = c(0, 0.1, 0.2, 0.3), labels = c(0, 0.1, 0.2, 0.3)) +
#   labs(color = "sd_RBL") +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") 
# rbl_s_map

##Map of predicted precip
# est_precip <- ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void() + 
#   geom_sf(data = sf_points, aes(color = modmax3), size = 2, pch = 16) +
#   scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Estimated \nPrecipitation \n Log(mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
#   labs(color = "LogPrecip") +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
#   theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
#         legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
# est_precip

##Anomaly map
# anomlimits <- range(points.rest3$precipanom, na.rm = T)
# anombreaks <- seq(anomlimits[1], anomlimits[2], 3)
# anomlab <- round(anombreaks, 2)
# #anombreaks <- c(-6.6, -4.4, -2.28, -0.12, 2.04, 4.2) #breaks from range(points.rest3$precipanom, na.rm=T)
# #anomlab <- c(-6.6, -4.4, -2.28, -0.12, 2.04, 4.2)
# coloranom <- colorRampPalette(c("aquamarine","lightyellow","steelblue"))(round(maxValue(r) - minValue(r)))
# anommap <- ggplot() +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
#   theme_void() +
#   geom_sf(data = sf_points, aes(color = points.rest3$precipanom), size = 1.6, pch = 16) +
#   scale_color_gradientn(colors = coloranom, name = "Precipitation \nAnomaly Log (mm)", na.value = "transparent", breaks = anombreaks, labels = anomlab) +
#   geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
# anommap

# Fossil analysis precipitation ----
fossil <- read.csv("inputs/Koobi_fora_full.csv", header = T)
dim(fossil)
fossil <- fossil[1:66, 1:9]

continents1 <- continents
continents1 <- st_crs(4326)

#Fossil map
fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))
fossmap <- act_precip + geom_sf(data = fossil.sf, col = "black", size = 2, pch = 17)
fossmap

#aggregating fossil data
fossilmean <- aggregate(fossil[ , 9], list(fossil$Time, fossil$Lat, fossil$Long), mean, na.rm = TRUE)
fossilsd <- aggregate(fossil[ , 9], list(fossil$Time), sd, na.rm = TRUE)
colnames(fossilmean) <- c("Time", "Lat", "Long", "Mean")
fossilmean <- merge(fossilmean, fossilsd, by.x = "Time", by.y = "Group.1")
colnames(fossilmean)[5] <- "SD"
fossildata <- fossilmean
head(fossildata)

fossildata$fossilmbc1 <- .bincode(fossildata$Mean, breaks = mbrks3)
fossildata$fossilsdbc1 <- .bincode(fossildata$SD, breaks = sdbrks3)

#calculate fossil precip estimates with 5% limits
fossilmodmax <- list() 
fossilmod <- list()
for (i in 1:length(fossildata[ , 1])) {
  if(!(is.na(fossildata$fossilmbc1[i]) | is.na(fossildata$fossilsdbc1[i]))) {
    if(length(which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])) < 1) {
      fossilmod[[i]] <- NA 
    } else {
      dat <- points.rest3$temp[which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])]
      fossilmod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      fossilmodmax$tempest[i] <- fossilmod[[i]]$x[which.max(fossilmod[[i]]$y)]
      maxlength <- length(which(fossilmod[[i]]$x < fossilmodmax$tempest[i]))
      bound.05 <- (length(fossilmod[[i]]$x)) * 0.05
      fossilmodmax$minlimit[i] <- fossilmod[[i]]$x[maxlength - bound.05]
      fossilmodmax$maxlimit[i] <- fossilmod[[i]]$x[maxlength + bound.05]
    }
  }
}

fossildata$tempest <- c(fossilmodmax$tempest, NA)
fossildata$minlimit <- c(fossilmodmax$minlimit, NA)
fossildata$maxlimit <- c(fossilmodmax$maxlimit, NA)

# fossildata$precipestUN <- c(exp(fossilmodmax$precipest)-(1), NA)
# fossildata$minlimitUN <- c(exp(fossilmodmax$minlimit)-(1), NA)
# fossildata$maxlimitUN <- c(exp(fossilmodmax$maxlimit)-(1), NA)

#plot fossil sites in ecometric space
df1 <- data.frame(xmin = fossildata$fossilmbc1 - 1,
                  xmax = fossildata$fossilmbc1,
                  ymin = fossildata$fossilsdbc1 - 1,
                  ymax = fossildata$fossilsdbc1)

text <- data_frame(point = c("A","B","C","D","E","F"),
                   x = fossildata$fossilmbc1 -.5,
                   y = fossildata$fossilsdbc1 + .5)
text$x[c(2)] <- text$x[c(2)] - 1
text$x[c(5)] <- text$x[c(5)] +1

fossilEMS <- temp_ecospace + 
  geom_rect(data = df1, aes(xmin = as.numeric(xmin), 
                            xmax = as.numeric(xmax), ymin = as.numeric(ymin),
                            ymax = as.numeric(ymax)), colour = "black", 
            alpha = 0, linewidth = 0.75) + geom_text(text, mapping = aes(x = x,
                                                                         y = y,
                                                                         label = point))
fossilEMS

#select sampling point and observed precipitation nearest to fossil site
fossil.sf <- st_as_sf(fossildata, coords = c(3, 2), crs = crs(continents1))

nearest <- list()
near.point <- list()
for(i in 1:dim(fossildata)[1]) {
  nearest[[i]] <- sf_points[which.min(st_distance(sf_points, fossil.sf[i,])), ]
  near.point <- rbind(near.point, unlist(nearest[[i]], recursive = FALSE)[1])
}

fossildata$near.point <- as.numeric(near.point)
fossildata <- merge(fossildata, sf_points, by.x = "near.point", by.y = "GlobalID")
fossildata

write.csv(fossildata, "outputs/fossildataGlobal.csv", row.names=TRUE)

fossilEMSmod <- fossilEMS + geom_rect(data = fossildata, 
                                      aes(xmin = as.numeric(mbc3) - 1,
                                          xmax = as.numeric(mbc3),
                                          ymin = as.numeric(sdbc3) - 1,
                                          ymax = as.numeric(sdbc3)), 
                                      colour = "#c44536", alpha = 0, size = 1)
fossilEMSmod

# Geographic space of fossil bin map
act_temp_fossil <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = temp), size = 2, pch = 16) +
  scale_color_gradientn(
    limits = c(-20, 30), 
    colors = color, 
    name = "Observed temperature (C)", 
    na.value = "transparent", 
    breaks = c(-20, -10, 0, 10, 20, 30), 
    labels = c(-20, -10, 0, 10, 20, 30)
  ) +
  labs(color = "LogPrecip") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  theme(
    legend.position = "top",  # ← move to top
    legend.direction = "horizontal",  # ← make horizontal
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", colour = "white")
  )

Time_A <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[1]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[1])
Time_BDE <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[2]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[2])
Time_C <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[3]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[3])

Time_A_points <- SpatialPointsDataFrame(Time_A[,2:3], Time_A)
sf_A_modern_points <- st_as_sf(Time_A_points)

Time_BDE_points <- SpatialPointsDataFrame(Time_BDE[,2:3], Time_BDE)
sf_BDE_modern_points <- st_as_sf(Time_BDE_points)

Time_C_points <- SpatialPointsDataFrame(Time_C[,2:3], Time_C)
sf_C_modern_points <- st_as_sf(Time_C_points)


Time_A_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_A_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  # scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Observed \nPrecipitation \n Log (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
  labs(title = "Time A") + theme(legend.position="none") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
# theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
#       legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))

Time_BDE_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_BDE_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  labs(title = "Times B D E") + theme(legend.position="none") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")

Time_C_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_C_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  labs(title = "Time C") + theme(legend.position="none")+
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")

combined_plot_temp <- act_temp_fossil + Time_A_fig + Time_BDE_fig+ Time_C_fig

ggsave("figures/combined_plot.png", plot = combined_plot_temp, width = 12, height = 8, dpi = 300)


# Africa models

points.af <- points.rest3 %>% filter(cont == "Africa") 

#Bin distributions 25x25
mRBLaf3 <- range(points.af$mean_RBL, na.rm=T)
sdRBLaf3 <- range(points.af$sd_RBL, na.rm = T)

#break points for mean & sd
mbrksaf3 <- seq(mRBLaf3[1]-0.001, mRBLaf3[2]+0.001, diff(mRBLaf3)/25)
sdbrksaf3 <- seq(sdRBLaf3[1]-0.001, sdRBLaf3[2]+0.001, diff(sdRBLaf3)/25)

#Assign bin codes
mbcaf3 <- .bincode(points.af$mean_RBL, breaks = mbrksaf3)
sdbcaf3 <- .bincode(points.af$sd_RBL, breaks = sdbrksaf3)

#Assign mean & sd bin code to each community
points.af$mbcaf3 <- mbcaf3
points.af$sdbcaf3 <- sdbcaf3


#Africa: Precipitation
mod <- list()
obj_af3 <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcaf3 == i & sdbcaf3 == j, na.rm = T) != 0) {
      dat <- points.af$temp[which(mbcaf3 == i & sdbcaf3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_af3[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_af3 <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_af3 <- setValues(r_af3, obj_af3)
rdf_af3 <- as.data.frame(r_af3, xy=TRUE)
names(rdf_af3) <- c("Mean", "SD", "Precip")

breaksm <- c((mbrks3[1]+(mRBL3[2]-mRBL3[1])/25),
             (mbrks3[13]+(mRBL3[2]-mRBL3[1])/25),
             (mbrks3[24]+(mRBL3[2]-mRBL3[1])/25))
breaksm <- trunc(breaksm*10^2)/10^2

breaksd <- c((sdbrks3[1]+(sdbrks3[2]-sdbrks3[1])/25), 
             (sdbrks3[13]+(sdbrks3[2]-sdbrks3[1])/25),
             (sdbrks3[24]+(sdbrks3[2]-sdbrks3[1])/25))
breaksd <- trunc(breaksd*10^2)/10^2

color <- c("#457b9d", "#fefae0", "#e63946")
precip_ecospace_af <-ggplot(data = rdf_af3) +
  geom_raster(mapping = aes(x = Mean, y = SD, fill = Precip)) +
  scale_fill_gradientn(colors = color, name = "Precipitation (loge(mm))", 
                       na.value = "transparent", 
                       breaks = waiver(), labels = waiver()) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.44), labels = breaksd) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.5), labels = breaksm) +
  coord_fixed() +
  theme_bw()
precip_ecospace_af


######Calculate maximum likelihood for all bins
modmax_af3 <- array(NA, dim = length(points.af[ , 1]))
mod <- list()

for (i in 1:length(points.af[ , 1])) {
  if(!(is.na(points.af$mbcaf3[i]) | is.na(points.af$sdbcaf3[i]))) {
    dat <- points.af$temp[which(points.af$mbcaf3 == points.af$mbcaf3[i] & points.af$sdbcaf3 == points.af$sdbcaf3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_af3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_af3
hist(modmax_af3)

#calculate anomaly
anom_precip_af3 <- points.af3$LogPrecip - modmax_af3
anom_precip_NL_af <- exp(anom_precip_af3) - 1
plot(anom_precip_af3, anom_precip_NL_af)

points.af3$precipanom <- as.numeric(anom_precip_af3)
summary(points.af3$precipanom, na.rm = T)
points.af3$precipest <- as.numeric(modmax_af3)
summary(anom_precip_NL_af, na.rm = T)

#Correlation of observed ~ expected 
precipcor_af <- lm(points.af3$precipest ~ points.af3$LogPrecip)
precipcor_af
summary(precipcor_af)
cor_precipest_af <- cor.test(points.af3$precipest, points.af3$LogPrecip, method = "pearson")
#cor=0.7679115 
cor_precipest_af















#display vectors with marginal boxplots
precip_ecospaceNA <- ggplot(data = rdf) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = NA, name = "Precipitation (mm)", na.value = "transparent", breaks = c(-3, 5, 13, 21, 29), labels = c(0, 21, 200, 755, 1950)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospaceNA

dfvectors <- data.frame(site = fossildata$Time,
                        fxmin = fossildata$fossilmbc1 - 1,
                        fxmax = fossildata$fossilmbc1,
                        fymin = fossildata$fossilsdbc1 - 1,
                        fymax = fossildata$fossilsdbc1,
                        mxmin = fossildata$mbc - 1,
                        mxmax = fossildata$mbc,
                        mymin = fossildata$sdbc - 1,
                        mymax = fossildata$sdbc)

vector <- data.frame("site" = rep(fossildata$Time, 2),
                     "xf" = rowMeans(dfvectors[, 2:3], na.rm = T),
                     "yf" = rowMeans(dfvectors[, 4:5], na.rm = T),
                     "xm" = rowMeans(dfvectors[, 6:7], na.rm = T),
                     "ym" = rowMeans(dfvectors[, 8:9], na.rm = T))

vector2 <- data.frame("Time" = rep(fossildata$Time, 2), "group" = c(rep("fossil", length(fossildata$Time)), rep("modern", length(fossildata$Time))), "x" = c(rowMeans(dfvectors[, 2:3], na.rm = T), rowMeans(dfvectors[, 6:7], na.rm = T)), "y" = c(rowMeans(dfvectors[, 4:5], na.rm = T), rowMeans(dfvectors[, 8:9], na.rm = T)))

vectorspace <- precip_ecospaceNA + 
  geom_rect(data = dfvectors, aes(xmin = as.numeric(mxmin),
                                  xmax = as.numeric(mxmax), 
                                  ymin = as.numeric(mymin), 
                                  ymax = as.numeric(mymax)), 
            colour = "#0868ac", alpha = 0, size = 1) +
  geom_segment(data = vector, aes(x = xf, y = yf, xend = xm, yend = ym), size = 0.8) +
  geom_point(data = vector2, aes(x = x, y = y, col = group), alpha = 0) +
  scale_color_manual(values = c("#8c510a", "#0868ac")) +
  theme(legend.position = "none")
vectorspace

library(ggExtra)

p <- ggMarginal(p = vectorspace, type = "histogram", 
                margins = "both", 
                size = 6,
                groupColour = TRUE,
                groupFill = TRUE)
p





##########
###################
#Africa ONLY
#Polygon shape files for continental shape
africa <- continents[continents$CONTINENT=="Africa",]
africa.sf <- st_as_sf(africa)

##convert sampling points to spatial points and sample species at each sampling locality
sp_af <- SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(africa))
sp_af <- raster ::intersect(sp_af, africa)
head(sp_af)
length(sp_af)

summary(sp_af$BIO12)
sd(sp_af$BIO12, na.rm=T)

summary(sp_af$BIO1)
sd(sp_af$BIO1, na.rm=T)

#Sample species at each sampling location.
# o_af <- over(sp_af, geography, returnList = T)
#save(o_af, file="o_af.Rdata")
load("o_af.Rdata")
head(o_af)
length(o_af)

points.af <- points.rest3 %>% filter(cont == "Africa") 
dim(points.af)

#Log transform precip
# points.af$LogPrecip <- log(points.af$BIO12+1)

#Summarize traits for community-level distribution by calculating mean and standard deviation
# points.af$mean_RBL <- unlist(lapply(o_af, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
# points.af$sd_RBL <- unlist(lapply(o_af, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))
# 
# ##calculate mammal richness at each point
# points.af$RBL_richness <- unlist(lapply(o_af, function(x) length(traits[x$binomial, "RBL"])))
# hist(points.af$RBL_richness)
# 
# #Calculate richness at each point with traits data
# points.af$count <- unlist(lapply(o_af, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))
# 
#Richness of just carnivoran species
# carn_rich_af <- array(NA, dim = length(o_af))
# for(i in 1:length(o_af)){
#   carn_rich_af[i]<-sum(o_af[[i]]$order_=="CARNIVORA")
# }
# hist(carn_rich_af)
# summary(carn_rich_af)
# points.af$carn_rich_af <-carn_rich_af
# summary(points.af$carn_rich_af)

#overall carnivoran richness (excluding zero-species communities)
# carn_no_zero_af<-points.af$carn_rich_af[points.af$carn_rich_af>0]
# summary(carn_no_zero_af)
# carn_three_plus_af<-points.af$carn_rich_af[points.af$carn_rich_af>2]
# #Global summary of carnivoran species richness in communities w/ at least 3 carnivoran species
# summary(carn_three_plus_af)
# sd(carn_three_plus_af)
# 
# points.af3<-points.af[points.af$count>2,]
# head(points.af3)
# dim(points.af3)
# 
# summary(points.af$mean_RBL)
# summary(points.af$sd_RBL)
points.af3 <- points.af
#Bin distributions 25x25
mRBLaf3 <- range(points.af3$mean_RBL, na.rm=T)
sdRBLaf3 <- range(points.af3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrksaf3 <- seq(mRBLaf3[1]-0.001, mRBLaf3[2]+0.001, diff(mRBLaf3)/25)
sdbrksaf3 <- seq(sdRBLaf3[1]-0.001, sdRBLaf3[2]+0.001, diff(sdRBLaf3)/25)

#Assign bin codes
mbcaf3 <- .bincode(points.af3$mean_RBL, breaks = mbrksaf3)
sdbcaf3 <- .bincode(points.af3$sd_RBL, breaks = sdbrksaf3)

#Assign mean & sd bin code to each community
points.af3$mbcaf3 <- mbcaf3
points.af3$sdbcaf3 <- sdbcaf3


#Africa: Precipitation
mod <- list()
obj_af3 <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcaf3 == i & sdbcaf3 == j, na.rm = T) != 0) {
      dat <- points.af3$LogPrecip[which(mbcaf3 == i & sdbcaf3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_af3[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_af3 <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_af3 <- setValues(r_af3, obj_af3)
rdf_af3 <- as.data.frame(r_af3, xy=TRUE)
names(rdf_af3) <- c("Mean", "SD", "Precip")

breaksm <- c((mbrks3[1]+(mRBL3[2]-mRBL3[1])/25),
             (mbrks3[13]+(mRBL3[2]-mRBL3[1])/25),
             (mbrks3[24]+(mRBL3[2]-mRBL3[1])/25))
breaksm <- trunc(breaksm*10^2)/10^2

breaksd <- c((sdbrks3[1]+(sdbrks3[2]-sdbrks3[1])/25), 
             (sdbrks3[13]+(sdbrks3[2]-sdbrks3[1])/25),
             (sdbrks3[24]+(sdbrks3[2]-sdbrks3[1])/25))
breaksd <- trunc(breaksd*10^2)/10^2

color <- c("#bc6c25", "#fefae0", "#606c38")
precip_ecospace_af <-ggplot(data = rdf_af3) +
  geom_raster(mapping = aes(x = Mean, y = SD, fill = Precip)) +
  scale_fill_gradientn(colors = rev(color), name = "Precipitation (loge(mm))", 
                       na.value = "transparent", 
                       breaks = waiver(), labels = waiver()) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.44), labels = breaksd) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), 
                     breaks = c(0.5, 13.5, 24.5), labels = breaksm) +
  coord_fixed() +
  theme_bw()
precip_ecospace_af


######Calculate maximum likelihood for all bins
modmax_af3 <- array(NA, dim = length(points.af3[ , 1]))
mod <- list()

for (i in 1:length(points.af3[ , 1])) {
  if(!(is.na(points.af3$mbcaf3[i]) | is.na(points.af3$sdbcaf3[i]))) {
    dat <- points.af3$LogPrecip[which(points.af3$mbcaf3 == points.af3$mbcaf3[i] & points.af3$sdbcaf3 == points.af3$sdbcaf3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_af3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_af3
hist(modmax_af3)

#calculate anomaly
anom_precip_af3 <- points.af3$LogPrecip - modmax_af3
anom_precip_NL_af <- exp(anom_precip_af3) - 1
plot(anom_precip_af3, anom_precip_NL_af)

points.af3$precipanom <- as.numeric(anom_precip_af3)
summary(points.af3$precipanom, na.rm = T)
points.af3$precipest <- as.numeric(modmax_af3)
summary(anom_precip_NL_af, na.rm = T)

#Correlation of observed ~ expected 
precipcor_af <- lm(points.af3$precipest ~ points.af3$LogPrecip)
precipcor_af
summary(precipcor_af)
cor_precipest_af <- cor.test(points.af3$precipest, points.af3$LogPrecip, method = "pearson")
#cor=0.7679115 
cor_precipest_af




#######Africa fossil
####Fossil analysis (precip)
fossil <- read.csv("inputs/Koobi_fora_full.csv", header = T)
dim(fossil)
fossil <- fossil[1:100, 1:9]

continents1 <- continents
continents1 <- st_crs(4326)

#Fossil map
fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))
fossmap <- act_precip + geom_sf(data = fossil.sf, col = "black", size = 2, pch = 17)
fossmap

#aggregating fossil data
fossilmean <- aggregate(fossil[ , 9], list(fossil$Time, fossil$Lat, fossil$Long), mean, na.rm = TRUE)
fossilsd <- aggregate(fossil[ , 9], list(fossil$Time), sd, na.rm = TRUE)
colnames(fossilmean) <- c("Time", "Lat", "Long", "Mean")
fossilmean <- merge(fossilmean, fossilsd, by.x = "Time", by.y = "Group.1")
colnames(fossilmean)[5] <- "SD"
fossildata <- fossilmean
head(fossildata)

fossildata$fossilmbc1 <- .bincode(fossildata$Mean, breaks = mbrks3)
fossildata$fossilsdbc1 <- .bincode(fossildata$SD, breaks = sdbrks3)

###
fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))

#calculate fossil precip estimates with 5% limits
fossilmodmax <- list() 
fossilmod <- list()
for (i in 1:length(fossildata[ , 1])) {
  if(!(is.na(fossildata$fossilmbc1[i]) | is.na(fossildata$fossilsdbc1[i]))) {
    if(length(which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])) < 1) {
      fossilmod[[i]] <- NA 
    } else {
      dat <- points.af3$LogPrecip[which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])]
      fossilmod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      fossilmodmax$precipest[i] <- fossilmod[[i]]$x[which.max(fossilmod[[i]]$y)]
      maxlength <- length(which(fossilmod[[i]]$x < fossilmodmax$precipest[i]))
      bound.05 <- (length(fossilmod[[i]]$x)) * 0.05
      fossilmodmax$minlimit[i] <- fossilmod[[i]]$x[maxlength - bound.05]
      fossilmodmax$maxlimit[i] <- fossilmod[[i]]$x[maxlength + bound.05]
    }
  }
}

fossildata$precipest <- c(fossilmodmax$precipest, NA)
fossildata$minlimit <- c(fossilmodmax$minlimit, NA)
fossildata$maxlimit <- c(fossilmodmax$maxlimit, NA)

fossildata$precipestUN <- c(exp(fossilmodmax$precipest)-(1), NA)
fossildata$minlimitUN <- c(exp(fossilmodmax$minlimit)-(1), NA)
fossildata$maxlimitUN <- c(exp(fossilmodmax$maxlimit)-(1), NA)

#plot fossil sites in ecometric space
df1 <- data.frame(xmin = fossildata$fossilmbc1 - 1,
                  xmax = fossildata$fossilmbc1,
                  ymin = fossildata$fossilsdbc1 - 1,
                  ymax = fossildata$fossilsdbc1)

text <- data_frame(point = c(1:6),
                   x = fossildata$fossilmbc1 -.5,
                   y = fossildata$fossilsdbc1 + .5)

text$x[c(2)] <- text$x[c(2)] - 1
text$x[c(5)] <- text$x[c(5)] +1

fossilEMS <- precip_ecospace_af + 
  geom_rect(data = df1, aes(xmin = as.numeric(xmin), 
                            xmax = as.numeric(xmax), ymin = as.numeric(ymin),
                            ymax = as.numeric(ymax)), colour = "black", 
            alpha = 0, linewidth = 0.75) + geom_text(text, mapping = aes(x = x,
                                                                         y = y,
                                                                         label = point))
fossilEMS
#select sampling point and observed precipitation nearest to fossil site
fossil.sf <- st_as_sf(fossildata, coords = c(3, 2), crs = crs(continents1))

nearest <- list()
near.point <- list()
for(i in 1:dim(fossildata)[1]) {
  nearest[[i]] <- sf_points[which.min(st_distance(sf_points, fossil.sf[i,])), ]
  near.point <- rbind(near.point, unlist(nearest[[i]], recursive = FALSE)[1])
}

fossildata$near.point <- as.numeric(near.point)
fossildata <- merge(fossildata, sf_points, by.x = "near.point", by.y = "GlobalID")
fossildata

write.csv(fossildata, "/Users/leila/Desktop/fossildata.csv", row.names=TRUE)

fossilEMSmod <- fossilEMS + geom_rect(data = fossildata, 
                                      aes(xmin = as.numeric(mbc3) - 1,
                                          xmax = as.numeric(mbc3),
                                          ymin = as.numeric(sdbc3) - 1,
                                          ymax = as.numeric(sdbc3)),
                                      colour = "#c44536", alpha = 0, size = 1)
fossilEMSmod

#display vectors with marginal boxplots
precip_ecospaceNA <- ggplot(data = rdf) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = NA, name = "Precipitation (mm)", na.value = "transparent", breaks = c(-3, 5, 13, 21, 29), labels = c(0, 21, 200, 755, 1950)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospaceNA

dfvectors <- data.frame(site = fossildata$Time,
                        fxmin = fossildata$fossilmbc1 - 1,
                        fxmax = fossildata$fossilmbc1,
                        fymin = fossildata$fossilsdbc1 - 1,
                        fymax = fossildata$fossilsdbc1,
                        mxmin = fossildata$mbc - 1,
                        mxmax = fossildata$mbc,
                        mymin = fossildata$sdbc - 1,
                        mymax = fossildata$sdbc)

vector <- data.frame("site" = rep(fossildata$Time, 2),
                     "xf" = rowMeans(dfvectors[, 2:3], na.rm = T),
                     "yf" = rowMeans(dfvectors[, 4:5], na.rm = T),
                     "xm" = rowMeans(dfvectors[, 6:7], na.rm = T),
                     "ym" = rowMeans(dfvectors[, 8:9], na.rm = T))

vector2 <- data.frame("Time" = rep(fossildata$Time, 2), "group" = c(rep("fossil", length(fossildata$Time)), rep("modern", length(fossildata$Time))), "x" = c(rowMeans(dfvectors[, 2:3], na.rm = T), rowMeans(dfvectors[, 6:7], na.rm = T)), "y" = c(rowMeans(dfvectors[, 4:5], na.rm = T), rowMeans(dfvectors[, 8:9], na.rm = T)))

vectorspace <- precip_ecospaceNA + 
  geom_rect(data = dfvectors, aes(xmin = as.numeric(mxmin),
                                  xmax = as.numeric(mxmax), 
                                  ymin = as.numeric(mymin), 
                                  ymax = as.numeric(mymax)), 
            colour = "#0868ac", alpha = 0, size = 1) +
  geom_segment(data = vector, aes(x = xf, y = yf, xend = xm, yend = ym), size = 0.8) +
  geom_point(data = vector2, aes(x = x, y = y, col = group), alpha = 0) +
  scale_color_manual(values = c("#8c510a", "#0868ac")) +
  theme(legend.position = "none")
vectorspace

library(ggExtra)

p <- ggMarginal(p = vectorspace, type = "histogram", 
                margins = "both", 
                size = 6,
                groupColour = TRUE,
                groupFill = TRUE)
p