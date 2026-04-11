#############RBL Analyses
#Worldwide & N. America
#Temp & Precip
#Sensitivity analysis
#Seven paleo sites
#Updated 12/8/2023

#setwd("/Users/leila/Desktop/Carn_RBL")
#setwd("C:/Users/rachel.short/OneDrive - South Dakota State University - SDSU/Colleague Work/Leila/RBLcode/Code _ data")

#install.packages(c("sp", "raster", "sf", "ggplot2"))
library(sp)
library(raster)
library(sf)
library(ggplot2)

#Polygon shape files for continental shape
continents <- shapefile("continents.shp")
continents.sf <- st_as_sf(continents)

#All global environmental data at 50km equidistant points
points <- read.csv("S_L_allpointsdata2.csv")
head(points)
dim(points)

#Convert sampling points to spatial points and sample species at each sampling locality
sp <- SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(continents))
sp <- raster::intersect(sp, continents)
head(sp)
length(sp)

#Save continent to points
points$cont<- sp@data$CONTINENT

##Load RBL trait data & format file
traits <- read.csv("RBL_global_Aug2023.csv", header = T, na= ".", stringsAsFactors = FALSE)
traits$TaxonName <- gsub("_", " ", traits$TaxonName)
rownames(traits) <- as.character(traits$TaxonName)
colnames(traits) <- c("binomial", "RBL", "IUCN", "Percent_vert", "Percent_nonvert", "Diet", "Diet_vert", "Body_size", "Diet_breadth", "Home_Range", "Geographic_Range", "Publication", "Sample_size", "Notes", "x")
names(traits)
table(traits$Diet)
traits$Diet <- as.factor(traits$Diet)
traits$Percent_vert <- as.factor(traits$Percent_vert)

#########################################
#this is an ANOVA that shows that there is a functional signal with diet category
model_1 <- summary(lm(traits$RBL~traits$Diet))
model_1
model_1$r.squared
#########################################

#########################################
#this narrows down the signal to percent vertebrate prey in the diet
model_2 <- summary(lm(traits$RBL~traits$Percent_vert))
model_2
model_2$r.squared
#########################################

se_rbl <- tapply(traits$RBL, traits$Diet, function(x) sd(x) / sqrt(length(x)))
sd_rbl <- tapply(traits$RBL, traits$Diet, sd)
mean_rbl <- tapply(traits$RBL, traits$Diet, mean)
diet <- levels(traits$Diet)
traits$Diet_se <- NA
for(i in 1:length(unique(traits$Diet))){
  traits$Diet_se[traits$Diet == unique(traits$Diet)[i]] <- se_rbl[i]
}
  
ggplot() +
  geom_bar(aes(x = diet, y = mean_rbl), stat = "identity", fill = "skyblue", alpha = 0.7) + 
  geom_errorbar(aes(x = diet, ymin = mean_rbl - se_rbl, ymax = mean_rbl + se_rbl), width = 0.4, colour = "orange", alpha = 0.9, size = 1.3)

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

#Rank order plot color-coded by diet type
traits$Diet <- factor(traits$Diet, levels=c('vert', 'omn', 'invert', 'plant', 'fruit'))
rankorderspecies1<-ggplot(traits, aes(x = binomial, y = RBL, col = Diet)) +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T), linetype = "dashed", color = "black") +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T) - sd(traits$RBL), linetype = "solid", color = "black") +
  geom_hline(yintercept = mean(traits$RBL, na.rm = T) + sd(traits$RBL), linetype = "solid", color = "black") +
  geom_segment(aes(x = binomial, xend = binomial, y = 0, yend = RBL), lwd=0.3) +
  geom_point(aes(x = binomial, y = RBL), pch = 16, cex = 2) +
  labs(x = "Species", y = "RBL") +
  scale_color_manual(values=c("red4", "orangered", "gold", "yellowgreen", "mediumturquoise"), labels=c("Carnivore", "Omnivore", "Invertivore", "Herbivore", "Frugivore")) +
  scale_x_discrete(breaks = levels(traits$binomial)[c(T, rep(F, 3))]) +
  theme(axis.line = element_line(color = "black"), 
        panel.background = element_rect(color = "black", fill = NA), 
        panel.grid = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", angle = 90, size = 12, hjust = 1, vjust = 0.25, face="italic"),
        axis.title = element_text(color = "black", size = 14)) +
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12))
rankorderspecies1

########################
#RS to double check diet categories
########################

#geography <- shapefile("TERRESTRIAL_MAMMALS.shp")
#save(geography, file="geography.Rdata")
load("geography.Rdata")
head(geography)

#Community of mammals at each sampling location
#o <- over(sp, geography, returnList = T)
#save(o, file="o.Rdata")
load("o.Rdata")
length(o)

points.iucn <- points

#Log precip
points.iucn$LogPrecip <- log(points$BIO12+1)
#points.iucn$TempCube <- ((points$BIO1)^(3))

#Summarize traits for community-level distribution: calculating mean & standard deviation
points.iucn$mean_RBL <- unlist(lapply(o, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
hist(points.iucn$mean_RBL)
points.iucn$sd_RBL <- unlist(lapply(o, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))
hist(points.iucn$sd_RBL)

##calculate carn richness at each point
carn_rich<-array(NA, dim = length(o))
for(i in 1:length(o)){
  carn_rich[i]<-sum(o[[i]]$order_=="CARNIVORA")
}
summary(carn_rich)
points.iucn$carn_rich <- carn_rich

#Calculate richness at each point with traits data
points.iucn$RBL_richness <- unlist(lapply(o, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))
hist(points.iucn$RBL_richness)

####Global 
#Min 3 species per community
points.rest3 <- points.iucn[points.iucn$carn_rich > 2, ]
dim(points.rest3)

summary(points.rest3$mean_RBL)
sd(points.rest3$mean_RBL)

#Bin distributions 25x25
#***place to refine # of bins
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

#this ecometric is not linear relationship, leaving this code here for now, but need to remove
# plot(points.rest3$mean_RBL, points.rest3$LogPrecip)
# reg_precip <- lm(LogPrecip ~ mean_RBL + sd_RBL, data = points.rest3)
# summary(reg_precip)
# cor_precip <- cor.test(points.rest3$LogPrecip, points.rest3$mean_RBL, method = "pearson")
# cor_precip

# #precipitation and sd RBL, not linear relationship
# plot(points.rest3$sd_RBL, points.rest3$LogPrecip)
# reg_precip2 <- lm(LogPrecip ~ sd_RBL, data = points.rest3)
# reg_precip2
# summary(reg_precip2)
# cor_precip2 <- cor.test(points.rest3$LogPrecip, points.rest3$sd_RBL, method = "pearson")
# cor_precip2

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

##RS: WHY SO MANY COLORS? COLORRAMPPALETTE WILL MAKE THE GRADIENT.
##RS: YOU GAVE 7 COLORS AND TOLD IT TO MAKE A GRADIENT OF 8 ???
color <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r) - minValue(r)))
precip_ecospace <-ggplot(data = rdf) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color), name = "Precipitation (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "St. Dev. RBL", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean RBL", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
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

##histogram of max like estimate
hist(modmax3)

#calculate precipitation anomaly
anom_precip <- points.rest3$LogPrecip - modmax3
anom_precip_NL <- exp(anom_precip) - 1 #removing the log transform

points.rest3$precipanom <- as.numeric(anom_precip)
points.rest3$precipanomNL <- as.numeric(anom_precip_NL)
summary(points.rest3$precipanom, na.rm = T)
summary(points.rest3$precipanomNL, na.rm = T)
points.rest3$precipest <- as.numeric(modmax3)
summary(points.rest3$precipest, na.rm = T)

#Correlation of observed ~ expected precipitation
precipcor <- lm(points.rest3$precipest ~ points.rest3$LogPrecip)
summary(precipcor)
cor_precipest <- cor.test(points.rest3$precipest, points.rest3$LogPrecip, method = "pearson")
#cor=0.7007268, new 0.7147735
cor_precipest

###########################
#########Randomization
###########################
cor_precipest_rand <- array(NA, dim = 100)
for(i in 1:100){
  cor_precipest_rand[i] <- cor(points.rest3$precipest[sample(1:nrow(points.rest3))], points.rest3$LogPrecip, use = "pairwise.complete.obs")
}
#cor=0.7007268, new 0.7147735
cor_precipest
hist(cor_precipest_rand, xlab = "Pearson's Coefficient - Randomization Test", main = paste("Original Pearson's Coefficient:", round(cor_precipest$estimate, digits = 2)))
#the histogram above shows that the random correlation values are very low and near zero and our value is far above any of those.

###########################
#the other randomization, only ran it for 21 times, because it takes a long time to run. 21 times are enough to see the distribution of correlation coefficients
###########################
# iteration<-100
# modmax_rand <- array(NA, dim = c(iteration,length(points.rest3[ , 1])))
# mod_rand <- list()
# rand_cor <- array(NA, dim = iteration)
# for(j in 22:iteration){
#   for (i in 1:length(points.rest3[ , 1])) {
#     if(!(is.na(points.rest3$mbc3[i]) | is.na(points.rest3$sdbc3[i]))) {
#       community_length<-length(which(points.rest3$mbc3 == 
#               points.rest3$mbc3[i] & points.rest3$sdbc3 == points.rest3$sdbc3[i]))
#       dat <- points.rest3$LogPrecip[sample(1:length(points.rest3$LogPrecip), size = community_length)]
#       mod_rand[[i]] <- density(dat[!is.na(dat)], bw = 1)
#       modmax_rand[j, i] <- mod_rand[[i]]$x[which.max(mod_rand[[i]]$y)]
#     }
#   }
#   #Save out correlation of observed ~ expected precipitation
#   rand_cor[j] <- cor(modmax_rand[j,], points.rest3$LogPrecip, use = "pairwise.complete.obs")
# }
# 
# save(rand_cor, file="rand_cor.Rdata")

load("rand_cor.Rdata")

rand_cor_hist <- hist(rand_cor)
hist(rand_cor, xlab = "Pearson's Coefficient - Randomization Test", main = paste("Original Pearson's Coefficient:", round(cor_precipest$estimate, digits = 2)))
#the histogram above shows that the random correlation values are very low and near zero and our value is far above any of those.


#################

#Maps
sp_points<-SpatialPointsDataFrame(points.rest3[,2:3],points.rest3)
sf_points <- st_as_sf(sp_points)

#Maps of RBL
color2 <- colorRampPalette(c("red4","orangered", "gold", "yellowgreen", "mediumturquoise", "aquamarine"))(round(maxValue(r) - minValue(r)))
color_RBL <- colorRampPalette(c("darkorchid4", "darkorchid4","darkorchid","violetred", "deeppink", "lightcoral","lightsalmon", "navajowhite","#ffffcc"))(round(maxValue(r) - minValue(r)))
#Mean RBL map
rbl_map<-ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = mean_RBL), size = 1.5, pch = 16)  +
  scale_color_gradientn(limits=as.numeric(c(0.55, 0.95)), colors = rev(color_RBL), name = "RBL Mean", na.value = "transparent", breaks = c(0.55, 0.75, 0.95), labels = c(0.55, 0.75, 0.95)) +
  labs(color = "mean_RBL") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") 
rbl_map

#RBL variation map
rbl_s_map<-ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = sd_RBL), size = 2, pch = 16)  +
  scale_color_gradientn(limits=as.numeric(c(0, 0.30)), colors = rev(color_RBL), name = "RBL \nStandard deviation", na.value = "transparent", breaks = c(0, 0.1, 0.2, 0.3), labels = c(0, 0.1, 0.2, 0.3)) +
  labs(color = "sd_RBL") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") 
rbl_s_map

##Actual precip map
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
##I SUGGEST LETTING R DETERMINE THE LIMITS AND BREAKS. IF YOU CHANGE ANYTHING, YOU HAVE TO REMEMBER TO
##COME BACK HERE AND CHANGE YOUR MANUALLY ENTERED NUMBERS.
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
##IF YOU WANT THE SCALE CENTERED ON 0, WE CAN MAKE THAT HAPPEN.

##THIS IS WHERE RACHEL STOPPED

##########
####Fossil analysis (precip)
fossil <- read.csv("fossil_RBL.csv", header = T)
dim(fossil)
fossil <- fossil[1:71, 1:13]

continents1 <- continents
continents1 <- st_crs(4326)

#Fossil map
fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))
fossmap <- act_precip + geom_sf(data = fossil.sf, col = "black", size = 2, pch = 17)
fossmap

#aggregating fossil data
fossilmean <- aggregate(fossil[ , 11], list(fossil$Site, fossil$Lat, fossil$Long), mean, na.rm = TRUE)
fossilsd <- aggregate(fossil[ , 11], list(fossil$Site), sd, na.rm = TRUE)
colnames(fossilmean) <- c("Site", "Lat", "Long", "Mean")
fossilmean <- merge(fossilmean, fossilsd, by.x = "Site", by.y = "Group.1")
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

fossildata$precipest <- fossilmodmax$precipest
fossildata$minlimit <- fossilmodmax$minlimit
fossildata$maxlimit <- fossilmodmax$maxlimit

fossildata$precipestUN <- exp(fossilmodmax$precipest)-(1)
fossildata$minlimitUN <- exp(fossilmodmax$minlimit)-(1)
fossildata$maxlimitUN <- exp(fossilmodmax$maxlimit)-(1)

#plot fossil sites in ecometric space
df1 <- data.frame(xmin = fossildata$fossilmbc1 - 1,
                  xmax = fossildata$fossilmbc1,
                  ymin = fossildata$fossilsdbc1 - 1,
                  ymax = fossildata$fossilsdbc1)
fossilEMS <- precip_ecospace + geom_rect(data = df1, aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax)),colour = "black", alpha = 0, linewidth = 0.75)
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

fossilEMSmod <- fossilEMS + geom_rect(data = fossildata, aes(xmin = as.numeric(mbc3)-1,xmax = as.numeric(mbc3), ymin = as.numeric(sdbc3)-1,ymax = as.numeric(sdbc3)), colour = "red", alpha = 0, linewidth = 1)
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

dfvectors <- data.frame(site = fossildata$Site,
                        fxmin = fossildata$fossilmbc1 - 1,
                        fxmax = fossildata$fossilmbc1,
                        fymin = fossildata$fossilsdbc1 - 1,
                        fymax = fossildata$fossilsdbc1,
                        mxmin = fossildata$mbc - 1,
                        mxmax = fossildata$mbc,
                        mymin = fossildata$sdbc - 1,
                        mymax = fossildata$sdbc)

vector <- data.frame("site" = rep(fossildata$Site, 2),
                     "xf" = rowMeans(dfvectors[, 2:3], na.rm = T),
                     "yf" = rowMeans(dfvectors[, 4:5], na.rm = T),
                     "xm" = rowMeans(dfvectors[, 6:7], na.rm = T),
                     "ym" = rowMeans(dfvectors[, 8:9], na.rm = T))

vector2 <- data.frame("site" = rep(fossildata$Site, 2), "group" = c(rep("fossil", length(fossildata$Site)), rep("modern", length(fossildata$Site))), "x" = c(rowMeans(dfvectors[, 2:3], na.rm = T), rowMeans(dfvectors[, 6:7], na.rm = T)), "y" = c(rowMeans(dfvectors[, 4:5], na.rm = T), rowMeans(dfvectors[, 8:9], na.rm = T)))

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


####Sensitivity Code (precip)
library(reshape2)
library(ggplot2)
library(ggpubr)

#number of bins
bins <- 25

#function for returning bin codes
bincodes <- function(traits) {
  trait_range <- range(traits, na.rm = T) #take the range of the means
  trait_range[1] <- trait_range[1]-0.000001 #add a little amount on each side of the range
  trait_range[2] <- trait_range[2]+0.000001
  brks <- seq(trait_range[1], trait_range[2], diff(trait_range)/bins) #split the range into the number of specified equal parts
  bc <- .bincode(traits, breaks = brks) #assign a bin code to each of the bins
  return(bc)
}

#function to retrieve bin breaks
binbreaks <- function(traits) {
  trait_range <- range(traits, na.rm = T) #take the range of the means
  trait_range[1] <- trait_range[1]-0.000001 #add a little amount on each side of the range
  trait_range[2] <- trait_range[2]+0.000001
  brks <- seq(trait_range[1], trait_range[2], diff(trait_range)/bins) #split the range into equal parts
  return(brks)
}

#define sample sizes
samplesizes <- seq(100, 10000, 1000) 

#define number of bootstrap iterations for each
sens_it <- 20

#define split of testing and training data
test_split <- 0.2

#set vars to collect data
#sensitivity is within the model and testing is transferability
anomp_sens <- array(NA, dim = c(length(samplesizes),sens_it))
anomp_test <- array(NA, dim = c(length(samplesizes),sens_it))
cor_sens <- array(NA, dim = c(length(samplesizes),sens_it))
cor_test <- array(NA, dim = c(length(samplesizes),sens_it))

#loop to downsample and split dataset into training and testing parts (80/20)
# for(samp_sizes in 8:length(samplesizes)){
#   testsample <- sample(points.rest3$GlobalID, size = samplesizes[samp_sizes], replace = F)
#   downdata <- points.rest3[points.rest3$GlobalID %in% testsample,]
#   
#   for (sens in 1:sens_it){
#     training <- sample(x = 1:length(downdata[,1]), size = length(downdata[,1])*(1 - test_split))
#     testing <- (1:length(downdata[,1]))[!(1:length(downdata[,1]) %in% training)]
#     training_data <- downdata[training,]
#     testing_data <- downdata[testing,]
#     
#     ## Ecometric Spaces ##
#     #Group the trait distributions into bins x bins trait bins & calculate the range
#     training_data$mbc.tr <- bincodes(training_data$mean_RBL) #ecometric space calculations, means
#     training_data$sdbc.tr <- bincodes(training_data$sd_RBL)  #ecometric space calculations, standard deviations
#     
#     #Calculate difference between observed and trait-predicted environment
#     maxlike <- array(NA, dim = length(training_data[ , 1])) 
#     obs_env <- array(NA, dim = length(training_data[ , 1])) 
#     anomp <- array(NA, dim = length(training_data[ , 1]))
#     
#     mod <- list()
#     
#     for (i in 1:length(training_data[ , 1])) {
#       if(!(is.na(training_data$mbc.tr[i]) | is.na(training_data$sdbc.tr[i]) | is.na(training_data$LogPrecip[i]))) {
#         if(length(which(training_data$mbc.tr == training_data$mbc.tr[i] & training_data$sdbc.tr == training_data$sdbc.tr[i])) < 1) {
#           mod[[i]] <- NA 
#         } else {
#           dat <- training_data$LogPrecip[which((training_data$mbc.tr == training_data$mbc.tr[i] & training_data$sdbc.tr == training_data$sdbc.tr[i]))]
#           mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
#           maxlike[i] <- mod[[i]]$x[which.max(mod[[i]]$y)] 
#           obs_env[i] <- training_data$LogPrecip[i]
#           anomp[i] <- obs_env[i] - maxlike[i]
#         }
#       }
#     }
#     
#     #hold anomaly for each testing training split
#     anomp <- anomp[!is.infinite(anomp)]
#     anomp <- anomp[!is.na(anomp)]
#     anomp_sens[samp_sizes,sens]<- mean(abs(anomp), na.rm = T)
#     cor_sens[samp_sizes,sens] <- cor(maxlike, obs_env, use = "pairwise.complete.obs")
#     
#     #testing data, get bin breaks from training data
#     mbrks <- binbreaks(training_data$mean_RBL)
#     sdbrks <- binbreaks(training_data$sd_RBL)
#     
#     #get bin codes for testing data
#     testing_data$mbc.te <- .bincode(testing_data$mean_RBL, breaks = mbrks)
#     testing_data$sdbc.te <-.bincode(testing_data$sd_RBL, breaks = sdbrks)
#     
#     #reset the values
#     maxlike <- array(NA, dim = length(training_data[ , 1])) 
#     obs_env <- array(NA, dim = length(training_data[ , 1])) 
#     anomp <- array(NA, dim = length(training_data[ , 1]))
#     
#     mod <- list()
#     
#     for (j in 1:length(testing_data[ , 1])) {
#       if(!(is.na(testing_data$mbc.te[j]) | is.na(testing_data$sdbc.te[j]) | is.na(testing_data$LogPrecip[j]))) {
#         if(length(which(training_data$mbc.tr == testing_data$mbc.te[j] & training_data$sdbc.tr == testing_data$sdbc.te[j])) < 1) {
#           mod[[j]] <- NA 
#         } else {
#           dat <- training_data$LogPrecip[which(training_data$mbc.tr == testing_data$mbc.te[j] & training_data$sdbc.tr == testing_data$sdbc.te[j])]
#           mod[[j]] <- density(dat[!is.na(dat)], bw = 1) 
#           maxlike[j] <- mod[[j]]$x[which.max(mod[[j]]$y)] 
#           obs_env[j] <- testing_data$LogPrecip[j]
#           anomp[j] <- obs_env[j] - maxlike[j]
#         }
#       }
#     }
#     anomp <- anomp[!is.infinite(anomp)]
#     anomp <- anomp[!is.na(anomp)]
#     anomp_test[samp_sizes,sens] <- mean(abs(anomp), na.rm = T)
#     cor_test[samp_sizes,sens] <- cor(maxlike, obs_env, use = "pairwise.complete.obs")
#   }
# }
# 
# training_cor <- data.frame(samplesizes = samplesizes, cor_sens)
# training_cor <- melt(training_cor, id.vars = 1)
# 
# testing_cor <- data.frame(samplesizes = samplesizes, cor_test)
# testing_cor <- melt(testing_cor, id.vars = 1)
# 
# save(training_cor, file = "training_cor.Rdata")
# save(testing_cor, file = "testing_cor.Rdata")

load("training_cor.Rdata")
load("testing_cor.Rdata")

testing_plot <- ggplot(testing_cor, aes(samplesizes, testing_cor[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  ylim(0.4, 1) +
  xlab("Communities") +
  ylab("Testing - Cor") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic() 

training_plot <- ggplot(training_cor, aes(samplesizes, training_cor[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  ylim(0.4, 1) +
  xlab("Communities") +
  ylab("Training - Cor") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic()

ggarrange(training_plot, testing_plot,
          labels = c("A", "B"), ncol = 1, nrow = 2)

# training_anom <- data.frame(samplesizes = samplesizes, anomp_sens)
# training_anom <- melt(training_anom, id.vars = 1)
# 
# testing_anom <- data.frame(samplesizes = samplesizes, anomp_test)
# testing_anom <- melt(testing_anom, id.vars = 1)
# 
# save(training_anom, file = "training_anom.Rdata")
# save(testing_anom, file = "testing_anom.Rdata")

load("training_anom.Rdata")
load("testing_anom.Rdata")

testing_plot <- ggplot(testing_anom, aes(samplesizes, testing_anom[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  ylim(0, 1) +
  xlab("Communities") +
  ylab("Testing - Anom") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic() 

training_plot <- ggplot(training_anom, aes(samplesizes, training_anom[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  ylim(0, 1) +
  xlab("Communities") +
  ylab("Training - Anom") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic()

ggarrange(training_plot, testing_plot, labels = c("A", "B"), ncol = 1, nrow = 2)





##########################Temperature, min species=3
mod <- list()
obj_t <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbc3 == i & sdbc3 == j, na.rm = T) != 0) {
      dat <- points.rest3$BIO1[which(mbc3 == i & sdbc3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_t[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}


##make raster for body mass and temperature
r_t <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_t <- setValues(r_t, obj_t)
rdf_t <- as.data.frame(r_t, xy=TRUE)
names(rdf_t) <- c("Mean", "SD", "Temp")
head(rdf_t)

color_t <- colorRampPalette(c("darkorchid4", "darkorchid4","darkorchid","violetred", "deeppink", "lightcoral","lightsalmon"))(round(maxValue(r_t) - minValue(r_t)))
temp_ecospace <-ggplot(data = rdf_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp/10)) +
  scale_fill_gradientn(colors = rev(color_t), name = "Temperature (C)", na.value = "transparent") +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
temp_ecospace

######Calculate maximum likelihood for all bins
modmax3_t <- array(NA, dim = length(points.rest3[ , 1]))
mod_temp <- list()

for (i in 1:length(points.rest3[ , 1])) {
  if(!(is.na(points.rest3$mbc3[i]) | is.na(points.rest3$sdbc3[i]))) {
    dat <- points.rest3$BIO1[which(points.rest3$mbc3 == points.rest3$mbc3[i] & points.rest3$sdbc3 == points.rest3$sdbc3[i])]
    mod_temp[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax3_t[i] <- mod_temp[[i]]$x[which.max(mod_temp[[i]]$y)]
  }
}
hist(modmax3_t)

#calculate temperature anomaly
anom_temp <- points.rest3$BIO1/10 - modmax3_t/10
summary(anom_temp)
hist(anom_temp)

# #Inverse of cubed temp (cube root)
# #R has a bit of trouble with cube root for datasets w/ negative numbers, so a new function is needed 
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }

# #testing the function, the TempCube data should be the same as BIO1 if the function works
# uncube<-Math.cbrt(points.rest3$TempCube)
# summary(uncube)
# summary(points.rest3$BIO1)

points.rest3$tempanom <- as.numeric(anom_temp)
summary(points.rest3$tempanom, na.rm = T)
points.rest3$tempest <- as.numeric(modmax3_t/10)
summary(points.rest3$tempest, na.rm = T)

#Correlation of observed ~ expected temperature
tempcor <- lm(points.rest3$tempest ~ points.rest3$BIO1)
summary(tempcor)
cor_tempest <- cor.test(points.rest3$tempest, points.rest3$BIO1, method = "pearson")
#cor=0.5430605, 0.6772141
cor_tempest

#########Randomization
iteration<-100
modmax_rand_temp <- array(NA, dim = c(iteration,length(points.rest3[ , 1])))
mod_rand_temp <- list()
rand_cor_temp <- array(NA, dim = iteration)
for(j in 23:iteration){
  for (i in 1:length(points.rest3[ , 1])) {
    if(!(is.na(points.rest3$mbc3[i]) | is.na(points.rest3$sdbc3[i]))) {
      community_length<-length(which(points.rest3$mbc3 ==
              points.rest3$mbc3[i] & points.rest3$sdbc3 == points.rest3$sdbc3[i]))
      dat <- points.rest3$BIO1[sample(1:length(points.rest3$BIO1), size = community_length)]
      mod_rand_temp[[i]] <- density(dat[!is.na(dat)], bw = 1)
      modmax_rand_temp[j, i] <- mod_rand_temp[[i]]$x[which.max(mod_rand_temp[[i]]$y)]
    }
  }
  #Save out correlation of observed ~ expected precipitation
  rand_cor_temp[j] <- cor(modmax_rand_temp[j,], points.rest3$BIO1, use = "pairwise.complete.obs")
}

save(rand_cor_temp, file="rand_cor_temp.Rdata")

load("rand_cor_temp.Rdata")

rand_cor_hist_temp <- hist(rand_cor_temp)
hist(rand_cor_temp, xlab = "Pearson's Coefficient - Randomization Test", main = paste("Original Pearson's Coefficient:", round(cor_tempest$estimate, digits = 2)))

#the other randomization
cor_tempest_rand <- array(NA, dim = 100)
for(i in 1:100){
  cor_tempest_rand[i] <- cor(points.rest3$tempest[sample(1:nrow(points.rest3))], points.rest3$BIO1, use = "pairwise.complete.obs")
}
hist(cor_tempest_rand, xlab = "Pearson's Coefficient - Randomization Test", main = paste("Original Pearson's Coefficient:", round(cor_tempest$estimate, digits = 2)))

###########

#Maps
sp_points<-SpatialPointsDataFrame(points.rest3[,2:3],points.rest3)
sf_points <- st_as_sf(sp_points)

##Actual temp map
act_temp<-ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void()  +
  geom_sf(data = sf_points, aes(color = BIO1/10), size = 0.5, pch = 16) +
  scale_color_gradientn(limits=as.numeric(c(-30, 30)), colors = rev(color_t), name = "Observed \nTemperature", na.value = "transparent") +
  labs(color = "BIO1") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
act_temp

##Map of predicted temp
est_temp <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() + 
  geom_sf(data = sf_points, aes(color = tempest), size = 0.5, pch = 16) +
  scale_color_gradientn(limits = as.numeric(c(-30, 30)), colors = rev(color_t), name = "Estimated \nTemp (C)", na.value = "transparent") +
  labs(color = "Temp") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))
est_temp

##Anomaly map
anomlimits_t <- range(points.rest3$tempanom, na.rm = T)
anombreaks_t <- seq(anomlimits_t[1], anomlimits_t[2], 5)
anomlab_t <- round(anombreaks_t)
coloranom_t <- colorRampPalette(c("lightsalmon","lightyellow","darkorchid4"))(round(maxValue(r_t) - minValue(r_t)))
anommap_t <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = "darkgray", color = "black") +
  theme_void() +
  geom_sf(data = sf_points, aes(color = points.rest3$tempanom), size = 0.5, pch = 16) +
  scale_color_gradientn(colors = coloranom_t, name = "Temperature \nAnomaly (C)", na.value = "transparent") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
anommap_t



##########
####Fossil analysis (temp)
#Fossil map
fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))
fossmap <- act_temp + geom_sf(data = fossil.sf, col = "black", size = 2, pch = 17)
fossmap

#aggregating fossil data
fossilmean <- aggregate(fossil[ , 11], list(fossil$Site, fossil$Lat, fossil$Long), mean, na.rm = TRUE)
fossilsd <- aggregate(fossil[ , 11], list(fossil$Site), sd, na.rm = TRUE)
colnames(fossilmean) <- c("Site", "Lat", "Long", "Mean")
fossilmean <- merge(fossilmean, fossilsd, by.x = "Site", by.y = "Group.1")
colnames(fossilmean)[5] <- "SD"
fossildata <- fossilmean
head(fossildata)

fossildata$fossilmbc1 <- .bincode(fossildata$Mean, breaks = mbrks3)
fossildata$fossilsdbc1 <- .bincode(fossildata$SD, breaks = sdbrks3)

fossil.sf <- st_as_sf(fossil[which(!is.na(fossil$Lat) | !is.na(fossil$Long)), ], coords = c(3, 2), crs = crs(continents1))

#calculate fossil temp estimates with 5% limits
fossilmodmax_t <- list() 
fossilmod_t <- list()
for (i in 1:length(fossildata[ , 1])) {
  if(!(is.na(fossildata$fossilmbc1[i]) | is.na(fossildata$fossilsdbc1[i]))) {
    if(length(which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])) < 1) {
      fossilmod_t[[i]] <- NA 
    } else {
      dat <- points.rest3$BIO1[which(mbc3 == fossildata$fossilmbc1[i] & sdbc3 == fossildata$fossilsdbc1[i])]
      fossilmod_t[[i]] <- density(dat[!is.na(dat)], bw = 1)
      fossilmodmax_t$temptest[i] <- fossilmod_t[[i]]$x[which.max(fossilmod_t[[i]]$y)]
      maxlength <- length(which(fossilmod_t[[i]]$x < fossilmodmax_t$temptest[i]))
      bound.05 <- (length(fossilmod_t[[i]]$x)) * 0.05
      fossilmodmax_t$minlimit[i] <- fossilmod_t[[i]]$x[maxlength - bound.05]
      fossilmodmax_t$maxlimit[i] <- fossilmod_t[[i]]$x[maxlength + bound.05]
    }
  }
}

fossildata$temptest <- fossilmodmax_t$temptest
fossildata$minlimit <- fossilmodmax_t$minlimit
fossildata$maxlimit <- fossilmodmax_t$maxlimit

#plot fossil sites in ecometric space
df1 <- data.frame(xmin = fossildata$fossilmbc1 - 1,
                  xmax = fossildata$fossilmbc1,
                  ymin = fossildata$fossilsdbc1 - 1,
                  ymax = fossildata$fossilsdbc1)
fossilEMS <- temp_ecospace + geom_rect(data = df1, aes(xmin = as.numeric(xmin), xmax = as.numeric(xmax), ymin = as.numeric(ymin), ymax = as.numeric(ymax)),colour = "black", alpha = 0, linewidth = 0.75)
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

#write.csv(fossildata, "/Users/leila/Desktop/fossil_results_temp.csv", row.names=TRUE)

fossilEMSmod1 <- fossilEMS + geom_rect(data = fossildata, aes(xmin = as.numeric(mbc3)-1,xmax = as.numeric(mbc3), ymin = as.numeric(sdbc3)-1,ymax = as.numeric(sdbc3)), colour = "green", alpha = 0, size = 1)
#export as 600 x maintain aspect ratio
fossilEMSmod1

#display vectors with marginal boxplots
temp_ecospaceNA <- ggplot(data = rdf_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = NA, name = "Temperature (C)", na.value = "transparent") +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
temp_ecospaceNA

dfvectors <- data.frame(site = fossildata$Site,
                        fxmin = fossildata$fossilmbc1 - 1,
                        fxmax = fossildata$fossilmbc1,
                        fymin = fossildata$fossilsdbc1 - 1,
                        fymax = fossildata$fossilsdbc1,
                        mxmin = fossildata$mbc - 1,
                        mxmax = fossildata$mbc,
                        mymin = fossildata$sdbc - 1,
                        mymax = fossildata$sdbc)

vector <- data.frame("site" = rep(fossildata$Site, 2),
                     "xf" = rowMeans(dfvectors[, 2:3], na.rm = T),
                     "yf" = rowMeans(dfvectors[, 4:5], na.rm = T),
                     "xm" = rowMeans(dfvectors[, 6:7], na.rm = T),
                     "ym" = rowMeans(dfvectors[, 8:9], na.rm = T))

vector2 <- data.frame("site" = rep(fossildata$Site, 2), "group" = c(rep("fossil", length(fossildata$Site)), rep("modern", length(fossildata$Site))), "x" = c(rowMeans(dfvectors[, 2:3], na.rm = T), rowMeans(dfvectors[, 6:7], na.rm = T)), "y" = c(rowMeans(dfvectors[, 4:5], na.rm = T), rowMeans(dfvectors[, 8:9], na.rm = T)))

vectorspace <- temp_ecospaceNA + 
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

p <- ggMarginal(p = vectorspace, type = "histogram", 
                margins = "both", 
                size = 6,
                groupColour = TRUE,
                groupFill = TRUE)
p


####Sensitivity Code (precip)
#number of bins
bins <- 25

#function for returning bin codes
bincodes <- function(traits) {
  trait_range <- range(traits, na.rm = T) #take the range of the means
  trait_range[1] <- trait_range[1]-0.000001 #add a little amount on each side of the range
  trait_range[2] <- trait_range[2]+0.000001
  brks <- seq(trait_range[1], trait_range[2], diff(trait_range)/bins) #split the range into the number of specified equal parts
  bc <- .bincode(traits, breaks = brks) #assign a bin code to each of the bins
  return(bc)
}

#function to retrieve bin breaks
binbreaks <- function(traits) {
  trait_range <- range(traits, na.rm = T) #take the range of the means
  trait_range[1] <- trait_range[1]-0.000001 #add a little amount on each side of the range
  trait_range[2] <- trait_range[2]+0.000001
  brks <- seq(trait_range[1], trait_range[2], diff(trait_range)/bins) #split the range into equal parts
  return(brks)
}

#define sample sizes
samplesizes <- seq(100, 10000, 1000) 

#define number of bootstrap iterations for each
sens_it <- 20

#define split of testing and training data
test_split <- 0.2

#set vars to collect data
#sensitivity is within the model and testing is transferability
anomp_sens_temp <- array(NA, dim = c(length(samplesizes),sens_it))
anomp_test_temp <- array(NA, dim = c(length(samplesizes),sens_it))
cor_sens_temp <- array(NA, dim = c(length(samplesizes),sens_it))
cor_test_temp <- array(NA, dim = c(length(samplesizes),sens_it))

#loop to downsample and split dataset into training and testing parts (80/20)
for(samp_sizes in 1:length(samplesizes)){
  testsample <- sample(points.rest3$GlobalID, size = samplesizes[samp_sizes], replace = F)
  downdata <- points.rest3[points.rest3$GlobalID %in% testsample,]
  
  for (sens in 1:sens_it){
    training <- sample(x = 1:length(downdata[,1]), size = length(downdata[,1])*(1 - test_split))
    testing <- (1:length(downdata[,1]))[!(1:length(downdata[,1]) %in% training)]
    training_data <- downdata[training,]
    testing_data <- downdata[testing,]
    
    ## Ecometric Spaces ##
    #Group the trait distributions into bins x bins trait bins & calculate the range
    training_data$mbc.tr <- bincodes(training_data$mean_RBL) #ecometric space calculations, means
    training_data$sdbc.tr <- bincodes(training_data$sd_RBL)  #ecometric space calculations, standard deviations
    
    ################################
    #Calculate difference between observed and trait-predicted environment
    maxlike <- array(NA, dim = length(training_data[ , 1])) 
    obs_env <- array(NA, dim = length(training_data[ , 1])) 
    anomp <- array(NA, dim = length(training_data[ , 1]))
    
    mod <- list()
    
    for (i in 1:length(training_data[ , 1])) {
      if(!(is.na(training_data$mbc.tr[i]) | is.na(training_data$sdbc.tr[i]) | is.na(training_data$BIO1[i]))) {
        if(length(which(training_data$mbc.tr == training_data$mbc.tr[i] & training_data$sdbc.tr == training_data$sdbc.tr[i])) < 1) {
          mod[[i]] <- NA 
        } else {
          dat <- training_data$BIO1[which((training_data$mbc.tr == training_data$mbc.tr[i] & training_data$sdbc.tr == training_data$sdbc.tr[i]))]
          mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
          maxlike[i] <- mod[[i]]$x[which.max(mod[[i]]$y)] 
          obs_env[i] <- training_data$BIO1[i]
          anomp[i] <- obs_env[i] - maxlike[i]
        }
      }
    }
    
    #hold anomaly for each testing training split
    anomp <- anomp[!is.infinite(anomp)]
    anomp <- anomp[!is.na(anomp)]
    anomp_sens[samp_sizes,sens]<- mean(abs(anomp), na.rm = T)
    cor_sens[samp_sizes,sens] <- cor(maxlike, obs_env, use = "pairwise.complete.obs")
    
    #testing data, get bin breaks from training data
    mbrks <- binbreaks(training_data$mean_RBL)
    sdbrks <- binbreaks(training_data$sd_RBL)
    
    #get bin codes for testing data
    testing_data$mbc.te <- .bincode(testing_data$mean_RBL, breaks = mbrks)
    testing_data$sdbc.te <-.bincode(testing_data$sd_RBL, breaks = sdbrks)
    
    #reset the values
    maxlike <- array(NA, dim = length(training_data[ , 1])) 
    obs_env <- array(NA, dim = length(training_data[ , 1])) 
    anomp <- array(NA, dim = length(training_data[ , 1]))
    
    mod <- list()
    
    for (j in 1:length(testing_data[ , 1])) {
      if(!(is.na(testing_data$mbc.te[j]) | is.na(testing_data$sdbc.te[j]) | is.na(testing_data$BIO1[j]))) {
        if(length(which(training_data$mbc.tr == testing_data$mbc.te[j] & training_data$sdbc.tr == testing_data$sdbc.te[j])) < 1) {
          mod[[j]] <- NA 
        } else {
          dat <- training_data$BIO1[which(training_data$mbc.tr == testing_data$mbc.te[j] & training_data$sdbc.tr == testing_data$sdbc.te[j])]
          mod[[j]] <- density(dat[!is.na(dat)], bw = 1) 
          maxlike[j] <- mod[[j]]$x[which.max(mod[[j]]$y)] 
          obs_env[j] <- testing_data$BIO1[j]
          anomp[j] <- obs_env[j] - maxlike[j]
        }
      }
    }
    anomp <- anomp[!is.infinite(anomp)]
    anomp <- anomp[!is.na(anomp)]
    anomp_test[samp_sizes,sens] <- mean(abs(anomp), na.rm = T)
    cor_test[samp_sizes,sens] <- cor(maxlike, obs_env, use = "pairwise.complete.obs")
  }
}

###################################################
training_cor <- data.frame(samplesizes = samplesizes, cor_sens)
training_cor <- melt(training_cor, id.vars = 1)

testing_cor <- data.frame(samplesizes = samplesizes, cor_test)
testing_cor <- melt(testing_cor, id.vars = 1)

testing_plot <- ggplot(testing_cor, aes(samplesizes, testing_cor[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  ylim(0.4, 1) +
  xlab("Communities") +
  ylab("Testing - Cor") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic() 

training_plot <- ggplot(training_cor, aes(samplesizes, training_cor[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  ylim(0.4, 1) +
  xlab("Communities") +
  ylab("Training - Cor") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic()

ggarrange(training_plot, testing_plot,
          labels = c("A", "B"), ncol = 1, nrow = 2)

###################################################

training_anom <- data.frame(samplesizes = samplesizes, anomp_sens)
training_anom <- melt(training_anom, id.vars = 1)

testing_anom <- data.frame(samplesizes = samplesizes, anomp_test)
testing_anom <- melt(testing_anom, id.vars = 1)
#  ylim(0, 1) +

testing_plot <- ggplot(testing_anom, aes(samplesizes, testing_anom[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  xlab("Communities") +
  ylab("Testing - Anom") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic() 

training_plot <- ggplot(training_anom, aes(samplesizes, training_anom[,3])) +   # Draw ggplot2 scatterplot with smooth curve
  geom_point(col = "gray") +
  xlab("Communities") +
  ylab("Training - Anom") +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  theme_classic()

ggarrange(training_plot, testing_plot, labels = c("A", "B"), ncol = 1, nrow = 2)




#########Continental trends
###################
#North America ONLY
#Polygon shape files for continental shape
continents<-shapefile("continents.shp")
head(continents)
n_america<-continents[continents$CONTINENT=="North America",]
n_america.sf<-st_as_sf(n_america)

##convert sampling points to spatial points and sample species at each sampling locality
sp_nAm<-SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(n_america))
sp_nAm<- raster ::intersect(sp, n_america)
head(sp_nAm)
length(sp_nAm)

#geography <- shapefile("TERRESTRIAL_MAMMALS.shp")
#save(geography, file="geography.Rdata")
#load("geography.Rdata")
#head(geography)

#Sample species at each sampling location.
#o_nAm <- over(sp_nAm, geography, returnList = T)
#save(o_nAm, file="o_nAm.Rdata")
load("o_nAm.Rdata")
head(o_nAm)
length(o_nAm)

points.nAm<-points[points$cont=="North America",]
dim(points.nAm)

#Transform precip 
points.nAm$LogPrecip <- log(points.nAm$BIO12+1)
#points.nAm$TempCube <- ((points.nAm$BIO1/10) ^ 3)

#Summarize traits for community-level distribution by calculating mean and standard deviation
points.nAm$mean_RBL <- unlist(lapply(o_nAm, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
points.nAm$sd_RBL <- unlist(lapply(o_nAm, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))

##calculate richness at each point
points.nAm$RBL_richness <- unlist(lapply(o_nAm, function(x) length(traits[x$binomial, "RBL"])))
hist(points.nAm$RBL_richness)

#Calculate richness at each point with traits data
points.nAm$count <- unlist(lapply(o_nAm, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))

####North America
#Min 3 species per community
points.nAm3<-points.nAm[points.nAm$count>2,]
head(points.nAm3)
dim(points.nAm3)

summary(points.nAm3$mean_RBL)
summary(points.nAm3$sd_RBL)

#Bin distributions 25x25
mRBLnAm3<-range(points.nAm3$mean_RBL, na.rm=T)
sdRBLnAm3 <- range(points.nAm3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrksnAm3 <- seq(mRBLnAm3[1]-0.001, mRBLnAm3[2]+0.001, diff(mRBLnAm3)/25)
sdbrksnAm3 <- seq(sdRBLnAm3[1]-0.001, sdRBLnAm3[2]+0.001, diff(sdRBLnAm3)/25)

#Assign bin codes
mbcnAm3<-.bincode(points.nAm3$mean_RBL, breaks = mbrksnAm3)
sdbcnAm3 <-.bincode(points.nAm3$sd_RBL, breaks = sdbrksnAm3)

#Assign mean & sd bin code to each community
points.nAm3$mbcnAm3<-mbcnAm3
points.nAm3$sdbcnAm3<-sdbcnAm3

###################
#North America: Precipitation
mod <- list()
obj_nAm3 <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcnAm3 == i & sdbcnAm3 == j, na.rm = T) != 0) {
      dat <- points.nAm3$LogPrecip[which(mbcnAm3 == i & sdbcnAm3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_nAm3[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_nAm3 <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_nAm3 <- setValues(r_nAm3, obj_nAm3)
rdf_nAm3 <- as.data.frame(r_nAm3, xy=TRUE)
names(rdf_nAm3) <- c("Mean", "SD", "Precip")

color_nAm3 <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r_nAm3) - minValue(r_nAm3)))
precip_ecospace_nAm <-ggplot(data = rdf_nAm3) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color_nAm3), name = "Precipitation (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_nAm

######Calculate maximum likelihood for all bins
modmax_nAm3 <- array(NA, dim = length(points.nAm3[ , 1]))
mod <- list()

for (i in 1:length(points.nAm3[ , 1])) {
  if(!(is.na(points.nAm3$mbcnAm3[i]) | is.na(points.nAm3$sdbcnAm3[i]))) {
    dat <- points.nAm3$LogPrecip[which(points.nAm3$mbcnAm3 == points.nAm3$mbcnAm3[i] & points.nAm3$sdbcnAm3 == points.nAm3$sdbcnAm3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_nAm3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_nAm3

hist(modmax_nAm3)

#calculate anomaly
anom_precip_nAm3 <- points.nAm3$LogPrecip - modmax_nAm3
anom_precip_NL_nAm <- exp(anom_precip_nAm3) - 1
plot(anom_precip_nAm3, anom_precip_NL_nAm)

points.nAm3$precipanom <- as.numeric(anom_precip_nAm3)
summary(points.nAm3$precipanom, na.rm = T)
points.nAm3$precipest <- as.numeric(modmax_nAm3)
summary(anom_precip_NL_nAm, na.rm = T)

#Correlation of observed ~ expected 
precipcor_nAm <- lm(points.nAm3$precipest ~ points.nAm3$LogPrecip)
precipcor_nAm
summary(precipcor_nAm)
cor_precipest_nAm <- cor.test(points.nAm3$precipest, points.nAm3$LogPrecip, method = "pearson")
#cor=0.6729513 
cor_precipest_nAm


######
#North America: Temperature
mod <- list()
obj_nAm3_t <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcnAm3 == i & sdbcnAm3 == j, na.rm = T) != 0) {
      dat <- points.nAm3$BIO1[which(mbcnAm3 == i & sdbcnAm3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_nAm3_t[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_nAm3_t <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_nAm3_t <- setValues(r_nAm3_t, obj_nAm3_t)
rdf_nAm3_t <- as.data.frame(r_nAm3_t, xy=TRUE)
names(rdf_nAm3_t) <- c("Mean", "SD", "Temp")

color_nAm3_t <- colorRampPalette(c("darkorchid4","darkorchid","violetred", "deeppink","lightcoral", "lightsalmon"))(round(maxValue(r_nAm3_t) - minValue(r_nAm3_t)))
precip_ecospace_nAm_t <-ggplot(data = rdf_nAm3_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = rev(color_nAm3_t), name = "Temperature (C)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_nAm_t

######Calculate maximum likelihood for all bins
modmax_nAm3_t <- array(NA, dim = length(points.nAm3[ , 1]))
mod <- list()

for (i in 1:length(points.nAm3[ , 1])) {
  if(!(is.na(points.nAm3$mbcnAm3[i]) | is.na(points.nAm3$sdbcnAm3[i]))) {
    dat <- points.nAm3$BIO1[which(points.nAm3$mbcnAm3 == points.nAm3$mbcnAm3[i] & points.nAm3$sdbcnAm3 == points.nAm3$sdbcnAm3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_nAm3_t[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_nAm3_t

hist(modmax_nAm3_t)


#calculate temperature anomaly
anom_temp_nAm <- points.nAm3$BIO1 - modmax_nAm3_t

# #Inverse of cubed temp (cube root)
# #R has a bit of trouble with cube root w/ negative numbers, so a new function is needed 
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }

# #Non-cube temp anomaly
# anom_temp_NC_nAm <- Math.cbrt(anom_temp_nAm)
# summary(anom_temp_NC_nAm)
# plot(anom_temp_nAm, anom_temp_NC_nAm)

points.nAm3$tempanom <- as.numeric(anom_temp_nAm)
summary(points.nAm3$tempanom, na.rm = T)
points.nAm3$tempest <- as.numeric(modmax_nAm3_t)
summary(anom_temp_NC_nAm, na.rm = T)

#Correlation of observed ~ expected 
tempcor_nAm <- lm(points.nAm3$tempest ~ points.nAm3$BIO1)
tempcor_nAm
summary(tempcor_nAm)
cor_tempest_nAm <- cor.test(points.nAm3$tempest, points.nAm3$BIO1, method = "pearson")
#cor=0.7629494 
cor_tempest_nAm



###################
#South America ONLY
#Polygon shape files for continental shape
continents<-shapefile("continents.shp")
head(continents)
s_america<-continents[continents$CONTINENT=="South America",]
s_america.sf<-st_as_sf(s_america)

##convert sampling points to spatial points and sample species at each sampling locality
sp_sAm<-SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(s_america))
sp_sAm<- raster ::intersect(sp, s_america)
head(sp_sAm)
length(sp_sAm)

#Sample species at each sampling location.
#o_sAm <- over(sp_sAm, geography, returnList = T)
#save(o_sAm, file="o_sAm.Rdata")
load("o_sAm.Rdata")
head(o_sAm)
length(o_sAm)

points.sAm<-points[points$cont=="South America",]
dim(points.sAm)

#Transform log 
points.sAm$LogPrecip <- log(points.sAm$BIO12+1)
#points.sAm$TempCube <- ((points.sAm$BIO1/10) ^ 3)

#Summarize traits for community-level distribution by calculating mean and standard deviation
points.sAm$mean_RBL <- unlist(lapply(o_sAm, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
points.sAm$sd_RBL <- unlist(lapply(o_sAm, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))

##calculate richness at each point
points.sAm$RBL_richness <- unlist(lapply(o_sAm, function(x) length(traits[x$binomial, "RBL"])))
hist(points.sAm$RBL_richness)

#Calculate richness at each point with traits data
points.sAm$count <- unlist(lapply(o_sAm, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))

####South America
#Min 3 species per community
points.sAm3<-points.sAm[points.sAm$count>2,]
head(points.sAm3)
dim(points.sAm3)

summary(points.sAm3$mean_RBL)
summary(points.sAm3$sd_RBL)

#Bin distributions 25x25
mRBLsAm3<-range(points.sAm3$mean_RBL, na.rm=T)
sdRBLsAm3 <- range(points.sAm3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrkssAm3 <- seq(mRBLsAm3[1]-0.001, mRBLsAm3[2]+0.001, diff(mRBLsAm3)/25)
sdbrkssAm3 <- seq(sdRBLsAm3[1]-0.001, sdRBLsAm3[2]+0.001, diff(sdRBLsAm3)/25)

#Assign bin codes
mbcsAm3<-.bincode(points.sAm3$mean_RBL, breaks = mbrkssAm3)
sdbcsAm3 <-.bincode(points.sAm3$sd_RBL, breaks = sdbrkssAm3)

#Assign mean & sd bin code to each community
points.sAm3$mbcsAm3<-mbcsAm3
points.sAm3$sdbcsAm3<-sdbcsAm3

###################
#South America: Precipitation
mod <- list()
obj_sAm3 <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcsAm3 == i & sdbcsAm3 == j, na.rm = T) != 0) {
      dat <- points.sAm3$LogPrecip[which(mbcsAm3 == i & sdbcsAm3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_sAm3[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_sAm3 <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_sAm3 <- setValues(r_sAm3, obj_sAm3)
rdf_sAm3 <- as.data.frame(r_sAm3, xy=TRUE)
names(rdf_sAm3) <- c("Mean", "SD", "Precip")

color_sAm3 <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r_sAm3) - minValue(r_sAm3)))
precip_ecospace_sAm <-ggplot(data = rdf_sAm3) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color_sAm3), name = "Precipitation (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_sAm

######Calculate maximum likelihood for all bins
modmax_sAm3 <- array(NA, dim = length(points.sAm3[ , 1]))
mod <- list()

for (i in 1:length(points.sAm3[ , 1])) {
  if(!(is.na(points.sAm3$mbcsAm3[i]) | is.na(points.sAm3$sdbcsAm3[i]))) {
    dat <- points.sAm3$LogPrecip[which(points.sAm3$mbcsAm3 == points.sAm3$mbcsAm3[i] & points.sAm3$sdbcsAm3 == points.sAm3$sdbcsAm3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_sAm3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_sAm3

hist(modmax_sAm3)

#calculate anomaly
anom_precip_sAm3 <- points.sAm3$LogPrecip - modmax_sAm3
anom_precip_NL_sAm <- exp(anom_precip_sAm3) - 1
plot(anom_precip_sAm3, anom_precip_NL_sAm)

points.sAm3$precipanom <- as.numeric(anom_precip_sAm3)
summary(points.sAm3$precipanom, na.rm = T)
points.sAm3$precipest <- as.numeric(modmax_sAm3)
summary(anom_precip_NL_sAm, na.rm = T)

#Correlation of observed ~ expected 
precipcor_sAm <- lm(points.sAm3$precipest ~ points.sAm3$LogPrecip)
precipcor_sAm
summary(precipcor_sAm)
cor_precipest_sAm <- cor.test(points.sAm3$precipest, points.sAm3$LogPrecip, method = "pearson")
#cor=0.7993706 
cor_precipest_sAm


######
#South America: Temperature
mod <- list()
obj_sAm3_t <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcsAm3 == i & sdbcsAm3 == j, na.rm = T) != 0) {
      dat <- points.sAm3$BIO1[which(mbcsAm3 == i & sdbcsAm3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_sAm3_t[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_sAm3_t <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_sAm3_t <- setValues(r_sAm3_t, obj_sAm3_t)
rdf_sAm3_t <- as.data.frame(r_sAm3_t, xy=TRUE)
names(rdf_sAm3_t) <- c("Mean", "SD", "Temp")

color_sAm3_t <- colorRampPalette(c("darkorchid4","darkorchid","violetred", "deeppink","lightcoral", "lightsalmon"))(round(maxValue(r_sAm3_t) - minValue(r_sAm3_t)))
precip_ecospace_sAm_t <-ggplot(data = rdf_sAm3_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = rev(color_sAm3_t), name = "Temperature (C)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_sAm_t

######Calculate maximum likelihood for all bins
modmax_sAm3_t <- array(NA, dim = length(points.sAm3[ , 1]))
mod <- list()

for (i in 1:length(points.sAm3[ , 1])) {
  if(!(is.na(points.sAm3$mbcsAm3[i]) | is.na(points.sAm3$sdbcsAm3[i]))) {
    dat <- points.sAm3$BIO1[which(points.sAm3$mbcsAm3 == points.sAm3$mbcsAm3[i] & points.sAm3$sdbcsAm3 == points.sAm3$sdbcsAm3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_sAm3_t[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_sAm3_t
hist(modmax_sAm3_t)

#calculate temperature anomaly
anom_temp_sAm <- points.sAm3$BIO1 - modmax_sAm3_t

# #Inverse of cubed temp (cube root)
# #R has a bit of trouble with cube root w/ negative numbers, so a new function is needed 
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }

# #Non-cube temp anomaly
# anom_temp_NC_sAm <- Math.cbrt(anom_temp_sAm)
# summary(anom_temp_NC_sAm)
# plot(anom_temp_sAm, anom_temp_NC_sAm)

points.sAm3$tempanom <- as.numeric(anom_temp_sAm)
summary(points.sAm3$tempanom, na.rm = T)
points.sAm3$tempest <- as.numeric(modmax_sAm3_t)
summary(anom_temp_NC_sAm, na.rm = T)

#Correlation of observed ~ expected 
tempcor_sAm <- lm(points.sAm3$tempest ~ points.sAm3$BIO1)
tempcor_sAm
summary(tempcor_sAm)
cor_tempest_sAm <- cor.test(points.sAm3$tempest, points.sAm3$BIO1, method = "pearson")
#cor=0.7546363 
cor_tempest_sAm



###################
#Europe ONLY
#Polygon shape files for continental shape
europe<-continents[continents$CONTINENT=="Europe",]
europe.sf<-st_as_sf(europe)

##convert sampling points to spatial points and sample species at each sampling locality
sp_eu<-SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(europe))
sp_eu<- raster ::intersect(sp_eu, europe)
head(sp_eu)
length(sp_eu)

#Sample species at each sampling location.
#o_eu <- over(sp_eu, geography, returnList = T)
#save(o_eu, file="o_eu.Rdata")
load("o_eu.Rdata")
head(o_eu)
length(o_eu)

points.eu<-points[points$cont=="Europe",]
dim(points.eu)

#Transform log
points.eu$LogPrecip <- log(points.eu$BIO12+1)
#points.eu$TempCube <- ((points.eu$BIO1/10) ^ 3)

#Summarize traits for community-level distribution by calculating mean and standard deviation
points.eu$mean_RBL <- unlist(lapply(o_eu, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
points.eu$sd_RBL <- unlist(lapply(o_eu, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))

##calculate richness at each point
points.eu$RBL_richness <- unlist(lapply(o_eu, function(x) length(traits[x$binomial, "RBL"])))
hist(points.eu$RBL_richness)

#Calculate richness at each point with traits data
points.eu$count <- unlist(lapply(o_eu, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))

####Europe
#Min 3 species per community
points.eu3<-points.eu[points.eu$count>2,]
head(points.eu3)
dim(points.eu3)

summary(points.eu3$mean_RBL)
summary(points.eu3$sd_RBL)

#Bin distributions 25x25
mRBLeu3<-range(points.eu3$mean_RBL, na.rm=T)
sdRBLeu3 <- range(points.eu3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrkseu3 <- seq(mRBLeu3[1]-0.001, mRBLeu3[2]+0.001, diff(mRBLeu3)/25)
sdbrkseu3 <- seq(sdRBLeu3[1]-0.001, sdRBLeu3[2]+0.001, diff(sdRBLeu3)/25)

#Assign bin codes
mbceu3<-.bincode(points.eu3$mean_RBL, breaks = mbrkseu3)
sdbceu3 <-.bincode(points.eu3$sd_RBL, breaks = sdbrkseu3)

#Assign mean & sd bin code to each community
points.eu3$mbceu3<-mbceu3
points.eu3$sdbceu3<-sdbceu3


#Europe: Precipitation
mod <- list()
obj_eu3 <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbceu3 == i & sdbceu3 == j, na.rm = T) != 0) {
      dat <- points.eu3$LogPrecip[which(mbceu3 == i & sdbceu3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_eu3[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_eu3 <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_eu3 <- setValues(r_eu3, obj_eu3)
rdf_eu3 <- as.data.frame(r_eu3, xy=TRUE)
names(rdf_eu3) <- c("Mean", "SD", "Precip")

color_eu3 <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r_eu3) - minValue(r_eu3)))
precip_ecospace_eu <-ggplot(data = rdf_eu3) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color_eu3), name = "Precipitation (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_eu

######Calculate maximum likelihood for all bins
modmax_eu3 <- array(NA, dim = length(points.eu3[ , 1]))
mod <- list()

for (i in 1:length(points.eu3[ , 1])) {
  if(!(is.na(points.eu3$mbceu3[i]) | is.na(points.eu3$sdbceu3[i]))) {
    dat <- points.eu3$LogPrecip[which(points.eu3$mbceu3 == points.eu3$mbceu3[i] & points.eu3$sdbceu3 == points.eu3$sdbceu3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_eu3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_eu3
hist(modmax_eu3)

#calculate anomaly
anom_precip_eu3 <- points.eu3$LogPrecip - modmax_eu3
anom_precip_NL_eu <- exp(anom_precip_eu3) - 1
plot(anom_precip_eu3, anom_precip_NL_eu)

points.eu3$precipanom <- as.numeric(anom_precip_eu3)
summary(points.eu3$precipanom, na.rm = T)
points.eu3$precipest <- as.numeric(modmax_eu3)
summary(anom_precip_NL_eu, na.rm = T)

#Correlation of observed ~ expected 
precipcor_eu <- lm(points.eu3$precipest ~ points.eu3$LogPrecip)
precipcor_eu
summary(precipcor_eu)
cor_precipest_eu <- cor.test(points.eu3$precipest, points.eu3$LogPrecip, method = "pearson")
#cor=0.682787 
cor_precipest_eu


######
#Europe: Temperature
mod <- list()
obj_eu3_t <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbceu3 == i & sdbceu3 == j, na.rm = T) != 0) {
      dat <- points.eu3$BIO1[which(mbceu3 == i & sdbceu3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_eu3_t[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_eu3_t <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_eu3_t <- setValues(r_eu3_t, obj_eu3_t)
rdf_eu3_t <- as.data.frame(r_eu3_t, xy=TRUE)
names(rdf_eu3_t) <- c("Mean", "SD", "Temp")

color_eu3_t <- colorRampPalette(c("darkorchid4","darkorchid","violetred", "deeppink","lightcoral", "lightsalmon"))(round(maxValue(r_eu3_t) - minValue(r_eu3_t)))
precip_ecospace_eu_t <-ggplot(data = rdf_eu3_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = rev(color_eu3_t), name = "Temperature (C)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_eu_t

######Calculate maximum likelihood for all bins
modmax_eu3_t <- array(NA, dim = length(points.eu3[ , 1]))
mod <- list()

for (i in 1:length(points.eu3[ , 1])) {
  if(!(is.na(points.eu3$mbceu3[i]) | is.na(points.eu3$sdbceu3[i]))) {
    dat <- points.eu3$BIO1[which(points.eu3$mbceu3 == points.eu3$mbceu3[i] & points.eu3$sdbceu3 == points.eu3$sdbceu3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_eu3_t[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_eu3_t
hist(modmax_eu3_t)

#calculate temperature anomaly
anom_temp_eu <- points.eu3$BIO1 - modmax_eu3_t

# #Inverse of cubed temp (cube root)
# #R has a bit of trouble with cube root w/ negative numbers, so a new function is needed 
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }

# #Non-cube temp anomaly
# anom_temp_NC_eu <- Math.cbrt(anom_temp_eu)
# summary(anom_temp_NC_eu)
# plot(anom_temp_eu, anom_temp_NC_eu)

points.eu3$tempanom <- as.numeric(anom_temp_eu)
summary(points.eu3$tempanom, na.rm = T)
points.eu3$tempest <- as.numeric(modmax_eu3_t)
summary(anom_temp_NC_eu, na.rm = T)

#Correlation of observed ~ expected 
tempcor_eu <- lm(points.eu3$tempest ~ points.eu3$BIO1)
tempcor_eu
summary(tempcor_eu)
cor_tempest_eu <- cor.test(points.eu3$tempest, points.eu3$BIO1, method = "pearson")
#cor=0.5512926 
cor_tempest_eu



###################
#Asia ONLY
#Polygon shape files for continental shape
asia<-continents[continents$CONTINENT=="Asia",]
asia.sf<-st_as_sf(asia)

##convert sampling points to spatial points and sample species at each sampling locality
sp_as<-SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(asia))
sp_as<- raster ::intersect(sp_as, asia)
head(sp_as)
length(sp_as)

#Sample species at each sampling location.
#o_as <- over(sp_as, geography, returnList = T)
#save(o_as, file="o_as.Rdata")
load("o_as.Rdata")
head(o_as)
length(o_as)

points.as<-points[points$cont=="Asia",]
dim(points.as)

#Transform precip
points.as$LogPrecip <- log(points.as$BIO12+1)
#points.as$TempCube <- ((points.as$BIO1/10) ^ 3)

#Summarize traits for community-level distribution by calculating mean and standard deviation
points.as$mean_RBL <- unlist(lapply(o_as, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
points.as$sd_RBL <- unlist(lapply(o_as, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))

##calculate richness at each point
points.as$RBL_richness <- unlist(lapply(o_as, function(x) length(traits[x$binomial, "RBL"])))
hist(points.as$RBL_richness)

#Calculate richness at each point with traits data
points.as$count <- unlist(lapply(o_as, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))

points.as3<-points.as[points.as$count>2,]
head(points.as3)
dim(points.as3)

summary(points.as3$mean_RBL)
summary(points.as3$sd_RBL)

#Bin distributions 25x25
mRBLas3<-range(points.as3$mean_RBL, na.rm=T)
sdRBLas3 <- range(points.as3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrksas3 <- seq(mRBLas3[1]-0.001, mRBLas3[2]+0.001, diff(mRBLas3)/25)
sdbrksas3 <- seq(sdRBLas3[1]-0.001, sdRBLas3[2]+0.001, diff(sdRBLas3)/25)

#Assign bin codes
mbcas3<-.bincode(points.as3$mean_RBL, breaks = mbrksas3)
sdbcas3 <-.bincode(points.as3$sd_RBL, breaks = sdbrksas3)

#Assign mean & sd bin code to each community
points.as3$mbcas3<-mbcas3
points.as3$sdbcas3<-sdbcas3


#Asia: Precipitation
mod <- list()
obj_as3 <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcas3 == i & sdbcas3 == j, na.rm = T) != 0) {
      dat <- points.as3$LogPrecip[which(mbcas3 == i & sdbcas3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_as3[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_as3 <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_as3 <- setValues(r_as3, obj_as3)
rdf_as3 <- as.data.frame(r_as3, xy=TRUE)
names(rdf_as3) <- c("Mean", "SD", "Precip")

color_as3 <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r_as3) - minValue(r_as3)))
precip_ecospace_as <-ggplot(data = rdf_as3) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color_as3), name = "Precipitation (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_as

######Calculate maximum likelihood for all bins
modmax_as3 <- array(NA, dim = length(points.as3[ , 1]))
mod <- list()

for (i in 1:length(points.as3[ , 1])) {
  if(!(is.na(points.as3$mbcas3[i]) | is.na(points.as3$sdbcas3[i]))) {
    dat <- points.as3$LogPrecip[which(points.as3$mbcas3 == points.as3$mbcas3[i] & points.as3$sdbcas3 == points.as3$sdbcas3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_as3[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_as3
hist(modmax_as3)

#calculate anomaly
anom_precip_as3 <- points.as3$LogPrecip - modmax_as3
anom_precip_NL_as <- exp(anom_precip_as3) - 1
plot(anom_precip_as3, anom_precip_NL_as)

points.as3$precipanom <- as.numeric(anom_precip_as3)
summary(points.as3$precipanom, na.rm = T)
points.as3$precipest <- as.numeric(modmax_as3)
summary(anom_precip_NL_as, na.rm = T)

#Correlation of observed ~ expected 
precipcor_as <- lm(points.as3$precipest ~ points.as3$LogPrecip)
precipcor_as
summary(precipcor_as)
cor_precipest_as <- cor.test(points.as3$precipest, points.as3$LogPrecip, method = "pearson")
#cor=0.8069793 
cor_precipest_as


######
#Asia: Temperature
mod <- list()
obj_as3_t <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcas3 == i & sdbcas3 == j, na.rm = T) != 0) {
      dat <- points.as3$BIO1[which(mbcas3 == i & sdbcas3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_as3_t[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_as3_t <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_as3_t <- setValues(r_as3_t, obj_as3_t)
rdf_as3_t <- as.data.frame(r_as3_t, xy=TRUE)
names(rdf_as3_t) <- c("Mean", "SD", "Temp")

color_as3_t <- colorRampPalette(c("darkorchid4","darkorchid","violetred", "deeppink","lightcoral", "lightsalmon"))(round(maxValue(r_as3_t) - minValue(r_as3_t)))
precip_ecospace_as_t <-ggplot(data = rdf_as3_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = rev(color_as3_t), name = "Temperature (C)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_as_t

######Calculate maximum likelihood for all bins
modmax_as3_t <- array(NA, dim = length(points.as3[ , 1]))
mod <- list()

for (i in 1:length(points.as3[ , 1])) {
  if(!(is.na(points.as3$mbcas3[i]) | is.na(points.as3$sdbcas3[i]))) {
    dat <- points.as3$BIO1[which(points.as3$mbcas3 == points.as3$mbcas3[i] & points.as3$sdbcas3 == points.as3$sdbcas3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_as3_t[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_as3_t
hist(modmax_as3_t)

#calculate temperature anomaly
anom_temp_as <- points.as3$BIO1 - modmax_as3_t
# #Inverse of cubed temp (cube root)
# #R has a bit of trouble with cube root w/ negative numbers, so a new function is needed 
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }
# #Non-cube temp anomaly
# anom_temp_NC_as <- Math.cbrt(anom_temp_as)
# summary(anom_temp_NC_as)
# plot(anom_temp_as, anom_temp_NC_as)

points.as3$tempanom <- as.numeric(anom_temp_as)
summary(points.as3$tempanom, na.rm = T)
points.as3$tempest <- as.numeric(modmax_as3_t)
summary(anom_temp_NC_as, na.rm = T)

#Correlation of observed ~ expected 
tempcor_as <- lm(points.as3$tempest ~ points.as3$BIO1)
tempcor_as
summary(tempcor_as)
cor_tempest_as <- cor.test(points.as3$tempest, points.as3$BIO1, method = "pearson")
#cor=0.4771841 
cor_tempest_as


###################
#Africa ONLY
#Polygon shape files for continental shape
africa<-continents[continents$CONTINENT=="Africa",]
africa.sf<-st_as_sf(africa)

##convert sampling points to spatial points and sample species at each sampling locality
sp_af<-SpatialPointsDataFrame(points[ , 2:3], data = points, proj4string = crs(africa))
sp_af<- raster ::intersect(sp_af, africa)
head(sp_af)
length(sp_af)

#Sample species at each sampling location.
#o_af <- over(sp_af, geography, returnList = T)
#save(o_af, file="o_af.Rdata")
load("o_af.Rdata")
head(o_af)
length(o_af)

points.af<-points[points$cont=="Africa",]
dim(points.af)

#Transform precip
points.af$LogPrecip <- log(points.af$BIO12+1)
#points.af$TempCube <- ((points.af$BIO1/10) ^ 3)

#Summarize traits for community-level distribution by calculating mean and standard deviation
points.af$mean_RBL <- unlist(lapply(o_af, function(x) mean(traits[x$binomial, "RBL"], na.rm = T)))
points.af$sd_RBL <- unlist(lapply(o_af, function(x) sd(traits[x$binomial, "RBL"], na.rm = T)))

##calculate richness at each point
points.af$RBL_richness <- unlist(lapply(o_af, function(x) length(traits[x$binomial, "RBL"])))
hist(points.af$RBL_richness)

#Calculate richness at each point with traits data
points.af$count <- unlist(lapply(o_af, function(x) sum(!is.na(traits[x$binomial, "RBL"]))))

points.af3<-points.af[points.af$count>2,]
head(points.af3)
dim(points.af3)

summary(points.af$mean_RBL)
summary(points.af$sd_RBL)

#Bin distributions 25x25
mRBLaf3<-range(points.af3$mean_RBL, na.rm=T)
sdRBLaf3 <- range(points.af3$sd_RBL, na.rm = T)

#break points for mean & sd
mbrksaf3 <- seq(mRBLaf3[1]-0.001, mRBLaf3[2]+0.001, diff(mRBLaf3)/25)
sdbrksaf3 <- seq(sdRBLaf3[1]-0.001, sdRBLaf3[2]+0.001, diff(sdRBLaf3)/25)

#Assign bin codes
mbcaf3<-.bincode(points.af3$mean_RBL, breaks = mbrksaf3)
sdbcaf3 <-.bincode(points.af3$sd_RBL, breaks = sdbrksaf3)

#Assign mean & sd bin code to each community
points.af3$mbcaf3<-mbcaf3
points.af3$sdbcaf3<-sdbcaf3


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

color_af3 <- colorRampPalette(c("steelblue","cornflowerblue","mediumaquamarine", "cyan3","aquamarine2", "aquamarine", "palegreen"))(round(maxValue(r_af3) - minValue(r_af3)))
precip_ecospace_af <-ggplot(data = rdf_af3) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Precip)) +
  scale_fill_gradientn(colors = rev(color_af3), name = "Precipitation (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
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


######
#Africa: Temperature
mod <- list()
obj_af3_t <- array(NA, dim = c(25, 25))
for (i in 1:25) {
  for (j in 1:25) {
    if(sum(mbcaf3 == i & sdbcaf3 == j, na.rm = T) != 0) {
      dat <- points.af3$BIO1[which(mbcaf3 == i & sdbcaf3 == j)]
      mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
      obj_af3_t[26 - j, i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
    }
  }
}

##make raster
r_af3_t <- raster(extent(0, 25, 0, 25), resolution = 1)

##set the values to the obj
r_af3_t <- setValues(r_af3_t, obj_af3_t)
rdf_af3_t <- as.data.frame(r_af3_t, xy=TRUE)
names(rdf_af3_t) <- c("Mean", "SD", "Temp")

color_af3_t <- colorRampPalette(c("darkorchid4","darkorchid","violetred", "deeppink","lightcoral", "lightsalmon"))(round(maxValue(r_af3_t) - minValue(r_af3_t)))
precip_ecospace_af_t <-ggplot(data = rdf_af3_t) +
  geom_raster(mapping=aes(x=Mean, y=SD, fill=Temp)) +
  scale_fill_gradientn(colors = rev(color_af3_t), name = "Temperature (C)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8), labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous(name = "SD", limits = c(0,25), expand = c(0,0), breaks = c(0, 8.3, 16.6, 25), labels = c(0, 0.4, 0.8, 1.2)) +
  scale_x_continuous(name = "Mean", limits = c(0,25), expand = c(0,0), breaks = c(0, 12.5, 25), labels = c(1, 2, 3)) +
  coord_fixed() +
  theme_bw()
precip_ecospace_af_t

######Calculate maximum likelihood for all bins
modmax_af3_t <- array(NA, dim = length(points.af3[ , 1]))
mod <- list()

for (i in 1:length(points.af3[ , 1])) {
  if(!(is.na(points.af3$mbcaf3[i]) | is.na(points.af3$sdbcaf3[i]))) {
    dat <- points.af3$BIO1[which(points.af3$mbcaf3 == points.af3$mbcaf3[i] & points.af3$sdbcaf3 == points.af3$sdbcaf3[i])]
    mod[[i]] <- density(dat[!is.na(dat)], bw = 1)
    modmax_af3_t[i] <- mod[[i]]$x[which.max(mod[[i]]$y)]
  }
}
modmax_af3_t
hist(modmax_af3_t)

#calculate temperature anomaly
anom_temp_af <- points.af3$BIO1 - modmax_af3_t
#Inverse of cubed temp (cube root)
# #R has a bit of trouble with cube root w/ negative numbers, so a new function is needed 
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }
# #Non-cube temp anomaly
# anom_temp_NC_af <- Math.cbrt(anom_temp_af)
# summary(anom_temp_NC_af)
# plot(anom_temp_af, anom_temp_NC_af)

points.af3$tempanom <- as.numeric(anom_temp_af)
summary(points.af3$tempanom, na.rm = T)
points.af3$tempest <- as.numeric(modmax_af3_t)
summary(anom_temp_NC_af, na.rm = T)

#Correlation of observed ~ expected 
tempcor_af <- lm(points.af3$tempest ~ points.af3$BIO1)
tempcor_af
summary(tempcor_af)
cor_tempest_af <- cor.test(points.af3$tempest, points.af3$BIO1, method = "pearson")
#cor=0.3118604 
cor_tempest_af