devtools::install_github("mariahm1995/commecometrics", force = TRUE)

# Download the files needed for the assessment (only one time)
options(timeout = 600)
download.file("https://ndownloader.figshare.com/files/55632824", destfile = "data.zip", mode = "wb")
unzip("data.zip")
unlink("data.zip")


library(commecometrics)
data("geoPoints", package = "commecometrics")
data("traits", package = "commecometrics")
data("spRanges", package = "commecometrics")


traitsByPoint <- summarize_traits_by_point(
  points_df = geoPoints,
  trait_df = traits,
  species_polygons = spRanges,
  summary_trait_1 = function(x) mean(x),
  trait_column = "RBL",
  species_name_col = "sci_name",
  continent = TRUE,
  parallel = TRUE
)

inspect_point_species(traitsByPoint)

ecoModel <- ecometric_model(
  points_df = traitsByPoint$points,
  env_var = "precip",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  min_species = 3
)

ecoPlot <- ecometric_space(ecoModel, env_name = "Precipitation (loge mm)")
print(ecoPlot)
head(traitsByPoint$points)

# Run ecometric model
ecoModel <- ecometric_model_qual(
  points_df = traitsByPoint$points,
  category_col = "precip"
)

# Plot trait–environment surface
ecoPlot <- ecometric_space_qual(ecoModel)
print(ecoPlot$probability_maps[[5]])
ls("package:pkg")

points <- read.csv("inputs/S_L_allpointsdata2.csv")
continents <- raster::shapefile("inputs/continents.shp")
traits <- read.csv("inputs/RBL_global_Dec2023.csv", header = TRUE, na = ".", stringsAsFactors = FALSE)
traits$TaxonName <- gsub("_", " ", traits$TaxonName)

geometry <- sf::st_read("inputs/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
geography_carnivora <- geometry[geometry$order_ == "CARNIVORA", ]
rm(geometry)