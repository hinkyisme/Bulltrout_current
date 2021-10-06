require(raster)
# create temp based on location, makes map of temperature in river basin
temp <- raster("BT_dams_rstLL.tif")
temp_v <- as.data.frame(values(temp))
colnames(temp_v) <- c("temp")
# get spatial points data frame
sites <- read.csv(file = "latlong_pts.csv", sep = ",", header = TRUE)

# extract points of temperature based off stocking points
require(ggmap)
require(gstudio)
map <- ggmap(population_map(sites, map.type = "terrain", zoom = 10), stratum = site)
map + geom_point(aes(x = Longitude, y = Latitude, color = Temp), data = sites, size = 3)
temp_v <- extract(temp, sites, method = "bilinear")


# make map of temperature at points in model area
# make data frame for locations and temperature
df <- data.frame(x = 1:566, y = 1:566, temp = 1:566)
# make x and y there own vectors
x <- xy.dat$Longitude
y <- xy.dat$Latitude
# place data into data frame
df$x <- x
df$y <- y
df$temp <- temp_v
# make graph of stocking points vs temperature at that point
require(ggplot2)
p <- ggplot(df, aes(x = x, y = y)) + geom_point(aes(color = temp))
p

