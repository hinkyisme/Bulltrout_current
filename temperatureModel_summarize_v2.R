# --------------------------------------------------
# temperatureModel_summarize.R
# v0 - Initial script from Zack Holden
# 2015-01-23
# Description: This script reads in daily temperature 
# values, queries by date ranges, queries by temp
# threshold to produce grow days and mean temperature. 
# WCT model: winter Oct 2 - May 31, summer June 1 - Oct 1
# BT/brook model: winter Oct 16 - July 31, summer August 1 - Oct 15
# --------------------------------------------------

require(raster)

# -------------
# Initial work
# -------------
# Directory location of raster stack
#site.dir <- "D:/projects/CDmetaPOP/Seattle/Data/Temperature/TempLoggers_ZHoldenPCA/"
site.dir <- "D:/projects/CDmetaPOP/RiverscapeSimProject/data/LargerExtent/daily_grids/"

# Points
#xy.dir <- "D:/projects/CDmetaPOP/Seattle/XY/"
xy.dir <- "D:/projects/CDmetaPOP/RiverscapeSimProject/XY/"
#ptfile = "PatchWCT406_deterministic.csv"
#ptfile = "GridCents1379.csv"
ptfile = "PatchXYcoords_centered3.csv"

# If already in stack
temp.stack <- brick(paste(site.dir, "sullivan_daily_streamtemp_stack.grd", sep=""))

# IF daily grids
site.dir <- "/Users/jamesonhinkle/Desktop/folders/vcu/landscape genetics/R code/landscape/daily_grids"
flist <-list.files(site.dir, full=T, pattern=".grd")
temp.stack <- stack(flist)

doy <- append(seq(274, 365, 1), seq(1,273, 1))
doy <- paste("D", doy, sep="")
names(temp.stack) <- doy


# BT and brook dates
startday.winter <- 289 # Oct 16
endday.winter <- 212 # july 31
outname.winter <- "winter_Oct16-July31"
startday.summer <- 213 # Aug 1
endday.summer <- 288 # Oct 15
outname.summer <- "summer_Aug1-Oct15"

# --------
# Winter
# --------
wint.doy <- append(seq(startday.winter, 365, 1), seq(1,endday.winter, 1))
wint.doy <- paste("D", wint.doy, sep="")

get.winter <- which(names(temp.stack) %in% wint.doy)
temp.stack.wint <- subset(temp.stack, get.winter)

# Mean in for just the date period
date_temp_wint <- calc(temp.stack.wint,fun=mean)
writeRaster(date_temp_wint, file=paste(site.dir, "sullivan_mean_temp_",outname.winter,".tif", sep=""), datatype="FLT4S", NAflag=-9999, overwrite=T )

# --------------
# Summer 
# --------------
summer.doy <- seq(startday.summer, endday.summer, 1)
summer.doy <- paste("D", summer.doy, sep="")

get.summer <- which(names(temp.stack) %in% summer.doy)
temp.stack.summer <- subset(temp.stack, get.summer)

# Mean in for just the date period
date_temp_summer <- calc(temp.stack.summer,fun=mean)
writeRaster(date_temp_summer, file=paste(site.dir, "sullivan_mean_temp_summer_aug1-oct15.tif", sep=""), datatype="FLT4S", NAflag=-9999, overwrite=T )

# -------------------
# Extract points here
# -------------------

# Read in points
sites <- read.csv(file = "PatchXYcoords_centered.csv", sep = ",", header = TRUE)

# Grab X,Y values
xvals <- sites$POINT_X
yvals <- sites$POINT_Y
xy <- cbind(xvals,yvals)

# Get a storage variable to append to
xyout <- xy

# Column headers
columnheaders <- c("X","Y")


# Extract X,Y values
tempvals <- extract(date_temp_wint, xy)
# Append to specific storage variable
xyout <- cbind(xyout,tempvals)
#Append to column name
columnheaders <- append(columnheaders,"MeanTempWinter")

tempvals <- extract(date_temp_summer, xy)
# Append to specific storage variable
xyout <- cbind(xyout,tempvals)
#Append to column name
columnheaders <- append(columnheaders,"MeanTempSummer")

# Write out specific information
write.table(xyout,file="Patch566_TempExtract_MeanTemp.csv",append=FALSE,sep=",",eol="\n",row.names=FALSE,col.names=columnheaders)

# ------------------------
# Plotting
# ------------------------
# Create a site object that is a spatial points data frame with XY coordinates
sites <- SpatialPointsDataFrame(coords=as.data.frame(xy),as.data.frame(xy))
plot(date_temp_wint)
points(sites)
#with(sites, text(sites, labels = row.names(sites), pos = 4))
