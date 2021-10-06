# --------------------------------------------------
# temperatureModel_summarize.R
# v0 - Initial script from Zack Holden
# 2015-01-23
# Description: This script reads in daily temperature 
# values, queries by date ranges, queries by temp
# threshold to produce grow days and mean temperature.
# ------------------------------

site.dir <- "/Users/jamesonhinkle/Desktop/folders/vcu/landscape genetics/R code/landscape/daily_grids"
flist <-list.files(site.dir, full=T, pattern=".grd")
temp.stack <- stack(flist)


doy <- append(seq(272, 365, 1), seq(1, 271, 1))
doy <- paste("D", doy, sep="")
names(temp.stack) <- doy

# number of grow days function 
gd_fun <- function(x) { sum(length(which(x > 5 & x < 21)))}
gd_mu_fun <- function(x) { mean(x)}

# Get grow days for Nov-May (winter) and June-Oct (summer)
aug1 <- strptime(as.Date("2013-11-01"), format="%Y-%m-%d")$yday
#212
sept30 <- strptime(as.Date("2014-05-01"), format="%Y-%m-%d")$yday
#272

wint.doy <- append(seq(304, 365, 1), seq(1, 120, 1))
wint.doy <- paste("D", wint.doy, sep="")
get.winter <- which(names(temp.stack) %in% wint.doy)
temp.stack.wint <- subset(temp.stack, get.winter)

gd_temp_wint <- calc(temp.stack.wint,fun=gd_fun)
writeRaster(gd_temp_wint, file=paste(site.dir, "mean_temp_largerextent_oct-sept.tif", sep=""), datatype="FLT4S", NAflag=-9999, overwrite=T )