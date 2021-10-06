#--------------------------------------------------------------------------------------------------------------------------------------------------
#  Extract_Matrix.R
#    By: Jacob J. Burkhart
#    Date: 26 February 2015
#    Updates: 
# 
# Description: Reads in a matrix and extracts the values for a row and column combination above and below the diagonal. Initially used to extract
#              cost values from an asymmetric cost distance matrix, where costs are different above and below diagonal. For example, it extacts all 
#              values in row "i" up to the diagonal and then all cost values below the diagonal in column "i" and combines them into a new dataset.
#              The new dataset is written to a file that can be used to plot cost gradients realized for each patch across the landscape.
#
#   v2 <- added in plotting functions to view cost gradients b/c ArcMAP was being finicky and fickle...
#---------------------------------------------------------------------------------------------------------------------------------------------------

# Extract the cost values from lower half and upper half, including diagonal, then write to new file
#------------------------------------------------------------------------------------------------------
mat <- read.csv("Patch566_20150225_speedwithbarr_costmat.csv", sep = ",", header=F)
mat <- as.matrix(mat)
mat_new <- matrix(NA,nrow=dim(mat)[1],ncol=(2*dim(mat)[2]))
dim(mat_new)
#str(mat_new)

for (i in 1:ncol(mat)){
  mat_new[,i] <- c(mat[i,1:i], mat[i:ncol(mat),i][-1])
  mat_new[,(ncol(mat)+i)] <- c(mat[1:i,i], mat[i,i:ncol(mat)][-1])
 }

write.csv(mat_new, file = "Patch566_20150225_speedwithbarr_costmat_extracted.csv")


# Plot the rasters with costs symbolized as the colors
#------------------------------------------------------------------------------------------------------
require(raster)
require(rgdal)
require(shapefiles)
require(maptools)
require(rgeos)
require(RColorBrewer)
require(ggplot2)
require(rasterVis)

#import data files
  cd <- read.csv("AsymCost_Extract_26Feb2015_Labeled.csv", header=T)     # imports the labeled extracted costs for the asymmetric matrix w/ XY-coords
  strm <- readShapeSpatial("streamnetwork.shp")                          # imports the streams shapefile
  hill <- raster("stmgrdntcls.tif")                                 # imports hillshade raster for background

#plot and prep data files for further analysis
  #plot(hill, col=gray(0:8 / 8)); plot(strm, add=T, col="blue", lwd=2)
  strm.buff <- gBuffer(strm, width=1500)
  hill.sub <- crop(hill, extent(strm.buff))                                                         # clips hill view extent to the 500m buffer on streams layer
  #plot(hill.sub, col=gray(0:8 / 8), main="hill sub"); plot(strm, add=T, col="blue", lwd=2)          # checks to make sure view exten is correct

  xy.dat <- cd[c("POINT_X","POINT_Y")]                                                              #extracts XY coords for the data 
  cd <- SpatialPointsDataFrame(coords=xy.dat,cd)

#set up which patch you want to call
  patch <- 431
  cd.pts.low <- paste("cd$X",patch,"_lower", sep="")
  cd.pts.up <- paste("cd$X",patch,"_upper", sep="")

#plot the hillshade, streams, and patches for lower and upper half of matrix using traditional plotting
  #par(mfrow=c(1,2))                                                                                                        #sets the plotting window to 1 row by 2 columns when not usig ggplot funcitonality
  
  # Plot lower matrix costs
    plot(strm, add=T, col="blue", lwd=5)                                                                        #adds stream layer
    points(x=cd$POINT_X, y=cd$POINT_Y, col=brewer.pal(dim(cd)[1], "YlGnBu"), add=T, pch=20)#, size=4)           #adds the cost gradient for lower half of matrix
    #points(x=cd$POINT_X, y=cd$POINT_Y, data=cd[,4+patch], col=heat.colors(dim(cd)[1], add=T, pch=20))          #matrix columns are not in numerical order
    #plot(cd[patch,], col="red", pch=3, size=8, add=TRUE)
    
  # Plot upper matrix costs 
    plot(hill.sub, col=gray(0:8 / 8), main=paste("Upper Matrix Cost -", patch), legend=F)                       #plots cropped hillshade extent
    plot(strm, add=T, col="blue", lwd=5)                                                                        #adds stream layer
    points(x=cd$POINT_X, y=cd$POINT_Y, data=noquote(cd.pts.up), col=heat.colors(dim(cd)[1]), add=T, pch=20)     #adds the cost gradient values for upper half of matrix


# plotting using ggpolt and rasterVis packages
  # Lower matrix costs
  p.low <- gplot(hill.sub, maxpixels=((dim(hill.sub)[1])*(dim(hill.sub)[2]))/4)  + geom_tile(aes(fill=value)) + coord_equal() +      #plots cropped hillshade extent
                 scale_fill_gradient(low="black", high="white") +
                 geom_line(strm, aes(size=5, colour="blue")) + 
                 geom_point(cd, fill=noquote(cd.pts.low)) 
  p.low

  #upper matrix costs
  p.up <- gplot(hill.sub, maxpixels=((dim(hill.sub)[1])*(dim(hill.sub)[2]))/4)  + geom_tile(aes(fill=value)) + coord_equal() +      #plots cropped hillshade extent
                scale_fill_gradient(low="black", high="white") +
                geom_line(strm, aes(size=5, colour="blue")) + 
                geom_point(cd, fill=noquote(cd.pts.low)) 
  p.up
  
  grid.arrange(p.low, p.up, ncol=2)                                                                           # plots the maps side by side for comparison