###########################################################################
## Erin L. Landguth
## MRMCRun.R											
##   Project Description: an example script to run  the function mantel.mc
## Libraries Needed:
##	1. Spatial
##	2. Ecodist
## Project Input:
##	1. Function parameter inputs...
##	2. MRMCFun.R must be loaded into workspace before this script is ran. 								
## Project Steps:										
##	1. User input function parameters
##	2. mantel.mc function run call
##	3. mantel.mc.analysis function run call
##	4. Some plotting commands for mantel.mc.analysis returned values	
## Project Date: December 10, 2008								
###########################################################################

###################
## Load Library	
###################
library(spatial)
library(ecodist)

###############################
## 1. Function parameter inputs	
###############################

batchstring <- ''

# Specify the number of batch runs used
batchno <- 1

# Specify the number of Monte Carlo runs used
mcrunno <- 5
mcstart <- 1

# The total individuals in each file
N <- 1763

# Total runtime length
looptime <- 200
		
# Uncomment one of these styles for nthfile analysis
#nthfile <- c(0,1,11,21,31,41)
nthfile <- c(0,1,2,3,4,5,10,20,30,40,50,100,150,200)
#nthfile <- c(0,1,2,3,4,5,10,20,30,40,50)
#nthfile <- seq(0,looptime,1)
#nthfile <- seq(0,250,50)
#nthfile <- c(1001)

# Genetic distance file location and names and number of
#gddir <- "D:/Code/CDPOP/CDPOP_v1.2.20_20140410/testdata/testmovenos_newscaled_1397588517/"
gddir <- "/Users/jamesonhinkle/Desktop/Max_Dispersal_Sims_GD/Dams"
gdfilename <- 'Gdmatrix'				 				

# Barrier distance file location and names and answer if you want to test barrier model
barrdir <- "D:/projects/CDFISH/Seattle/CDMats/"			
barrfilename <- c("BarrierDmat_406plusPanK.csv")	
barrans <- 'N'				

# Euclidean distance file location and names and answer if you want to test distanc model
distdir <- "/Users/jamesonhinkle/Desktop/Full_CDmat_andXY"		
distfilename <- c("CD_ibd_full.csv")
distans <- 'Y'								

# Landscape distance file location and names and answer if you want to test landscape model
landdir <- "/Users/jamesonhinkle/Desktop/Full_CDmat_andXY"	
landfilename <- c("CD_dams_full.csv")		
landans <- 'Y'					

# Now specify the sample style to use: 
#	known = a known set of indeces to be read from a file
#	random = a random n draw from the total
#	all = run analysis on all points
samplestyle <- 'all'
# If samplestyle = 'random'
sampleno <- 10
# Else if samplestyle = 'known'
sampledir <- "/Users/jamesonhinkle/Desktop/Full_CDmat_andXY/xy_dams_full.csv" 

# Here specify Mantel Test run information 9 of them:
#	Y or N to the specific simple and partial Mantel Tests
#	Define number of permutations for significance test
gentodist.ans <- 'Y'			# Simple genetic ~ distance
gentologdist.ans <- 'N'			# Simple genetic ~ log(distance)
gentobarr.ans <- 'N'			# Simple genetic ~ barrier
gentoland.ans <- 'Y'			# Simple genetic ~ landscape
gentodist.barr.ans <- 'N'		# Partial genetic ~ distance|barrier
gentodist.land.ans <- 'Y'		# Partial genetic ~ distance|landscape
gentobarr.dist.ans <- 'N'		# Partial genetic ~ barrier|distance
gentobarr.land.ans <- 'N'		# Partial genetic ~ barrier|landscape
gentoland.dist.ans <- 'Y'		# Partial genetic ~ landscape|distance
gentoland.barr.ans <- 'N'		# Partial genetic ~ landscape|barrier
mperms <- 99				# Mantel permutations

# Here specifiy Mantel Correlogram to run
#	Y or N to the specific Mantel correlogram
#	Define year to run test at
mgram.gentodist.ans <- 'N'		# Simple genetic ~ distance
mgram.gentobarr.ans <- 'N'		# Simple genetic ~ barrier
mgram.gentoland.ans <- 'N'		# Simple genetic ~ landscape
mgram.gentodist.barr.ans <- 'N'	# Partial genetic ~ distance|barrier
mgram.gentodist.land.ans <- 'N'	# Partial genetic ~ distance|landscape
mgram.gentobarr.dist.ans <- 'N'	# Partial genetic ~ barrier|distance
mgram.gentobarr.land.ans <- 'N'	# Partial genetic ~ barrier|landscape
mgram.gentoland.dist.ans <- 'N'	# Partial genetic ~ landscape|distance
mgram.gentoland.barr.ans <- 'N'	# Partial genetic ~ landscape|barrier
mgramruntime <- 100

###################
## Function code
###################
mantel.mc(batchno,mcrunno,N,nthfile,gddir,gdfilename,barrdir,barrfilename,barrans,
	distdir,distfilename,distans,landdir,landfilename,landans,samplestyle,sampleno,sampledir,gentodist.ans,gentologdist.ans,
	gentobarr.ans,gentoland.ans,gentodist.barr.ans,gentodist.land.ans,gentobarr.dist.ans,gentobarr.land.ans,
	gentoland.dist.ans,gentoland.barr.ans,mperms,mgram.gentodist.ans,mgram.gentobarr.ans,mgram.gentoland.ans,
	mgram.gentodist.barr.ans,mgram.gentodist.land.ans,mgram.gentobarr.dist.ans,mgram.gentobarr.land.ans,
	mgram.gentoland.dist.ans,mgram.gentoland.barr.ans,mgramruntime,batchstring,mcstart)

results <- mantel.mc.analysis(batchno,mcrunno,N,nthfile,gddir,gdfilename,barrdir,barrfilename,barrans,
	distdir,distfilename,distans,landdir,landfilename,landans,samplestyle,sampleno,sampledir,gentodist.ans,gentologdist.ans,
	gentobarr.ans,gentoland.ans,gentodist.barr.ans,gentodist.land.ans,gentobarr.dist.ans,gentobarr.land.ans,
	gentoland.dist.ans,gentoland.barr.ans,mperms,mgram.gentodist.ans,mgram.gentobarr.ans,mgram.gentoland.ans,
	mgram.gentodist.barr.ans,mgram.gentodist.land.ans,mgram.gentobarr.dist.ans,mgram.gentobarr.land.ans,
	mgram.gentoland.dist.ans,mgram.gentoland.barr.ans,mgramruntime)

###################################################
## Plot returned information from mantel.mc.anlysis
###################################################
# If nthfile has a length of 1, then create a vector
if (length(nthfile) == 1)
{
	nthfile <- seq(0,looptime,as.integer(nthfile))	
}

# Simple Mantel genetic~distance
if (gentologdist.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentologdist.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentologdist.mr.left,lty="dashed")
	lines(time,results$gentologdist.mr.right,lty="dashed")
	title("Mantel Test (genetic~log(distance))")
	axis(1,1:(length(nthfile)),nthfile)
}
# Simple Mantel genetic~distance
if (gentodist.ans=='Y')
{  
  x11()
  time <- seq(1,length(nthfile),1)
  plot(time,results$gentodist.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
  lines(time,results$gentodist.mr.left,lty="dashed")
  lines(time,results$gentodist.mr.right,lty="dashed")
  title("Mantel Test (genetic~distance)")
  axis(1,1:(length(nthfile)),nthfile)
}

# Simple Mantel genetic~barrier
if (gentobarr.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentobarr.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentobarr.mr.left,lty="dashed")
	lines(time,results$gentobarr.mr.right,lty="dashed")
	title("Mantel Test (genetic~barrier)")
	axis(1,1:(length(nthfile)),nthfile)
}

# Simple Mantel genetic~landscape
if (gentoland.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentoland.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentoland.mr.left,lty="dashed")
	lines(time,results$gentoland.mr.right,lty="dashed")
	title("Mantel Test (genetic~landscape")
	axis(1,1:(length(nthfile)),nthfile)
}

# Simple Mantel genetic~distance|barrier
if (gentodist.barr.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentodist.barr.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentodist.barr.mr.left,lty="dashed")
	lines(time,results$gentodist.barr.mr.right,lty="dashed")
	title("Mantel Test (genetic~distance|barrier)")
	axis(1,1:(length(nthfile)),nthfile)
}
# Simple Mantel genetic~distance|landscape
if (gentodist.land.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentodist.land.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentodist.land.mr.left,lty="dashed")
	lines(time,results$gentodist.land.mr.right,lty="dashed")
	title("Mantel Test (genetic~distance|landscape)")
	axis(1,1:(length(nthfile)),nthfile)
}
# Simple Mantel genetic~barrier|distance
if (gentobarr.dist.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentobarr.dist.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentobarr.dist.mr.left,lty="dashed")
	lines(time,results$gentobarr.dist.mr.right,lty="dashed")
	title("Mantel Test (genetic~barrier|distance)")
	axis(1,1:(length(nthfile)),nthfile)
}
# Simple Mantel genetic~barrier|landscape
if (gentobarr.land.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentobarr.land.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentobarr.land.mr.left,lty="dashed")
	lines(time,results$gentobarr.land.mr.right,lty="dashed")
	title("Mantel Test (genetic~barrier|landscape)")
	axis(1,1:(length(nthfile)),nthfile)
}
# Simple Mantel genetic~landscape|barrier
if (gentoland.barr.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentoland.barr.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentoland.barr.mr.left,lty="dashed")
	lines(time,results$gentoland.barr.mr.right,lty="dashed")
	title("Mantel Test (genetic~landscape|barrier)")
	axis(1,1:(length(nthfile)),nthfile)
}
# Simple Mantel genetic~landscape|distance
if (gentoland.dist.ans=='Y')
{	
	x11()
	time <- seq(1,length(nthfile),1)
	plot(time,results$gentoland.dist.mr.mean,type="b",ylab="Mantel r",xlab="Time",xaxt="n",ylim=c(-0.01,1.0),las=1)
	lines(time,results$gentoland.dist.mr.left,lty="dashed")
	lines(time,results$gentoland.dist.mr.right,lty="dashed")
	title("Mantel Test (genetic~landscape|distance)")
	axis(1,1:(length(nthfile)),nthfile)
}

# Simple Mgram genetic~distance
if (mgram.gentodist.ans=='Y')
{	
	x11()
	plot(results$gentodist.mg.lag,results$gentodist.mg.mean,type="b",ylab="Mantel r",xlab="Distance",las=1)
	lines(results$gentodist.mg.lag,results$gentodist.mg.left,lty="dashed")
	lines(results$gentodist.mg.lag,results$gentodist.mg.right,lty="dashed")
	title("Mantel Test (genetic~distance)")
	
}

# Simple Mgram genetic~barrier
if (mgram.gentobarr.ans=='Y')
{	
	x11()
	plot(results$gentobarr.mg.lag,results$gentobarr.mg.mean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(results$gentobarr.mg.lag,results$gentobarr.mg.left,lty="dashed")
	lines(results$gentobarr.mg.lag,results$gentobarr.mg.right,lty="dashed")
	title("Mantel Correlogram (genetic~barrier)")
	
}

# Simple Mgram genetic~landscape
if (mgram.gentoland.ans=='Y')
{	
	x11()
	lag <- results$gentoland.mg.lag
	mgmean <- results$gentoland.mg.mean
	mgright <- results$gentoland.mg.right
	mgleft <- results$gentoland.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~landscape)")
	
}

# Simple Mgram genetic~distance|barrier
if (mgram.gentodist.barr.ans=='Y')
{	
	x11()
	lag <- results$gentodist.barr.mg.lag
	mgmean <- results$gentodist.barr.mg.mean
	mgright <- results$gentodist.barr.mg.right
	mgleft <- results$gentodist.barr.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~distance|barrier)")
	
}

# Simple Mgram genetic~distance|landscape
if (mgram.gentodist.land.ans=='Y')
{	
	x11()
	lag <- results$gentodist.land.mg.lag
	mgmean <- results$gentodist.land.mg.mean
	mgright <- results$gentodist.land.mg.right
	mgleft <- results$gentodist.land.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~distance|landscape)")
	
}

# Simple Mgram genetic~barrier|distance
if (mgram.gentobarr.dist.ans=='Y')
{	
	x11()
	lag <- results$gentobarr.dist.mg.lag
	mgmean <- results$gentobarr.dist.mg.mean
	mgright <- results$gentobarr.dist.mg.right
	mgleft <- results$gentobarr.dist.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~barrier|distance)")
	
}
# Simple Mgram genetic~barrier|landscape
if (mgram.gentobarr.land.ans=='Y')
{	
	x11()
	lag <- results$gentobarr.land.mg.lag
	mgmean <- results$gentobarr.land.mg.mean
	mgright <- results$gentobarr.land.mg.right
	mgleft <- results$gentobarr.land.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~barrier|landscape)")
	
}
# Simple Mgram genetic~landscape|distance
if (mgram.gentoland.dist.ans=='Y')
{	
	x11()
	lag <- results$gentoland.dist.mg.lag
	mgmean <- results$gentoland.dist.mg.mean
	mgright <- results$gentoland.dist.mg.right
	mgleft <- results$gentoland.dist.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~landscape|distance)")
	
}
# Simple Mgram genetic~landscape|barrier
if (mgram.gentoland.barr.ans=='Y')
{	
	x11()
	lag <- results$gentoland.barr.mg.lag
	mgmean <- results$gentoland.barr.mg.mean
	mgright <- results$gentoland.barr.mg.right
	mgleft <- results$gentoland.barr.mg.left
	plot(lag,mgmean,type="b",ylab="Mantel r",xlab="Distance",xaxt="n",las=1)
	lines(lag,mgleft,lty="dashed")
	lines(lag,mgright,lty="dashed")
	title("Mantel Correlogram (genetic~landscape|barrier)")
	
}