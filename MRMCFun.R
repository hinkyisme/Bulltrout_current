###########################################################################
## Erin L. Landguth
## MRMCFun.R											
##   Project Description: 
##	1. mantel.mc: Function to run simple and partial mantel test on 
##   	distance, barrier, and genetic matrices extracting information from 
##   	mutliple folders that were created through a batch or Monte Carlo 
##   	process.  Results outputted to MRMCgentoXXXXXX.csv - the correspond-
#	ing Mantel test.
##	2. mantel.mc.analysis: Function to read in MRMCgentoXXX.csv, calculate
##	mean, sd, confidence intervals and plot these results.
## Libraries Needed:
##	1. Spatial
##	2. Ecodist
## Project Input:
##	1. Function parameter inputs... 								
## MRMCRun.R is an example script to run this function 																	
## Project Date: December 10, 2008								
###########################################################################

#################################
## Function code for Mantel tests
#################################
mantel.mc <- function(batchno,mcrunno,N,nthfile,gddir,gdfilename,barrdir,barrfilename,barrans,
	distdir,distfilename,distans,landdir,landfilename,landans,samplestyle,sampleno,sampledir,gentodist.ans,gentologdist.ans,
	gentobarr.ans,gentoland.ans,gentodist.barr.ans,gentodist.land.ans,gentobarr.dist.ans,gentobarr.land.ans,
	gentoland.dist.ans,gentoland.barr.ans,mperms,mgram.gentodist.ans,mgram.gentobarr.ans,mgram.gentoland.ans,
	mgram.gentodist.barr.ans,mgram.gentodist.land.ans,mgram.gentobarr.dist.ans,mgram.gentobarr.land.ans,
	mgram.gentoland.dist.ans,mgram.gentoland.barr.ans,mgramruntime,batchstring,mcstart)
{
	###############################################
	## 1. Read in sample information - random later
	###############################################
	# Known points draw
	if (samplestyle == 'known')
	{
		sampledraw <- read.table(paste(sampledir,sep=""),sep=",",header=TRUE)
		sampledraw <- sampledraw$SelectedID
	}
	# All points draw
	if (samplestyle == 'all')
	{
		sampledraw <- seq(1,N)
	}

	# If nthfile has a length of 1, then create a vector
	#if (length(nthfile) == 1 && nthfile !=0)
	#{
	#	nthfile <- seq(0,looptime,as.integer(nthfile))	
	#}
	
	######################
	## 2. Batch Loop Begin
	######################
	for (i in 1:1)
	{
		
		###########################################
		## 3. Read in Data - Cost Distance Matrices
		###########################################
		# Read in barrier matrix
		if (barrans == 'Y')
		{
			if (length(barrfilename) == 1)
			{
				barrierfile <- read.table(paste(barrdir,barrfilename[1],sep=""),sep=",",header=FALSE)
			}
			if (length(barrfilename) != 1)
			{
				barrierfile <- read.table(paste(barrdir,barrfilename[i],sep=""),sep=",",header=FALSE)
			}
			barrierfile <- as.matrix(barrierfile[,1:N])
		}
		# Read in distance matrix	
		if (distans == 'Y')
		{
			if (length(distfilename) == 1)
			{
				distancefile <- read.table(paste(distdir,distfilename[1],sep=""),sep=",",header=FALSE)
			}
			if (length(distfilename) != 1)
			{
				distancefile <- read.table(paste(distdir,distfilename[i],sep=""),sep=",",header=FALSE)
			}
			distancefile <- as.matrix(distancefile[,1:N])
		}
		# Read in landscape matrix
		if (landans == 'Y')
		{
			if (length(landfilename) == 1)
			{
				landscapefile <- read.table(paste(landdir,landfilename[1],sep=""),sep=",",header=FALSE)
			}
			if (length(landfilename) != 1)
			{
				landscapefile <- read.table(paste(landdir,landfilename[i],sep=""),sep=",",header=FALSE)
			}
			landscapefile <- as.matrix(landscapefile[,1:N])
		}

		# Create batch directory string
		batchfiledir <- paste(batchstring,'batchrun',as.character(batchno),sep="")				
		print(batchfiledir)	
    
		############################
		## 4. Monte Carlo Loop Begin
		############################
		for (j in mcstart:mcrunno)
		{
	
			# Create Monte Carlo directory string
			mcfiledir <- paste('mcrun',as.character(j-1),'/',sep="")
			print(mcfiledir)
      
      #####################################
      ## Sample here
      #####################################
			# Random draw
			if (samplestyle == 'random')
			{
			  sampledraw <- runif(sampleno,1,N)
			}
			# Sample in barrier matrix
			if (barrans == 'Y')
			{
			  barrier <- lower(barrierfile[sampledraw,sampledraw])
			}
			# Sample in distance matrix	
			if (distans == 'Y')
			{
			  distance <- lower(distancefile[sampledraw,sampledraw])
			}
			# Sample in landscape matrix
			if (landans == 'Y')
			{
			  landscape <- lower(landscapefile[sampledraw,sampledraw])
			}			
      
			#####################################
			## 5. Preliminary vector storage work
			#####################################
			# Create empty vectors to append to for mantelr,pval1,llim,ulim: check all cases
			# Simple genetic ~ distance
			if (gentodist.ans=='Y')
			{
				gentodist.mr <- c()
				gentodist.pv1 <- c()
				gentodist.pv2 <- c()
				gentodist.pv3 <- c()
			}
			# Simple genetic ~ log(distance)
			if (gentologdist.ans=='Y')
			{
				gentologdist.mr <- c()
				gentologdist.pv1 <- c()
				gentologdist.pv2 <- c()
				gentologdist.pv3 <- c()
			}
			# Simple genetic ~ barrier
			if (gentobarr.ans=='Y')
			{
				gentobarr.mr <- c()
				gentobarr.pv1 <- c()
				gentobarr.pv2 <- c()
				gentobarr.pv3 <- c()
			}
			# Simple genetic ~ landscape
			if (gentoland.ans=='Y')
			{
				gentoland.mr <- c()
				gentoland.pv1 <- c()
				gentoland.pv2 <- c()
				gentoland.pv3 <- c()
			}
			# Partial genetic ~ distance|barrier
			if (gentodist.barr.ans=='Y')
			{
				gentodist.barr.mr <- c()
				gentodist.barr.pv1 <- c()
				gentodist.barr.pv2 <- c()
				gentodist.barr.pv3 <- c()
			}
			# Partial genetic ~ distance|landscape
			if (gentodist.land.ans=='Y')
			{
				gentodist.land.mr <- c()
				gentodist.land.pv1 <- c()
				gentodist.land.pv2 <- c()
				gentodist.land.pv3 <- c()
			}
			# Partial genetic ~ barrier|landscape
			if (gentobarr.land.ans=='Y')
			{
				gentobarr.land.mr <- c()
				gentobarr.land.pv1 <- c()
				gentobarr.land.pv2 <- c()
				gentobarr.land.pv3 <- c()
			}
			# Partial genetic ~ barrier|distance
			if (gentobarr.dist.ans=='Y')
			{
				gentobarr.dist.mr <- c()
				gentobarr.dist.pv1 <- c()
				gentobarr.dist.pv2 <- c()
				gentobarr.dist.pv3 <- c()
			}
			# Partial genetic ~ landscape|distance
			if (gentoland.dist.ans=='Y')
			{
				gentoland.dist.mr <- c()
				gentoland.dist.pv1 <- c()
				gentoland.dist.pv2 <- c()
				gentoland.dist.pv3 <- c()
			}
			# Partial genetic ~ landscape|barrier
			if (gentoland.barr.ans=='Y')
			{
				gentoland.barr.mr <- c()
				gentoland.barr.pv1 <- c()
				gentoland.barr.pv2 <- c()
				gentoland.barr.pv3 <- c()
			}
	
			###################################################################
			## 6. Mantel Test: 
			##	Simple: genetic ~ distance
			##	Simple: genetic ~ barrier
			##	Partial: genetic ~ distance|barrier
			##	Partial: genetic ~ barrier|distance
			##	Partial: genetic ~ distance|distance
			##	Mantel Correlogram for each of the above or specified.
			###################################################################
			# Start for loop through each Gdmatrix
			for (k in 1:length(nthfile))
			{
				# Read in genetic distance matrix - piece it up for character read
				file1 <- as.character(nthfile[k])		# For specified nthfile
				file2 <- ".csv"
				genetic <- read.table(paste(gddir,batchfiledir,mcfiledir,gdfilename,file1,file2,sep=""),sep=",",header=FALSE)
				genetic <- as.matrix(genetic[,1:N])
				print(nthfile[k])
		
				# Make matrices lower
				genetic <- lower(genetic[sampledraw,sampledraw])
		    
				# Run Mantel appending results to empty vector: distance
				if (gentodist.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~distance, nperm=mperms, nboot=0)
					# Append loop information
					gentodist.mr <- append(gentodist.mr,mantelrun[1])
					gentodist.pv1 <- append(gentodist.pv1,mantelrun[2])
					gentodist.pv2 <- append(gentodist.pv2,mantelrun[3])
					gentodist.pv3 <- append(gentodist.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentodist.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentodist.mgram <- mgram(genetic,distance,nperm=mperms)
						}
					}
				}
				# Run Mantel appending results to empty vector: distance
				if (gentologdist.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~log(distance), nperm=mperms, nboot=0)
					# Append loop information
					gentologdist.mr <- append(gentologdist.mr,mantelrun[1])
					gentologdist.pv1 <- append(gentologdist.pv1,mantelrun[2])
					gentologdist.pv2 <- append(gentologdist.pv2,mantelrun[3])
					gentologdist.pv3 <- append(gentologdist.pv3,mantelrun[4])
				}

				# Run Mantel appending results to empty vector: barrier
				if (gentobarr.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~barrier, nperm=mperms, nboot=0)
					# Append loop information
					gentobarr.mr <- append(gentobarr.mr,mantelrun[1])
					gentobarr.pv1 <- append(gentobarr.pv1,mantelrun[2])
					gentobarr.pv2 <- append(gentobarr.pv2,mantelrun[3])
					gentobarr.pv3 <- append(gentobarr.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentobarr.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentobarr.mgram <- mgram(genetic,barrier,nperm=mperms)
						}
					}
				}

				# Run Mantel appending results to empty vector: landscape
				if (gentoland.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~landscape, nperm=mperms, nboot=0)
					# Append loop information
					gentoland.mr <- append(gentoland.mr,mantelrun[1])
					gentoland.pv1 <- append(gentoland.pv1,mantelrun[2])
					gentoland.pv2 <- append(gentoland.pv2,mantelrun[3])
					gentoland.pv3 <- append(gentoland.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentoland.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentoland.mgram <- mgram(genetic,landscape,nperm=mperms)
						}
					}
				}

				# Run Mantel appending results to empty vector: distance|barrier
				if (gentodist.barr.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~distance+barrier, nperm=mperms, nboot=0)
					# Append loop information
					gentodist.barr.mr <- append(gentodist.barr.mr,mantelrun[1])
					gentodist.barr.pv1 <- append(gentodist.barr.pv1,mantelrun[2])
					gentodist.barr.pv2 <- append(gentodist.barr.pv2,mantelrun[3])
					gentodist.barr.pv3 <- append(gentodist.barr.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentodist.barr.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentodist.barr.mgram <- mgram(genetic,distance+barrier,nperm=mperms)
						}
					}
				}
		
				# Run Mantel appending results to empty vector: distance|landscape
				if (gentodist.land.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~distance+landscape, nperm=mperms, nboot=0)
					# Append loop information
					gentodist.land.mr <- append(gentodist.land.mr,mantelrun[1])
					gentodist.land.pv1 <- append(gentodist.land.pv1,mantelrun[2])
					gentodist.land.pv2 <- append(gentodist.land.pv2,mantelrun[3])
					gentodist.land.pv3 <- append(gentodist.land.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentodist.land.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentodist.land.mgram <- mgram(genetic,distance+landscape,nperm=mperms)
						}
					}
				}
		
				# Run Mantel appending results to empty vector: barrier|distance
				if (gentobarr.dist.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~barrier+distance, nperm=mperms, nboot=0)
					# Append loop information
					gentobarr.dist.mr <- append(gentobarr.dist.mr,mantelrun[1])
					gentobarr.dist.pv1 <- append(gentobarr.dist.pv1,mantelrun[2])
					gentobarr.dist.pv2 <- append(gentobarr.dist.pv2,mantelrun[3])
					gentobarr.dist.pv3 <- append(gentobarr.dist.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentobarr.dist.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentobarr.dist.mgram <- mgram(genetic,barrier+distance,nperm=mperms)
						}
					}
				}

				# Run Mantel appending results to empty vector: barrier|landscape
				if (gentobarr.land.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~barrier+landscape, nperm=mperms, nboot=0)
					# Append loop information
					gentobarr.land.mr <- append(gentobarr.land.mr,mantelrun[1])
					gentobarr.land.pv1 <- append(gentobarr.land.pv1,mantelrun[2])
					gentobarr.land.pv2 <- append(gentobarr.land.pv2,mantelrun[3])
					gentobarr.land.pv3 <- append(gentobarr.land.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentobarr.land.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentobarr.land.mgram <- mgram(genetic,barrier+landscape,nperm=mperms)
						}
					}
				}

				# Run Mantel appending results to empty vector: landscape|distance
				if (gentoland.dist.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~landscape+distance, nperm=mperms, nboot=0)
					# Append loop information
					gentoland.dist.mr <- append(gentoland.dist.mr,mantelrun[1])
					gentoland.dist.pv1 <- append(gentoland.dist.pv1,mantelrun[2])
					gentoland.dist.pv2 <- append(gentoland.dist.pv2,mantelrun[3])
					gentoland.dist.pv3 <- append(gentoland.dist.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentoland.dist.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentoland.dist.mgram <- mgram(genetic,landscape+distance,nperm=mperms)
						}
					}
				}

				# Run Mantel appending results to empty vector: landscape|barrier
				if (gentoland.barr.ans == 'Y')
				{	
					mantelrun <- mantel(genetic~landscape+barrier, nperm=mperms, nboot=0)
					# Append loop information
					gentoland.barr.mr <- append(gentoland.barr.mr,mantelrun[1])
					gentoland.barr.pv1 <- append(gentoland.barr.pv1,mantelrun[2])
					gentoland.barr.pv2 <- append(gentoland.barr.pv2,mantelrun[3])
					gentoland.barr.pv3 <- append(gentoland.barr.pv3,mantelrun[4])
					# Check for correlogram
					if (mgram.gentoland.barr.ans == 'Y')
					{
						if (nthfile[k] == mgramruntime)
						{
							gentoland.barr.mgram <- mgram(genetic,landscape+barrier,nperm=mperms)
						}
					}

				}# Last Mantel Run End
			
			}# Mantel Loop End

			#######################################
			## 7. Output mantel information to file
			#######################################				
			## Transpose and data.frame vectors and then print and write to file
			if (gentodist.ans=='Y')
			{
				gentodist.mr <- t(data.frame(gentodist.mr))
				gentodist.pv1 <- t(data.frame(gentodist.pv1))
				gentodist.pv2 <- t(data.frame(gentodist.pv2))
				gentodist.pv3 <- t(data.frame(gentodist.pv3))
				# File name
				fileoutputname <- "MRMCgentodist.csv"
				write.table(gentodist.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentodist.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentodist.csv"
					write.table(t(data.frame(gentodist.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentodist.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentodist.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}			
			}
			## Transpose and data.frame vectors and then print and write to file
			if (gentologdist.ans=='Y')
			{
				gentologdist.mr <- t(data.frame(gentologdist.mr))
				gentologdist.pv1 <- t(data.frame(gentologdist.pv1))
				gentologdist.pv2 <- t(data.frame(gentologdist.pv2))
				gentologdist.pv3 <- t(data.frame(gentologdist.pv3))
				# File name
				fileoutputname <- "MRMCgentologdist.csv"
				write.table(gentologdist.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentologdist.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentologdist.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentologdist.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)			
			}
			# Simple genetic ~ barrier
			if (gentobarr.ans=='Y')
			{
				gentobarr.mr <- t(data.frame(gentobarr.mr))
				gentobarr.pv1 <- t(data.frame(gentobarr.pv1))
				gentobarr.pv2 <- t(data.frame(gentobarr.pv2))
				gentobarr.pv3 <- t(data.frame(gentobarr.pv3))
				# File name
				fileoutputname <- "MRMCgentobarr.csv"
				write.table(gentobarr.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentobarr.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentobarr.csv"
					write.table(t(data.frame(gentobarr.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentobarr.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentobarr.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Simple genetic ~ landscape
			if (gentoland.ans=='Y')
			{
				gentoland.mr <- t(data.frame(gentoland.mr))
				gentoland.pv1 <- t(data.frame(gentoland.pv1))
				gentoland.pv2 <- t(data.frame(gentoland.pv2))
				gentoland.pv3 <- t(data.frame(gentoland.pv3))
				# File name
				fileoutputname <- "MRMCgentoland.csv"
				write.table(gentoland.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentoland.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentoland.csv"
					write.table(t(data.frame(gentoland.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentoland.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentoland.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Partial genetic ~ distance|barrier
			if (gentodist.barr.ans=='Y')
			{
				gentodist.barr.mr <- t(data.frame(gentodist.barr.mr))
				gentodist.barr.pv1 <- t(data.frame(gentodist.barr.pv1))
				gentodist.barr.pv2 <- t(data.frame(gentodist.barr.pv2))
				gentodist.barr.pv3 <- t(data.frame(gentodist.barr.pv3))
				# File name
				fileoutputname <- "MRMCgentodist.barr.csv"
				write.table(gentodist.barr.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.barr.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.barr.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.barr.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentodist.barr.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentodist.barr.csv"
					write.table(t(data.frame(gentodist.barr.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentodist.barr.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentodist.barr.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Partial genetic ~ distance|landscape
			if (gentodist.land.ans=='Y')
			{
				gentodist.land.mr <- t(data.frame(gentodist.land.mr))
				gentodist.land.pv1 <- t(data.frame(gentodist.land.pv1))
				gentodist.land.pv2 <- t(data.frame(gentodist.land.pv2))
				gentodist.land.pv3 <- t(data.frame(gentodist.land.pv3))
				# File name
				fileoutputname <- "MRMCgentodist.land.csv"
				write.table(gentodist.land.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.land.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.land.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentodist.land.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentodist.land.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentodist.land.csv"
					write.table(t(data.frame(gentodist.land.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentodist.land.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentodist.land.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Partial genetic ~ barrier|landscape
			if (gentobarr.land.ans=='Y')
			{
				gentobarr.land.mr <- t(data.frame(gentobarr.land.mr))
				gentobarr.land.pv1 <- t(data.frame(gentobarr.land.pv1))
				gentobarr.land.pv2 <- t(data.frame(gentobarr.land.pv2))
				gentobarr.land.pv3 <- t(data.frame(gentobarr.land.pv3))
				# File name
				fileoutputname <- "MRMCgentobarr.land.csv"
				write.table(gentobarr.land.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.land.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.land.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.land.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentobarr.land.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentobarr.land.csv"
					write.table(t(data.frame(gentobarr.land.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentobarr.land.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentobarr.land.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Partial genetic ~ barrier|distance
			if (gentobarr.dist.ans=='Y')
			{
				gentobarr.dist.mr <- t(data.frame(gentobarr.dist.mr))
				gentobarr.dist.pv1 <- t(data.frame(gentobarr.dist.pv1))
				gentobarr.dist.pv2 <- t(data.frame(gentobarr.dist.pv2))
				gentobarr.dist.pv3 <- t(data.frame(gentobarr.dist.pv3))
				# File name
				fileoutputname <- "MRMCgentobarr.dist.csv"
				write.table(gentobarr.dist.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.dist.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.dist.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentobarr.dist.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentobarr.dist.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentobarr.dist.csv"
					write.table(t(data.frame(gentobarr.dist.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentobarr.dist.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentobarr.dist.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Partial genetic ~ landscape|distance
			if (gentoland.dist.ans=='Y')
			{
				gentoland.dist.mr <- t(data.frame(gentoland.dist.mr))
				gentoland.dist.pv1 <- t(data.frame(gentoland.dist.pv1))
				gentoland.dist.pv2 <- t(data.frame(gentoland.dist.pv2))
				gentoland.dist.pv3 <- t(data.frame(gentoland.dist.pv3))
				# File name
				fileoutputname <- "MRMCgentoland.dist.csv"
				write.table(gentoland.dist.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.dist.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.dist.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.dist.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentoland.dist.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentoland.dist.csv"
					write.table(t(data.frame(gentoland.dist.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentoland.dist.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentoland.dist.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}
			# Partial genetic ~ landscape|barrier
			if (gentoland.barr.ans=='Y')
			{
				gentoland.barr.mr <- t(data.frame(gentoland.barr.mr))
				gentoland.barr.pv1 <- t(data.frame(gentoland.barr.pv1))
				gentoland.barr.pv2 <- t(data.frame(gentoland.barr.pv2))
				gentoland.barr.pv3 <- t(data.frame(gentoland.barr.pv3))
				# File name
				fileoutputname <- "MRMCgentoland.barr.csv"
				write.table(gentoland.barr.mr,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.barr.pv1,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.barr.pv2,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol=",",
					row.names=TRUE,col.names=FALSE)
				write.table(gentoland.barr.pv3,file=paste(gddir,fileoutputname,sep=""),append=TRUE,sep=",",eol="\n",
					row.names=TRUE,col.names=FALSE)
				if (mgram.gentoland.barr.ans=='Y')
				{
					# File folder header
					fileoutputname1 <- "MGramMCgentoland.barr.csv"
					write.table(t(data.frame(gentoland.barr.mgram$mgram[,3])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentoland.barr.mgram$mgram[,1])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol=",",
						row.names=TRUE,col.names=FALSE)
					write.table(t(data.frame(gentoland.barr.mgram$mgram[,4])),file=paste(gddir,fileoutputname1,sep=""),append=TRUE,sep=",",eol="\n",
						row.names=TRUE,col.names=FALSE)
				}
			}

		}# Monte Carlo Loop End		

	}# Batch Loop End

}# Function End

#################################################
## Function code for Mantel analysis and plotting
#################################################
mantel.mc.analysis <- function(batchno,mcrunno,N,nthfile,gddir,gdfilename,barrdir,barrfilename,barrans,
	distdir,distfilename,distans,landdir,landfilename,landans,samplestyle,sampleno,sampledir,gentodist.ans,gentologdist.ans,
	gentobarr.ans,gentoland.ans,gentodist.barr.ans,gentodist.land.ans,gentobarr.dist.ans,gentobarr.land.ans,
	gentoland.dist.ans,gentoland.barr.ans,mperms,mgram.gentodist.ans,mgram.gentobarr.ans,mgram.gentoland.ans,
	mgram.gentodist.barr.ans,mgram.gentodist.land.ans,mgram.gentobarr.dist.ans,mgram.gentobarr.land.ans,
	mgram.gentoland.dist.ans,mgram.gentoland.barr.ans,mgramruntime)
{
	##########################################################################
	## For specified Mantel analysis (must be a file allready created for it):
	##	1. Create storage mean and sd vectors
	##	2. Read in file
	##	3. Loop the file, storing and then take the mean and sd and error.
	##	4. Calculate the SE, and left and right bars
	##	5. Write information back to file...appending mean,sd,error,left,right
	##	6. Plot results if specified.
	##########################################################################
	
	# Storage vectors to be return by function
	gentodist.mr.mean <- c()
	gentodist.mr.left <- c()
	gentodist.mr.right <- c()
	gentodist.mg.mean <- c()
	gentodist.mg.left <- c()
	gentodist.mg.right <- c()
	gentodist.mg.lag <- c()
	gentologdist.mr.mean <- c()
	gentologdist.mr.left <- c()
	gentologdist.mr.right <- c()
	gentobarr.mr.mean <- c()
	gentobarr.mr.left <- c()
	gentobarr.mr.right <- c()
	gentobarr.mg.mean <- c()
	gentobarr.mg.left <- c()
	gentobarr.mg.right <- c()
	gentobarr.mg.lag <- c()
	gentoland.mr.mean <- c()
	gentoland.mr.left <- c()
	gentoland.mr.right <- c()
	gentoland.mg.mean <- c()
	gentoland.mg.left <- c()
	gentoland.mg.right <- c()
	gentoland.mg.lag <- c()

	gentodist.barr.mr.mean <- c()
	gentodist.barr.mr.left <- c()
	gentodist.barr.mr.right <- c()
	gentodist.barr.mg.mean <- c()
	gentodist.barr.mg.left <- c()
	gentodist.barr.mg.right <- c()
	gentodist.barr.mg.lag <- c()
	gentobarr.dist.mr.mean <- c()
	gentobarr.dist.mr.left <- c()
	gentobarr.dist.mr.right <- c()
	gentobarr.dist.mg.mean <- c()
	gentobarr.dist.mg.left <- c()
	gentobarr.dist.mg.right <- c()
	gentobarr.dist.mg.lag <- c()
	gentoland.barr.mr.mean <- c()
	gentoland.barr.mr.left <- c()
	gentoland.barr.mr.right <- c()
	gentoland.barr.mg.mean <- c()
	gentoland.barr.mg.left <- c()
	gentoland.barr.mg.right <- c()
	gentoland.barr.mg.lag <- c()

	gentodist.land.mr.mean <- c()
	gentodist.land.mr.left <- c()
	gentodist.land.mr.right <- c()
	gentodist.land.mg.mean <- c()
	gentodist.land.mg.left <- c()
	gentodist.land.mg.right <- c()
	gentodist.land.mg.lag <- c()
	gentobarr.land.mr.mean <- c()
	gentobarr.land.mr.left <- c()
	gentobarr.land.mr.right <- c()
	gentobarr.land.mg.mean <- c()
	gentobarr.land.mg.left <- c()
	gentobarr.land.mg.right <- c()
	gentobarr.land.mg.lag <- c()
	gentoland.dist.mr.mean <- c()
	gentoland.dist.mr.left <- c()
	gentoland.dist.mr.right <- c()
	gentoland.dist.mg.mean <- c()
	gentoland.dist.mg.left <- c()
	gentoland.dist.mg.right <- c()
	gentoland.dist.mg.lag <- c()	

	# Create mean and sd storage vectors and open file
	# Simple genetic ~ distance
	if (gentodist.ans=='Y')
	{
		# Storage vectors
		gentodist.mr.mean <- c()
		gentodist.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentodist.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentodist.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentodist.mr.mean <- append(gentodist.mr.mean,mean(meantemp))
			gentodist.mr.sd <- append(gentodist.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentodist.mr.error <- qnorm(0.975)*gentodist.mr.sd/sqrt((mcrunno))
		gentodist.mr.left <- gentodist.mr.mean - gentodist.mr.error
		gentodist.mr.right <- gentodist.mr.mean + gentodist.mr.error

		# Write information back to file
		write.table(t(data.frame(gentodist.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentodist.ans == 'Y')
		{
			# Storage vectors
			gentodist.mg.mean <- c()
			gentodist.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentodist.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentodist.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentodist.mg.mean <- append(gentodist.mg.mean,mean(meantemp))
				gentodist.mg.sd <- append(gentodist.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentodist.mg.error <- qnorm(0.975)*gentodist.mg.sd/sqrt((mcrunno))
			gentodist.mg.left <- gentodist.mg.mean - gentodist.mg.error
			gentodist.mg.right <- gentodist.mg.mean + gentodist.mg.error

			# Write information back to file
			write.table(t(data.frame(gentodist.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)

			# Store lag intervals for plotting
			gentodist.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		
		}
	}
	
	# Simple genetic ~ distance
	if (gentologdist.ans=='Y')
	{
		# Storage vectors
		gentologdist.mr.mean <- c()
		gentologdist.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentologdist.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentologdist.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentologdist.mr.mean <- append(gentologdist.mr.mean,mean(meantemp))
			gentologdist.mr.sd <- append(gentologdist.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentologdist.mr.error <- qnorm(0.975)*gentologdist.mr.sd/sqrt((mcrunno))
		gentologdist.mr.left <- gentologdist.mr.mean - gentologdist.mr.error
		gentologdist.mr.right <- gentologdist.mr.mean + gentologdist.mr.error

		# Write information back to file
		write.table(t(data.frame(gentologdist.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentologdist.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentologdist.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentologdist.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentologdist.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	}
	
	# Simple genetic ~ barrier
	if (gentobarr.ans=='Y')
	{
		# Storage vectors
		gentobarr.mr.mean <- c()
		gentobarr.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentobarr.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentobarr.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentobarr.mr.mean <- append(gentobarr.mr.mean,mean(meantemp))
			gentobarr.mr.sd <- append(gentobarr.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentobarr.mr.error <- qnorm(0.975)*gentobarr.mr.sd/sqrt((mcrunno))
		gentobarr.mr.left <- gentobarr.mr.mean - gentobarr.mr.error
		gentobarr.mr.right <- gentobarr.mr.mean + gentobarr.mr.error

		# Write information back to file
		write.table(t(data.frame(gentobarr.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentobarr.ans == 'Y')
		{
			# Storage vectors
			gentobarr.mg.mean <- c()
			gentobarr.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentobarr.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentobarr.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentobarr.mg.mean <- append(gentobarr.mg.mean,mean(meantemp))
				gentobarr.mg.sd <- append(gentobarr.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentobarr.mg.error <- qnorm(0.975)*gentobarr.mg.sd/sqrt((mcrunno))
			gentobarr.mg.left <- gentobarr.mg.mean - gentobarr.mg.error
			gentobarr.mg.right <- gentobarr.mg.mean + gentobarr.mg.error

			# Write information back to file
			write.table(t(data.frame(gentobarr.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentobarr.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}
	}
	# Simple genetic ~ landscape
	if (gentoland.ans=='Y')
	{
		# Storage vectors
		gentoland.mr.mean <- c()
		gentoland.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentoland.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentoland.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentoland.mr.mean <- append(gentoland.mr.mean,mean(meantemp))
			gentoland.mr.sd <- append(gentoland.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentoland.mr.error <- qnorm(0.975)*gentoland.mr.sd/sqrt((mcrunno))
		gentoland.mr.left <- gentoland.mr.mean - gentoland.mr.error
		gentoland.mr.right <- gentoland.mr.mean + gentoland.mr.error

		# Write information back to file
		write.table(t(data.frame(gentoland.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentoland.ans == 'Y')
		{
			# Storage vectors
			gentoland.mg.mean <- c()
			gentoland.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentoland.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentoland.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentoland.mg.mean <- append(gentoland.mg.mean,mean(meantemp))
				gentoland.mg.sd <- append(gentoland.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentoland.mg.error <- qnorm(0.975)*gentoland.mg.sd/sqrt((mcrunno))
			gentoland.mg.left <- gentoland.mg.mean - gentoland.mg.error
			gentoland.mg.right <- gentoland.mg.mean + gentoland.mg.error

			# Write information back to file
			write.table(t(data.frame(gentoland.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentoland.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}
	}
	
	# Partial genetic ~ distance|barrier
	if (gentodist.barr.ans=='Y')
	{
		# Storage vectors
		gentodist.barr.mr.mean <- c()
		gentodist.barr.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentodist.barr.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentodist.barr.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentodist.barr.mr.mean <- append(gentodist.barr.mr.mean,mean(meantemp))
			gentodist.barr.mr.sd <- append(gentodist.barr.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentodist.barr.mr.error <- qnorm(0.975)*gentodist.barr.mr.sd/sqrt((mcrunno))
		gentodist.barr.mr.left <- gentodist.barr.mr.mean - gentodist.barr.mr.error
		gentodist.barr.mr.right <- gentodist.barr.mr.mean + gentodist.barr.mr.error

		# Write information back to file
		write.table(t(data.frame(gentodist.barr.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.barr.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.barr.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.barr.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.barr.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentodist.barr.ans == 'Y')
		{
			# Storage vectors
			gentodist.barr.mg.mean <- c()
			gentodist.barr.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentodist.barr.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentodist.barr.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentodist.barr.mg.mean <- append(gentodist.barr.mg.mean,mean(meantemp))
				gentodist.barr.mg.sd <- append(gentodist.barr.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentodist.barr.mg.error <- qnorm(0.975)*gentodist.barr.mg.sd/sqrt((mcrunno))
			gentodist.barr.mg.left <- gentodist.barr.mg.mean - gentodist.barr.mg.error
			gentodist.barr.mg.right <- gentodist.barr.mg.mean + gentodist.barr.mg.error

			# Write information back to file
			write.table(t(data.frame(gentodist.barr.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.barr.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.barr.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.barr.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.barr.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentodist.barr.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}

	}
	# Partial genetic ~ distance|landscape
	if (gentodist.land.ans=='Y')
	{
		# Storage vectors
		gentodist.land.mr.mean <- c()
		gentodist.land.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentodist.land.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentodist.land.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentodist.land.mr.mean <- append(gentodist.land.mr.mean,mean(meantemp))
			gentodist.land.mr.sd <- append(gentodist.land.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentodist.land.mr.error <- qnorm(0.975)*gentodist.land.mr.sd/sqrt((mcrunno))
		gentodist.land.mr.left <- gentodist.land.mr.mean - gentodist.land.mr.error
		gentodist.land.mr.right <- gentodist.land.mr.mean + gentodist.land.mr.error

		# Write information back to file
		write.table(t(data.frame(gentodist.land.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.land.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.land.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.land.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentodist.land.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentodist.land.ans == 'Y')
		{
			# Storage vectors
			gentodist.land.mg.mean <- c()
			gentodist.land.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentodist.land.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentodist.land.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentodist.land.mg.mean <- append(gentodist.land.mg.mean,mean(meantemp))
				gentodist.land.mg.sd <- append(gentodist.land.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentodist.land.mg.error <- qnorm(0.975)*gentodist.land.mg.sd/sqrt((mcrunno))
			gentodist.land.mg.left <- gentodist.land.mg.mean - gentodist.land.mg.error
			gentodist.land.mg.right <- gentodist.land.mg.mean + gentodist.land.mg.error

			# Write information back to file
			write.table(t(data.frame(gentodist.land.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.land.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.land.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.land.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentodist.land.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentodist.land.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}

	}
	# Partial genetic ~ barrier|landscape
	if (gentobarr.land.ans=='Y')
	{
		# Storage vectors
		gentobarr.land.mr.mean <- c()
		gentobarr.land.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentobarr.land.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentobarr.land.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentobarr.land.mr.mean <- append(gentobarr.land.mr.mean,mean(meantemp))
			gentobarr.land.mr.sd <- append(gentobarr.land.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentobarr.land.mr.error <- qnorm(0.975)*gentobarr.land.mr.sd/sqrt((mcrunno))
		gentobarr.land.mr.left <- gentobarr.land.mr.mean - gentobarr.land.mr.error
		gentobarr.land.mr.right <- gentobarr.land.mr.mean + gentobarr.land.mr.error

		# Write information back to file
		write.table(t(data.frame(gentobarr.land.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.land.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.land.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.land.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.land.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentobarr.land.ans == 'Y')
		{
			# Storage vectors
			gentobarr.land.mg.mean <- c()
			gentobarr.land.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentobarr.land.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentobarr.land.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentobarr.land.mg.mean <- append(gentobarr.land.mg.mean,mean(meantemp))
				gentobarr.land.mg.sd <- append(gentobarr.land.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentobarr.land.mg.error <- qnorm(0.975)*gentobarr.land.mg.sd/sqrt((mcrunno))
			gentobarr.land.mg.left <- gentobarr.land.mg.mean - gentobarr.land.mg.error
			gentobarr.land.mg.right <- gentobarr.land.mg.mean + gentobarr.land.mg.error

			# Write information back to file
			write.table(t(data.frame(gentobarr.land.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.land.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.land.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.land.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.land.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentobarr.land.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}

	}
	# Partial genetic ~ barrier|distance
	if (gentobarr.dist.ans=='Y')
	{
		# Storage vectors
		gentobarr.dist.mr.mean <- c()
		gentobarr.dist.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentobarr.dist.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentobarr.dist.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentobarr.dist.mr.mean <- append(gentobarr.dist.mr.mean,mean(meantemp))
			gentobarr.dist.mr.sd <- append(gentobarr.dist.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentobarr.dist.mr.error <- qnorm(0.975)*gentobarr.dist.mr.sd/sqrt((mcrunno))
		gentobarr.dist.mr.left <- gentobarr.dist.mr.mean - gentobarr.dist.mr.error
		gentobarr.dist.mr.right <- gentobarr.dist.mr.mean + gentobarr.dist.mr.error

		# Write information back to file
		write.table(t(data.frame(gentobarr.dist.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.dist.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.dist.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.dist.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentobarr.dist.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentobarr.dist.ans == 'Y')
		{
			# Storage vectors
			gentobarr.dist.mg.mean <- c()
			gentobarr.dist.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentobarr.dist.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentobarr.dist.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentobarr.dist.mg.mean <- append(gentobarr.dist.mg.mean,mean(meantemp))
				gentobarr.dist.mg.sd <- append(gentobarr.dist.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentobarr.dist.mg.error <- qnorm(0.975)*gentobarr.dist.mg.sd/sqrt((mcrunno))
			gentobarr.dist.mg.left <- gentobarr.dist.mg.mean - gentobarr.dist.mg.error
			gentobarr.dist.mg.right <- gentobarr.dist.mg.mean + gentobarr.dist.mg.error

			# Write information back to file
			write.table(t(data.frame(gentobarr.dist.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.dist.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.dist.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.dist.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentobarr.dist.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentobarr.dist.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}

	}
	# Partial genetic ~ landscape|distance
	if (gentoland.dist.ans=='Y')
	{
		# Storage vectors
		gentoland.dist.mr.mean <- c()
		gentoland.dist.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentoland.dist.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentoland.dist.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentoland.dist.mr.mean <- append(gentoland.dist.mr.mean,mean(meantemp))
			gentoland.dist.mr.sd <- append(gentoland.dist.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentoland.dist.mr.error <- qnorm(0.975)*gentoland.dist.mr.sd/sqrt((mcrunno))
		gentoland.dist.mr.left <- gentoland.dist.mr.mean - gentoland.dist.mr.error
		gentoland.dist.mr.right <- gentoland.dist.mr.mean + gentoland.dist.mr.error

		# Write information back to file
		write.table(t(data.frame(gentoland.dist.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.dist.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.dist.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.dist.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.dist.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentoland.dist.ans == 'Y')
		{
			# Storage vectors
			gentoland.dist.mg.mean <- c()
			gentoland.dist.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentoland.dist.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentoland.dist.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentoland.dist.mg.mean <- append(gentoland.dist.mg.mean,mean(meantemp))
				gentoland.dist.mg.sd <- append(gentoland.dist.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentoland.dist.mg.error <- qnorm(0.975)*gentoland.dist.mg.sd/sqrt((mcrunno))
			gentoland.dist.mg.left <- gentoland.dist.mg.mean - gentoland.dist.mg.error
			gentoland.dist.mg.right <- gentoland.dist.mg.mean + gentoland.dist.mg.error

			# Write information back to file
			write.table(t(data.frame(gentoland.dist.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.dist.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.dist.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.dist.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.dist.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentoland.dist.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}

	}
	# Partial genetic ~ landscape|barrier
	if (gentoland.barr.ans=='Y')
	{
		# Storage vectors
		gentoland.barr.mr.mean <- c()
		gentoland.barr.mr.sd <- c()
		
		# Read in
		filename <- "MRMCgentoland.barr.csv"
		mrmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
		# New filename to write to.
		filename <- "MRMCgentoland.barr.analysis.csv"

		# Loop through the length of the mantel r vectors
		for (i in 1:((length(mrmcoutput)/4)-1))
		{
			# Create temp variable
			meantemp <- c()
						
			# Loop through the number of samples to average over and grab and append
			for (j in 1:(mcrunno))
			{
				meantemp <- append(meantemp,mrmcoutput[j,(i+1)])
			}
		
			# Get the mean and sd
			gentoland.barr.mr.mean <- append(gentoland.barr.mr.mean,mean(meantemp))
			gentoland.barr.mr.sd <- append(gentoland.barr.mr.sd,sd(meantemp))
		}
		
		# Calculate error and left and right error bars
		gentoland.barr.mr.error <- qnorm(0.975)*gentoland.barr.mr.sd/sqrt((mcrunno))
		gentoland.barr.mr.left <- gentoland.barr.mr.mean - gentoland.barr.mr.error
		gentoland.barr.mr.right <- gentoland.barr.mr.mean + gentoland.barr.mr.error

		# Write information back to file
		write.table(t(data.frame(gentoland.barr.mr.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.barr.mr.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.barr.mr.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.barr.mr.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
		write.table(t(data.frame(gentoland.barr.mr.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
			row.names=TRUE,col.names=FALSE)
	
		# Mantel correlogram averaging...
		if (mgram.gentoland.barr.ans == 'Y')
		{
			# Storage vectors
			gentoland.barr.mg.mean <- c()
			gentoland.barr.mg.sd <- c()
		
			# Read in
			filename <- "MgramMCgentoland.barr.csv"
			mgmcoutput <- read.csv(paste(gddir,filename,sep=""),sep=",",header=FALSE)
			# New filename to write to.
			filename <- "MgramMCgentoland.barr.analysis.csv"

			# Loop through the length of the mgram vectors
			for (i in 1:((length(mgmcoutput)/3)-1))
			{
				# Create temp variable
				meantemp <- c()
						
				# Loop through the number of samples to average over and grab and append
				for (j in 1:(mcrunno))
				{
					meantemp <- append(meantemp,mgmcoutput[j,(i+1)])
				}
		
				# Get the mean and sd
				gentoland.barr.mg.mean <- append(gentoland.barr.mg.mean,mean(meantemp))
				gentoland.barr.mg.sd <- append(gentoland.barr.mg.sd,sd(meantemp))
			}
		
			# Calculate error and left and right error bars
			gentoland.barr.mg.error <- qnorm(0.975)*gentoland.barr.mg.sd/sqrt((mcrunno))
			gentoland.barr.mg.left <- gentoland.barr.mg.mean - gentoland.barr.mg.error
			gentoland.barr.mg.right <- gentoland.barr.mg.mean + gentoland.barr.mg.error

			# Write information back to file
			write.table(t(data.frame(gentoland.barr.mg.mean)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.barr.mg.sd)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.barr.mg.error)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.barr.mg.left)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			write.table(t(data.frame(gentoland.barr.mg.right)),file=paste(gddir,filename,sep=""),append=TRUE,sep=",",eol="\n",
				row.names=TRUE,col.names=FALSE)
			
			# Store lag intervals for plotting
			gentoland.barr.mg.lag <- t(data.frame(mgmcoutput[1,((length(mgmcoutput)/3)+2):(((length(mgmcoutput)/3)+2)+((length(mgmcoutput)/3)-2))]))
		}

	}

#Return values
list(gentodist.mr.mean=gentodist.mr.mean,
	gentodist.mr.left=gentodist.mr.left,
	gentodist.mr.right=gentodist.mr.right,
	gentodist.mg.mean=gentodist.mg.mean,
	gentodist.mg.left=gentodist.mg.left,
	gentodist.mg.right=gentodist.mg.right,
	gentodist.mg.lag=gentodist.mg.lag,
	gentobarr.mr.mean=gentobarr.mr.mean,
	gentobarr.mr.left=gentobarr.mr.left,
	gentobarr.mr.right=gentobarr.mr.right,
	gentobarr.mg.mean=gentobarr.mg.mean,
	gentobarr.mg.left=gentobarr.mg.left,
	gentobarr.mg.right=gentobarr.mg.right,
	gentobarr.mg.lag=gentobarr.mg.lag,
	gentoland.mr.mean=gentoland.mr.mean,
	gentoland.mr.left=gentoland.mr.left,
	gentoland.mr.right=gentoland.mr.right,
	gentoland.mg.mean=gentoland.mg.mean,
	gentoland.mg.left=gentoland.mg.left,
	gentoland.mg.right=gentoland.mg.right,
	gentoland.mg.lag=gentoland.mg.lag,
	gentodist.barr.mr.mean=gentodist.barr.mr.mean,
	gentodist.barr.mr.left=gentodist.barr.mr.left,
	gentodist.barr.mr.right=gentodist.barr.mr.right,
	gentodist.barr.mg.mean=gentodist.barr.mg.mean,
	gentodist.barr.mg.left=gentodist.barr.mg.left,
	gentodist.barr.mg.right=gentodist.barr.mg.right,
	gentodist.barr.mg.lag=gentodist.barr.mg.lag,
	gentobarr.dist.mr.mean=gentobarr.dist.mr.mean,
	gentobarr.dist.mr.left=gentobarr.dist.mr.left,
	gentobarr.dist.mr.right=gentobarr.dist.mr.right,
	gentobarr.dist.mg.mean=gentobarr.dist.mg.mean,
	gentobarr.dist.mg.left=gentobarr.dist.mg.left,
	gentobarr.dist.mg.right=gentobarr.dist.mg.right,
	gentobarr.dist.mg.lag=gentobarr.dist.mg.lag,
	gentoland.barr.mr.mean=gentoland.barr.mr.mean,
	gentoland.barr.mr.left=gentoland.barr.mr.left,
	gentoland.barr.mr.right=gentoland.barr.mr.right,
	gentoland.barr.mg.mean=gentoland.barr.mg.mean,
	gentoland.barr.mg.left=gentoland.barr.mg.left,
	gentoland.barr.mg.right=gentoland.barr.mg.right,
	gentoland.barr.mg.lag=gentoland.barr.mg.lag,
	gentodist.land.mr.mean=gentodist.land.mr.mean,
	gentodist.land.mr.left=gentodist.land.mr.left,
	gentodist.land.mr.right=gentodist.land.mr.right,
	gentodist.land.mg.mean=gentodist.land.mg.mean,
	gentodist.land.mg.left=gentodist.land.mg.left,
	gentodist.land.mg.right=gentodist.land.mg.right,
	gentodist.land.mg.lag=gentodist.land.mg.lag,
	gentobarr.land.mr.mean=gentobarr.land.mr.mean,
	gentobarr.land.mr.left=gentobarr.land.mr.left,
	gentobarr.land.mr.right=gentobarr.land.mr.right,
	gentobarr.land.mg.mean=gentobarr.land.mg.mean,
	gentobarr.land.mg.left=gentobarr.land.mg.left,
	gentobarr.land.mg.right=gentobarr.land.mg.right,
	gentobarr.land.mg.lag=gentobarr.land.mg.lag,
	gentoland.dist.mr.mean=gentoland.dist.mr.mean,
	gentoland.dist.mr.left=gentoland.dist.mr.left,
	gentoland.dist.mr.right=gentoland.dist.mr.right,
	gentoland.dist.mg.mean=gentoland.dist.mg.mean,
	gentoland.dist.mg.left=gentoland.dist.mg.left,
	gentoland.dist.mg.right=gentoland.dist.mg.right,
	gentoland.dist.mg.lag=gentoland.dist.mg.lag,
	gentologdist.mr.mean=gentologdist.mr.mean,
	gentologdist.mr.right=gentologdist.mr.right,
	gentologdist.mr.left=gentologdist.mr.left)	
	
}
