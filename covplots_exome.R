#!/usr/bin/env Rscript

#Plots coverage across chromosomes
#Run with "Rscript covplots_genome.R [infile.txt] [SampleName1] [xcoverage] [out_directory]" from bash
#©2011 Henrik Stranneheim

#1st argument infile
#2nd argument sample name 1
#3rd argument x scale for coverage
#4th argument out directory

args <- commandArgs(TRUE) #Collects arguments

infile <-args[1] #Collect infile

sampleName1 <- args[2] #Collects sample name 1

xcov <- as.numeric(args[3]) #How much to plot on the x and how many lines on the y axis

od <- args[4] #Collect location of out directory

#Read infile

covChr <- read.table(paste(infile), header=F)
#covChr <-read.table("10-7053.110713_AD035EACXX.1_sorted_genomeCoverageBed.txt", header=T)

#set working directory
setwd( paste(od, sep="") )

##########
#Collect start and stop for entries column 1 & Plot bases vs coverage
##########
pdf(paste( paste( sampleName1, "Chr_all_perBaseVScov", sep="_"),"pdf", sep=".") )

chr <-c(1:22,"X","Y","MT","genome")
chr_start_stop <- data.frame()#(1:26,1:26)
startcounter <- 0
for (pos in 1:length(chr) )  { # For all entries

	 for(i in 1:nrow(covChr)  ) { #Collect entries start_stop coordinates
	
		if (covChr[i,1] == chr[pos]) {
		
			if (startcounter == 0) {
			chr_start_stop[pos,1] <- i #Start
			}
			startcounter <- 1
			chr_start_stop[pos,2] <- i #Stop when end of loop
		}
		
	}
	startcounter <- 0 #Restart startcounter for next entry
}
k <- 0
for (pos in 1:length(chr) )  { #Plot for chr1-24, MT and Genome

start <- chr_start_stop[pos,1] #Find start coordinates
stop <- chr_start_stop[pos,2] #Find stop coordinates

plotcov <- covChr[start:stop,] #Create subset

ktrack <-0 #Enables plot of chr1-22

	if (pos == 1) { #Create empty plot
	max_bases<-max(covChr[0:chr_start_stop[24,2],3]) #Find the y-scale for chr1-24
	plot(c(0:xcov),ylim=c(0, max_bases), col="white", xlab="Coverage (x)", ylab="Nr bases ")
	title("Bases vs coverage")
	legend(x="topright", c("chrX", "chrY", "Chr1-22"), cex=0.8, 
	col=c("red","blue", "black"), pch=c(22,22,-1), lty= c(-1,-1,1) )
	}

	if ( pos == 23 ) {

	points(plotcov[,2], plotcov[,3], col="red", pch=22)
	ktrack <- 1
	k <- k+1
	}
	if ( pos == 24 ) {

	points(plotcov[,2], plotcov[,3], col="blue", pch=22)
	ktrack <- 1
	k <- k+1
	dev.off()
	}
	if ( pos == 25) { #no point in plotting MT
	ktrack <- 1
	}
	if ( pos == 26) {
	pdf(paste( paste( sampleName1, "Genome_perBaseVScov", sep="_"),"pdf", sep=".") )
	plot(plotcov[,2], plotcov[,3], , xlab="Coverage (x)", ylab="Nr bases ")
	title("Bases vs coverage - Genome")
	ktrack <- 1
	dev.off()
	}

	if(ktrack < 1) {
	lines(plotcov[,2], plotcov[,3], col=k, pch=21)
	k <- k+1
	}
}

# ==============================================================================
#		Plot Fraction of bases vs coverage
# ------------------------------------------------------------------------------

k <-0
for (pos in 1:length(chr) )  { # For all entries

start <- chr_start_stop[pos,1] #Find start coordinates
stop <- chr_start_stop[pos,2] #Find stop coordinates

	 for(i in start:stop  ) { #Run for rows per entry

		covChr[i,6] <- sum(as.numeric( covChr[(start+k):stop,3] ) ) #Sum up all bases startign from coverage 0
		covChr[i,7] <- covChr[i,6]/covChr[i,4] #Fraction of total
		k <- k+1
	}
	k <-0
}

pdf(paste( paste( sampleName1, "Chr_all_Fraction_BaseVScov", sep="_"),"pdf", sep=".") )

k <- 0
for (pos in 1:length(chr) )  { #Plot for chr1-24, MT and Genome

start <- chr_start_stop[pos,1] #Find start coordinates
stop <- chr_start_stop[pos,2] #Find stop coordinates

plotcov <- covChr[start:stop,] #Create subset

ktrack <-0 #Enables plot of chr1-22

	if (pos == 1) { #Create empty plot
	plot(c(0:xcov),ylim=c(0, 1), col="white", xlab="Coverage (x)", ylab="Fraction of bases ")
	title("Fraction of Bases vs coverage")
	legend(x="right", c("chrX", "chrY", "Chr1-22", "MT"), cex=0.8, 
	col=c("red","blue", "black", "Black"), pch=c(22,22,-1,21), lty= c(-1,-1,1, -1) )
	}

	if ( pos == 23 ) {

	points(plotcov[,2], plotcov[,7], col="red", pch=22)
	ktrack <- 1
	k <- k+1
	}
	if ( pos == 24 ) {

	points(plotcov[,2], plotcov[,7], col="blue", pch=22)
	ktrack <- 1
	k <- k+1
	}
	if ( pos == 25) { 
	points(plotcov[,2], plotcov[,7], col="black", pch=21)
	ktrack <- 1
	dev.off()
	}
	if ( pos == 26) {
	pdf(paste( paste( sampleName1, "Genome_Fraction_BaseVScov", sep="_"),"pdf", sep=".") )
	plot(plotcov[,2], plotcov[,7], , xlab="Coverage (x)", ylab="Fraction of bases ")
	title("Fraction of Bases vs coverage - Genome")
	ktrack <- 1
	dev.off()
	}

	if(ktrack < 1) {
	lines(plotcov[,2], plotcov[,7], col=k, pch=21)
	k <- k+1
	}
}

