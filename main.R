##################################
##################################
##
## INFLECTIONAL NETWORKS: RESOURCES FOR GRAPH-THEORETIC ANALYSIS OF LINGUISTIC MORPHOLOGY
##
## ANDREA D. SIMS
## v. 1.0
## version date: 24 December 2019
##
## For how to cite, see the accompanying README.txt
##
##################################
##################################
##
## This script contains the main function calls for graph-theoretic measures of inflection class systems
##
## network.fnc() takes in a (series of) plat(s) of inflectional exponents, calculates various graph-theoretic measures of each inflectional network, 
## makes network and other graphs, and writes out summary files of network properties. It only calculates measures for each language individually.
##
## comparison.fnc() graphs some limited comparisons across all plats -- e.g. plotting mean shortest path length against global clustering coefficient -- 
## and adds Monte Carlo simulations of each language.
##
## entropy.fnc() makes probability/entropy calculations: probability of exponents (unconditioned and conditioned on one other exponent), 
## unconditioned entropy of a cell, entropy of a cell conditioned on one other cell (in the style of Ackerman et al. 2009 and Ackerman and Malouf 2013). 
## It also calculates "entropy difference" -- the difference between the mean conditional entropy of paradigm cells in a full inflectional system and 
## the mean conditional entropy of paradigm cells of the same system with one class removed. It iterates through all classes in the plat and outputs 
## a summary file for each modified IC system.
##
## Other functions in this script are called by these three primary functions.
##
## Input format: A 'language plat' is .txt file containing a tab-delimited table of inflectional exponents. The script expects the file to be located
## in the folder ./language_plats/ . The functions above will, by default, (try to) perform calculations for *all* files in this directory. An unlimited 
## number of language plats is accommodated, although the usual caveats about run time apply.
## 
## Except for the header row, each row in the input file represents an inflection class. Except for the first two columns, each column represents a 
## paradigm cell (i.e. set of inflectional values). Each [row,column] combination is filled with an inflectional exponent for that paradigm cell in that 
## inflection class. The form should NOT be a full word-form. This script does not do segmentation into stems and affixes (or more accurately, 
## themes and distinguishers in the terminology of Stump and Finkel 2013). The assumption is that appropriate segmentation has already been done and 
## the input language plat reflects the result.
##
## The first column of the input file must contain a label for the inflection class. (These can be non-unique, even all "NA". Uniqueness is checked within the 
## script and random numbers are assigned if non-unique labels are found.) The second column of the input file must contain type frequency counts (if type 
## frequency is irrelevant to the analysis, these can all be 1.) The script assumes that paradigm cells begin in column 3. Any number of columns and rows
## is accommodated. In general, column names are not important (as long as they exist).
##
## An example of the input format is made available with this script, in the github repository. See the accompanying read-me.txt for more documentation.
##
#############################
#############################

library(igraph)
library(effects)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(rnetcarto) #For calculating network modularity
library(RColorBrewer)
library(scales) #Rescaling data
library(jmuOutlier) #For calculating laplace distribution
library(lawstat) #Goodness-of-fit for laplace

source("entropyCalc.R")
source("adjacencyMatrix.R")
source("calcDegree.R")
source("complexityIncrease.R")
source("graphs.R")

#################################
#################################

#List all files in directory language_plats -- data from which to make calculations and graphs
all_files <- list.files("language_plats/")

#Create directories into which to put results
dir.create("graphs")  #For generated graphs
dir.create("subtract_classes")  #For 'drop class' / 'entropy difference' calculations (entropy of individual plat with one class removed) + summaries
dir.create("network_properties")  #For network/node properties (degree, clustering, modularity, path length, betweenness centrality, etc.)
dir.create("entropy")  #For entropy calculations over plats + summaries

##################################
##################################

## Do all calculations and graphing

main.fnc = function(all_files){
	network.fnc(all_files)
	entropy.fnc(all_files)
	comparison.fnc(all_files)
}

##################################

## Make networks/graph-theoretic calculations and graphs.

network.fnc = function(all_files){
	
##Calculate network properties and produce graphs for individual languages
	for(j in 1:length(all_files)){

		cat("Starting network property calculations and graphing for ", all_files[j], "\n")

#Read in plat
		file_in <- paste0("language_plats/", all_files[j])
		plat <- read.delim(file_in, header=T, sep="\t")
		
		
		plat[,2] <- as.numeric(as.character(plat[,2]))  #Type frequency is in column 2

		language = unlist(strsplit(all_files[j], ".txt"))

#Check for duplicate rows (classes that are identical) and combine rows where found
		plat <- checkDuplicates.fnc(plat)
		
#Check for non-unique inflection class labels
		plat <- checkLabels.fnc(plat)

#Calculate adjacency matrix for the plat
		adjMatrix <- classMatrix.fnc(plat)	

#Create directory to put individual language graphs in
		dir.create(paste0("graphs/", language))
	
#Plot an undirected graph
#modularity passes a vector of modules. Not used here, but used by identifyModules.fnc, which also calls plotNetwork.fnc

		modules <- vector(length = nrow(adjMatrix))
		modules[] <- 0
		colorByModule <- F #Should nodes be colored depending on their group? (This only works if groups are passed via modules)
		moduleType <- "noGroups" #This is irrelevant unless colorByModule = T, but is used to name the resulting graph
		trim=T
		plotNetwork.fnc(language, plat, adjMatrix, modules, colorByModule, moduleType, trim)
		trim=F
		plotNetwork.fnc(language, plat, adjMatrix, modules, colorByModule, moduleType, trim)

#Calculate basic network properties by node and write out to table
#Only do this if the plat has more than one class.
		if(nrow(plat) > 1){
			
			networkProperties <- calculateDegree.fnc(adjMatrix, plat)
			temp <- unlist(strsplit(file_in, ".*/"))
			temp <- unlist(strsplit(temp, "language_plats/"))
			if(length(temp) > 1){
				temp <- temp[2]
			}
			file_out <- paste0("network_properties/networkProperties_", temp)
			write.table(networkProperties, file_out, row.names=F, col.names=T, quote=F, sep="\t")

#Plot degree distribution of network and related network properties
			degreeDistrib.fnc(language, plat, adjMatrix, networkProperties)
		
#Identify modules using simulated annealing and calculate modularity (incl. role-to-role connectivity profile-defining numbers).
#Plot network graph with coloring based on modules
			identifyModules.fnc(language, plat, adjMatrix)

#Make scatterplots to visualize relationship between degree, edge weight, and clustering
			networkPropertyScatters.fnc(language, plat, networkProperties)

#Plot distributions of mean shortest path length, clustering coefficient, and betweenness centrality, by node. 
			plotShortestPath.fnc(language, networkProperties)
		
#Make network graph with nodes colored according to betweenness centrality.	
			plotNetworkBetweenness.fnc(language, plat, adjMatrix, networkProperties)
		}else{
			cat("Plat ", language, " has only one class. Network properties not calculated or graphed.\n")
		}

		cat("Done with network property calculations and graphing for ", all_files[j],"\n")
	}
}

###############################

entropy.fnc = function(all_files){

	for(i in 1:length(all_files)){

#Read in plat
		file_in <- paste0("language_plats/",all_files[i])
		plat <- read.delim(file_in, header=T, sep="\t")

		plat[,2] <- as.numeric(as.character(plat[,2])) #Type frequency is in column 2

#Check for duplicate rows and remove any found (combine type frequency counts)
		plat <- checkDuplicates.fnc(plat)

#Identify language name, for use when writing out entropy calc files
		language = unlist(strsplit(all_files[i], ".txt"))
		
		dir.create(paste0("entropy/", language))

#Calculate entropy
# 'weighted' is boolean. Should the calculations be weighted by type frequency (i.e. by the number of lexical items in a class) 
# or unweighted (i.e. the type frequency column of the output reflects the number of unique classes in 
# which the relevant form occurs, but not the number of words in the class).
# Also adds to summary file with averages over pairwise numbers for each language, weighted and unweighted

		cat("Starting weighted entropy calculations for ", language, "\n")
		weighted=T
		calc <- doCalc.fnc(plat, language, weighted)
		calc <- as.data.frame(calc)

## Generate summary files and write out
		
		filename = paste0("entropy/", language, "/", language, "_entropy_weighted_Rcalc.txt") #Only used for column entry in summary file

#Find all *unique* combinations of A and B values (so that more entries for one than another doesn't skew calculation of mean)
		temp<-aggregate(data=calc, FUN=mean, entropy_A ~ msps_A+msps_B)
#Calculate mean entropy of A over unique values
		mean_ent_A <- mean(round(temp$entropy_A, digits=9))

#Calculate mean entropy of A given B over unique values		
		temp<-aggregate(data=calc, FUN=mean, entropy_AgB ~ msps_A+msps_B)
		mean_ent_AgB <- mean(round(temp$entropy_AgB, digits=9))

#Calculate mean mutual information over unique values			
		calc$MI <- as.numeric(as.character(calc$entropy_A)) - as.numeric(as.character(calc$entropy_AgB))
		temp<-aggregate(data=calc, FUN=mean, MI ~ msps_A+msps_B)		
		mean_MI <- mean(round(temp$MI, digits=9))
		
		summary <- cbind(filename, weighted, mean_ent_A, mean_ent_AgB, mean_MI)
		summary_file <- paste0("entropy/", language, "/", language, "_entropy_summary.txt")
		
		write.table(summary, summary_file, append=F, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding="UTF-8")

		cat("Starting unweighted entropy calculations for ", language, "\n")
		weighted=F
		calc <- doCalc.fnc(plat, language, weighted)
		calc <- as.data.frame(calc)
		
		filename = paste0("entropy/", language, "/", language, "_entropy_nonweighted_Rcalc.txt") #Only used for column entry in summary file

#Find all *unique* combinations of A and B values (so that more entries for one than another doesn't skew calculation of mean)
		temp<-aggregate(data=calc, FUN=mean, entropy_A ~ msps_A+msps_B)
#Calculate mean entropy of A over unique values
		mean_ent_A <- mean(round(temp$entropy_A, digits=9))

#Calculate mean entropy of A given B over unique values		
		temp<-aggregate(data=calc, FUN=mean, entropy_AgB ~ msps_A+msps_B)
		mean_ent_AgB <- mean(round(temp$entropy_AgB, digits=9))

#Calculate mean mutual information over unique values			
		calc$MI <- as.numeric(as.character(calc$entropy_A)) - as.numeric(as.character(calc$entropy_AgB))
		temp<-aggregate(data=calc, FUN=mean, MI ~ msps_A+msps_B)		
		mean_MI <- mean(round(temp$MI, digits=9))

		summary <- cbind(filename, weighted, mean_ent_A, mean_ent_AgB, mean_MI)
		
		write.table(summary, summary_file, append=T, quote=F, sep="\t", row.names=F, col.names=F, fileEncoding="UTF-8")

#Determine whether individual "entropy difference" (aka "drop class") files already exist. Calculate only if not.
#BEWARE the problems this will cause if the language plat is modified but not renamed, or if the script fails part of the way through.
		list_subtract_files <- list.files(paste0("subtract_classes/", language))
		drop_temp <- grep(language, list_subtract_files)
		if(length(drop_temp > 0)){
			cat("Drop class files found. Skipping calculation of entropy difference. \n")
		}else{
			cat("Starting calculation of entropy difference.\n")
			entropyDiff.fnc(all_files[i])
		}

#Make scatterplots to visualize contribution of network properties to the complexity of the system (entropy difference) as a whole
#Do only if the language plat has more than one class.
		if(nrow(plat) > 1){
			network_file <- unlist(strsplit(all_files[i], "language_plats/"))
			network_file <- paste0("network_properties/networkProperties_", network_file)
			networkProperties <- read.delim(network_file, header=T, sep="\t")
			adjMatrix <- classMatrix.fnc(plat)
			complexityScatters.fnc(file_in, summary_file, adjMatrix, networkProperties)		
		}else{
			cat("Plat ", language, " has only one class. Not graphing entropy difference.\n")
		}
	}
}
		
##############################

##Produce graphs for COMPARISONS across languages

comparison.fnc = function(all_files){
	all_in <- paste0("language_plats/", all_files)

#Plot global shortest path length x global clustering coefficient, plotting real languages against Monte Carlo simulations
	graphMonteCarlo.fnc(all_in)
	
#Graph t-values for different network properties as predictors for system complexity (entropy difference), with each language as a point on the graph
# graphTVals.fnc works, but is commented out here because it requires an input file to be hand-generated, based on fitted regression models for each plat.
#	graphTVals.fnc()
}

###############################

##This runs only if there is not already a set of "drop class" files in the folder subtract_classes. Takes a long time. Takes in filename for single plat.
##Outputs one file for each inflection class in the plat, in subtract_classes directory.

entropyDiff.fnc = function(plat_file){

#Entropy summary file (mean numbers for language plat) -- output by entropy.fnc()
	language = unlist(strsplit(plat_file, ".txt"))
#Where to put drop class files
	dir.create(paste0("subtract_classes/",language))
	summary_file <- paste0("entropy/", language, "/", language, "_entropy_summary.txt")
	plat_file <- paste0("language_plats/", plat_file)

##Calculate difference between the average conditional entropy of full system and average conditional entropy with one class removed,
##iterating through classes and writing out each version of the system to a separate "drop class" file.
	complexityIncrease.fnc(plat_file, summary_file)
}

###############################

#Check for duplicated nodes -- i.e., two classes (rows) with the exact same set of exponents. If found, delete duplicate row and combine type frequencies.
#Find (only) second instance of duplicated row

checkDuplicates.fnc = function(plat){

	plat_temp<- plat
	duplicates <- which(duplicated(plat[,3:ncol(plat)])) #Doesn't include columns 1 and 2 because these are class name and type frequency count
#Find FIRST instance of duplicated row by reversing order of search
	firsts <- which(duplicated(plat[,3:ncol(plat)], fromLast=T)) 
	if(length(duplicates) > 0){
		cat("Duplicate rows found. Rows:\n")
		for(s in 1:length(duplicates)){
			cat(duplicates[s], "\n")
		}
		cat("are duplicates of rows (not necessarily in order):\n")
		for (u in 1:length(firsts)){
			cat(firsts[u], "\n") 
		}
#		stop("Please fix the input file before continuing.\n")
		cat("Deleting duplicate lines and combining type frequency counts.\n")
		for(v in 1:length(firsts)){
			for(q in 1:length(duplicates)){
				duplicatesRow <- duplicates[q]
				firstsRow <- firsts[v]
				temp <- as.vector(plat[duplicatesRow,3:ncol(plat)] == plat[firstsRow,3:ncol(plat)])
#Any FALSE means non-identical. All TRUE means identical. So if (FALSE %in% temp) returns FALSE (i.e. all positions are TRUE -- identical), else condition --> merge rows.
				if(FALSE %in% temp){
				}else{
#Combine type frequencies
					plat_temp[firstsRow,2] <- as.numeric(as.character(plat_temp[firstsRow,2])) + as.numeric(as.character(plat[duplicatesRow,2]))
				}
			}
#Delete duplicate rows
			plat_temp2 <- plat_temp[-duplicates,]
		}		
		return(plat_temp2)
	}else{
		return(plat)
	}
}
	
#################################

#Check whether inflection class labels (class names) are unique and replace if not.

checkLabels.fnc = function(plat){
	plat_new <- plat
	if(length(unique(as.character(plat[,1]))) != length(as.character(plat[,1]))){
		cat("Inflection class labels are not unique. Substituting random codes.\n")
		plat_new[,1] = c(1:length(plat_new[,1]))
	} else if(length(as.character(plat[,1])) == 1){
		cat("Only one inflection class. Substituting random code.\n")
		plat_new[,1] = "1"
	}
	return(plat_new)
}
