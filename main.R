###########################
#
# Main function calls for graph-theoretic measures of inflection class systems
#
# Andrea D. Sims
# v. 8-6-19
#
#############################

library(igraph)
library(effects)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(rnetcarto) #For calculating network modularity
library(RColorBrewer)
library(jmuOutlier) #For calculating laplace distribution
library(lawstat) #Goodness-of-fit for laplace
library(scales) #Rescaling data

source("adjacencyMatrix.R")
source("calcDegree.R")
source("complexityIncrease.R")
source("graphs.R")

#################################

#List all files in directory language_plats -- data from which to make calculations and graphs
	all_files <- list.files("language_plats/")

#Create directories into which to put results
	dir.create("graphs")
#	dir.create("subtract_classes")
	dir.create("network_properties")
#	dir.create("entropy")

##################################

#Calling major functions for network property calculations and graphing for individual languages

main.fnc = function(all_files){

##Calculate network properties and produce graphs for individual languages
	for(j in 1:length(all_files)){
#	for(j in 1:1){

		cat("Starting", all_files[j], "\n")

#Read in plat
		file_in <- paste0("language_plats/",all_files[j])
		plat <- read.delim(file_in, header=T, sep="\t")
		
		plat$typeFreq <- as.numeric(as.character(plat$typeFreq))

		language = unlist(strsplit(all_files[j], "_.*"))
		language = unlist(strsplit(language, ".txt"))

#Check for duplicate rows (classes that are identical) and combine rows where found
		plat <- checkDuplicates.fnc(plat)
		
#		cat(nrow(plat))
		
#Check for non-unique inflection class labels
		plat <- checkLabels.fnc(plat)
		
#Create folders in which to put graphs and entropy difference calc files (will not overwrite if already exists)
		dir.create(paste0("graphs/",tolower(language)))
		
#Identify "entropies file" -- pre-calculated numbers for entropy of a full inflectional system, with no classes removed
#		entropies_file <- paste0("entropy/", unlist(strsplit(all_files[j], "\\.txt")), "_entropies.txt")

#Calculate adjacency matrix for the plat
		adjMatrix <- classMatrix.fnc(plat)	
	
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
		networkProperties <- calculateDegree.fnc(adjMatrix, plat)
		temp <- unlist(strsplit(file_in, ".*/"))
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

##Calculate difference between the average conditional entropy of full system and average conditional entropy with one class removed,
##iterating through classes and writing out each version of the system to a separate "drop class" file.
##Run this only if there is not already a set of "drop class" files in the folder subtract_classes.
##Takes a long time.
#		complexityIncrease.fnc(file_in, entropies_file)

#Make scatterplots to visualize contribution of network properties to the complexity of the system as a whole
#		complexityScatters.fnc(file_in, entropies_file, adjMatrix, networkProperties)

		cat("Done with", all_files[j],"\n")
	}
}

##############################

##Produce graphs for COMPARISONS across languages

comparison.fnc = function(all_files){
	all_in <- paste0("language_plats/", all_files)

#Plot global shortest path length x global clustering coefficient, plotting real languages against Monte Carlo simulations
	graphMonteCarlo.fnc(all_in)
	
#Graph t-values for different network properties as predictors for system complexity, with each language as a point on the graph
	graphTVals.fnc()
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
