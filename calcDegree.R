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
## This script calculates node size (type frequency), degree, edge weighting, and local clustering coefficient 
## measures for a graph and return them
## 
## The functions in this script are called by other scripts
##
###################################
###################################

library(igraph)

######################

calculateDegree.fnc = function(adjMatrix2, plat){
	
	plot <- graph_from_adjacency_matrix(adjMatrix2, mode="upper", weight=T)
	degree = degree(plot, v=V(plot))

#Calculate which connected component node belongs to
	component <- as.vector(components(plot, mode=c("strong"))[[1]])

	meanWeight = vector(length=0)
	summedWeight = vector(length=0)
#Each row in the adjacency matrix is a class.
	for(n in 1:nrow(adjMatrix2)){
		if(sum(as.numeric(as.character(adjMatrix2[n,]))) == 0){
			meanWeight = c(meanWeight, 0)
			summedWeight = c(summedWeight, 0)
		}else{
#Calculate mean weighting across all edges associated with the node (class). Crucially, the order of classes in 
#degree object has to be the same as in adjMatrix.
			meanWeight = c(meanWeight, sum(as.numeric(as.character(adjMatrix2[n,])))/degree[n])
			summedWeight = c(summedWeight, sum(as.numeric(as.character(adjMatrix2[n,]))))
		}
	}

#Calculate local clustering coefficient (i.e. by node)
	clusterCoef_local = transitivity(plot, type = "local", isolates="zero")

#Calculate shortest paths for all pairs of nodes, but only within the same connected component.
#This is effectively the same as setting shortest path length for elements in different components to NA.
	listComponents <- unique(component)
	
	meanShortestPath_local = vector(length=nrow(adjMatrix2))
	
	meanShortestPath_local_noWeight = vector(length=nrow(adjMatrix2))
	
	meanShortestPath_local_reverseWeight = vector(length=nrow(adjMatrix2))
	
	betweennessCentrality = vector(length=nrow(adjMatrix2))
	
	for(p in 1:length(listComponents)){
		temp <- which(component == listComponents[p])
#Skip components containing only one node, since there can be no path length for this component.
		if(length(temp) == 1){
			meanShortestPath_local[temp] <- NA
			meanShortestPath_local_noWeight[temp] <- NA
			meanShortestPath_local_reverseWeight[temp] <- NA
			betweennessCentrality[temp] <- NA
			cat("Component", listComponents[p], "has only one node. Assigning NA for shortest path length and betweenness centrality.\n")
#For components with more than one node...
		}else{
#Define a new adjacency matrix and then graph containing only the nodes that are in the same component
			adjMatrix_temp <- adjMatrix2[temp,temp]
			plot_temp <- graph_from_adjacency_matrix(adjMatrix_temp, mode="upper", weight=T)
#Calculate shortest path lengths between all pairs of nodes, calculate mean, and assign back to vector.
			shortestPath_local = distances(plot_temp, v=V(plot_temp), to = V(plot_temp)) 
#The following one removes edge weight in calculations -- same basis as the calculation over the entire network with mean_distance (see Monte Carlo simulation in graphs.R)
			shortestPath_local_noWeight = distances(plot_temp, v=V(plot_temp), to = V(plot_temp), weights=NA) 

#Calculate reverse weight. The issue here is that in networks, edge weight is usually interpreted as *cost* -- e.g., length of a route. So it is normally directly correlated to distance between nodes. But in my graphs it is a correlate to *closeness*. The goal in path length is to *minimize* weight. In calculating shortest path length, using weighting leads the system to look for paths wiht a bias towards going through minimally similar nodes, when in fact it makes more sense for the path to go through maximally similar nodes. So this calculation inverts edge weight for path length calculations; inverted edge weight is positively correlated to distance/non-similarity. Note: If there are 12 cells in paradigm, but the max overlap (edge weight) is 4, then 0 stays 0, and values 1:4 are inverted to 11:8. So reverseWeight becomes a measure of how many cells nodes do NOT have the same exponence in.  

			plot_temp2 <- reverseWeight.fnc(adjMatrix_temp, plat)
			shortestPath_local_reverseWeight = distances(plot_temp2, v=V(plot_temp2), to = V(plot_temp2))

##Calculate betweenness centrality of each node within connected component -- based on reverse weight graph
			betweennessVals <- as.vector(betweenness(plot_temp2, v=V(plot_temp2), directed=F))
				
#Get vector positions in full graph to assign values back to (required because paths calculated only within component)
			names <- colnames(shortestPath_local)
			rows <- vector(length=0)
			for(t in 1:length(names)){
				rows <- c(rows, which(as.character(plat[,1]) == names[t]))
			}
			
#Assign back betweenness centrality values
			betweennessCentrality[rows] <- betweennessVals
			for(m in 1:nrow(shortestPath_local)){
				row_temp <- rows[m]
#Remove 0 value resulting from node-to-itself calculation before calculating mean of shortest paths.
				vals_temp <- shortestPath_local[m,]
				vals_temp_noWeight <- shortestPath_local_noWeight[m,]
				vals_temp_reverse <- shortestPath_local_reverseWeight[m,]
				keep <- which(vals_temp != 0)
#Calculate mean shortest paths and assign back
				meanShortestPath_local[row_temp] <- mean(vals_temp[keep])
				meanShortestPath_local_noWeight[row_temp] <- mean(vals_temp_noWeight[keep])
				meanShortestPath_local_reverseWeight[row_temp] <- mean(vals_temp_reverse[keep])
			}	
		}
	}
	
#Calculating residuals -- weight residualized on degree and clustering coefficient
	weight_lm = lm(summedWeight ~ degree + clusterCoef_local)
	resid_weight = resid(weight_lm)
	
#Get names of classes
	class <- colnames(adjMatrix2)
	
	typeFreq <- vector(length=length(class))
#Match them to type frequencies listed in the plat
	for(s in 1:length(class)){
		freqRow <- which(as.character(plat[,1]) == class[s])  #column 1 contains class labels
		typeFreq[s] <- plat[freqRow,2]  #Type frequency is in column 2
	}
	
	degree_data = as.data.frame(cbind(class, typeFreq, component, degree, meanWeight, summedWeight, resid_weight, clusterCoef_local, meanShortestPath_local, meanShortestPath_local_noWeight, meanShortestPath_local_reverseWeight, betweennessCentrality))
	return(degree_data)
}

########################

reverseWeight.fnc = function(adjMatrix_temp, plat){
	maxWeight <- ncol(plat)-3 #Minus 2 because first two columns in the plat aren't cells, and one more because a cell can't overlap with itself (no loops)
	origWeight <- c(1:maxWeight)
	invertWeight <- c(maxWeight:1)
	adjMatrix_invert <- adjMatrix_temp

#Note: this for-loop will crash if two rows in a plat are identical, since they will have more than the allowed number of overlapping cells (i.e., all of them)
#Error message: Error in adjMatrix_invert[k, j] <- invertWeight[match] : replacement has length zero
#But that shouldn't happen if running from main.R, because checkDuplicates.fnc() is called first.
	for(k in 1:nrow(adjMatrix_temp)){
		for(j in 1:ncol(adjMatrix_temp)){
			if(adjMatrix_temp[k,j] > 0){
				match <- which(origWeight == adjMatrix_temp[k,j])
				adjMatrix_invert[k,j] <- invertWeight[match]
			}
		}
	}
	plot_temp2 <- graph_from_adjacency_matrix(adjMatrix_invert, mode="upper", weight=T)
	return(plot_temp2)
}

			