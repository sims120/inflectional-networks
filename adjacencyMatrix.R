##################################
##
## This script creates an adjacency matrix of relatedness between inflection classes, counting number of overlapping 
## exponents (given cell) between each class. Takes as input a data frame with id (inflection class labels) in the 
## first column, type frequency of classes in second column, and remaining columns are MSPS (exponents).
##
## Called by a variety of scripts.
##
## Andrea D. Sims
## v. 8-5-19
###################################

classMatrix.fnc = function(data2){
	
#Initialize
	numRows <- nrow(data2)
	adjMatrix_classes <- data.frame(matrix(0, ncol=numRows, nrow=numRows))
	rownames(adjMatrix_classes) = as.character(data2[,1])
	colnames(adjMatrix_classes) = as.character(data2[,1])
	nodes = as.character(data2[,1])

# Fill
	for(i in 1:nrow(data2)){
		row = which(rownames(adjMatrix_classes) == as.character(data2[i,1]))
		for(j in 3:ncol(data2)){
			temp = as.character(data2[i,j])	#Starting loop at 3 hard codes that first two columns are label and typeFreq; all others are cells
			match = which(as.character(data2[,j]) == temp)
			label_match = as.character(data2[match,1])
			if(length(label_match) > 1){
				for(k in 1:length(label_match)){
					if(label_match[k] == as.character(data2[i,1])){	
					} else {
						col = which(colnames(adjMatrix_classes) == label_match[k])
						adjMatrix_classes[row,col] = adjMatrix_classes[row,col] + 1
					}
				}
			}
			if(length(label_match) == 1){
				if(label_match[1] == as.character(data2[i,1])){
				} else {
					cat("Something went wrong here. Match isn't same as input.\n")
				}
			}
		}	
	}
	adjMatrix_classes <- as.matrix(adjMatrix_classes)
#	write.table(adjMatrix_classes, "adjMatrix_classes_temp.txt")
	return(adjMatrix_classes)
}