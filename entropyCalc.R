################################################################
##
## This script takes as input a tab-delimited file for which:
## each row represents an inflection class,
## - the first column designates the class
## - the second column is the type frequency of the class
## - all subsequent columns (3 to n) represent morphosyntactic property sets (named in any way)
##    and specify the exponent which relizes the property set in that class
## 
## It outputs various conditional entropy measures based on pairwise comparison of forms.
## 
## Andrea D. Sims
## v. 7-2-14
## 
## Modified by Jeff Parker to accomodate Russian (and other languages)
## 9-2014
##
## Streamlined by Andrea 4-2016

################################################################
################################################################

## set name of input file
#file_name<-"greek_49classes.txt"

#path = paste("1_matrices_with_typeFreq/", file_name, sep="")

# Set whether calculations should be weighted by type frequency of inflection classes
#weighted = T	

################################################################

# Read in raw data
#forms<-read.delim(path, header=T, sep="\t", stringsAsFactors=FALSE)

## to get the UTF-8 characters to display correctly; 
#forms<-as.matrix(forms)

######################################################################

## This function reads in a TXT file of inflectional formatives and outputs a TXT file of 
## calculations run over pairs of paradigm cells (A and B, realized by forms a and b): type frequency of a,
## type frequency of b, joint frequency of ab, probability of a, probability of b, joint probability,
## probability of a given b, entropy of the cell A, and conditional entropy of the cell A given cell B.
##
## The calculations are either weighted by type frequency (i.e. the type frequency of individual formatives, 
## and probabilities, is weighted by the number of lexical items in a class) or unweighted (i.e. each class
## is weighted equally, so the type frequency column of the output reflects the number of unique classes in 
## which the relevant form occurs). This is a parameter setting for the function call.
##
## WARNING: The script does not check for uniqueness of rows. For calculations to be accurate, each row
## needs to be unique. No duplicate inflection classes. (We will fix this someday...)

doCalc.fnc = function(file, weighted){
	file_name = strsplit(file, ".*/")
	file_name = file_name[[1]][2]
	forms<-read.delim(file, header=T, sep="\t", stringsAsFactors=FALSE)
##removing first column (class index) so input file is correctly formatted
	forms<-forms[,2:ncol(forms)]
#	pairwiseTable <- calcProb.fnc(forms, weighted)
	calcProb.fnc(forms, weighted)
	cat("Done with pairwise probabilities\n")
#Read in output of calcProb.fnc()
	pairwiseTable = read.delim("prob_calcs.txt", sep="\t", header=T)
	output <- calcEntropy.fnc(pairwiseTable)
	cat("Done with entropy calculations\n")
	if(weighted == T){
		write.table(output, file=paste0("2_pairwise_numbers/", strsplit(file_name,".txt"), "_entropy_weighted_Rcalc.txt"), append=F, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
		temp<-aggregate(data=output, FUN=mean, entropy_AgB ~ msps_A+msps_B)
		mean<- mean(round(temp$entropy_AgB, digits=9))
		cat("Weighted mean is", mean, "\n")
	}
	if(weighted == F){
		write.table(output, file=paste0("2_pairwise_numbers/", strsplit(file_name,".txt"), "_entropy_nonweighted_Rcalc.txt"), append=F, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
		temp<-aggregate(data=output, FUN=mean, entropy_AgB ~ msps_A+msps_B)
		mean<- mean(round(temp$entropy_AgB, digits=9))
		cat("Unweighted mean is", mean, "\n")
	}
#	save.image(paste0(getwd(),".RData"))	
}

######################################################################

## Calculate pairwise probabilities and write out to a file prob_calcs.txt. Called by doCalc.fnc().

calcProb.fnc = function(forms, weighted){

# initialize a data frame into which calculations will be dumped
	typeFreqTable<-data.frame(msps_A=character(), msps_B=character(), form_a=character(), form_b=character(), typeFreq_a=numeric(), typeFreq_b=numeric(), typeFreq_ab=numeric(), total_lexemes=numeric(), prob_a = numeric(), prob_b=numeric(), joint_prob_ab=numeric(), prob_a_given_b=numeric())

# initialize an output file	
write.table(typeFreqTable, file="prob_calcs.txt", quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
	

# find all unique forms of msps A and msps B
	colNames_temp<-as.character(names(forms))

# going through all columns but the first (where the type frequencies are)
	for(a in 2:ncol(forms)){
		tempA<-unique(as.character(forms[ , a] ))
		for (b in 2:ncol(forms)){
			tempB<-unique(as.character(forms[ , b]))

#Don't count ones where a and b represent the same msps
			if(a == b){
			}
			else{

#count joint type frequency of all (a,b) where a != b
				for(i in 1:length(tempA)){
					subsetA<-subset(forms, as.character(forms[ , a])==tempA[i]) #all rows with form a
					for(k in 1:length(tempB)){
						subsetAB<-subset(subsetA, as.character(subsetA[ , b])==tempB[k]) #joint distrib of a and b
#Only count ones where forms a and b co-occur in same lexeme
						if(nrow(subsetAB)>0){
## Version with type frequency weighting
							if(weighted == T){
								typeFreq_ab<-sum(as.numeric(subsetAB$typeFreq))

#count individual type frequencies of forms a and b
								typeFreq_a<-sum(as.numeric(subsetA$typeFreq))
								subsetB<-subset(forms, as.character(forms[ , b])==tempB[k]) #all rows with form b
								typeFreq_b<-sum(as.numeric(subsetB$typeFreq))

#total number of lexemes in data set = p(B) and p(A) b/c all lexemes assumed to have all forms
#this is where the script could be modified to take account of different token frequencies
#of different msps
								total_lexemes<-sum(as.numeric(forms$typeFreq))
							}
## Version without type frequency weighting
							if(weighted == F){
								typeFreq_ab<-nrow(subsetAB) #Number of inflection classes in which a,b co-occur

#count individual type frequencies of forms a and b
##This is the part that is changed to take out type frequency weighting
								typeFreq_a<-nrow(subsetA)
								subsetB<-subset(forms, as.character(forms[ , b])==tempB[k]) #all rows with form b
								typeFreq_b<-nrow(subsetB)
							
								total_lexemes<-nrow(forms)
							}	

#Calculate individual form probabilities: p(a|A), p(b|B)						
#p(a|A)
							prob_a<-typeFreq_a/total_lexemes

#p(b|B)
							prob_b<-typeFreq_b/total_lexemes

#p(a,b)
							joint_prob_ab<-typeFreq_ab/total_lexemes

#p(a|b)
							prob_a_given_b<-joint_prob_ab/prob_b

#Output to typeFreqTable
							msps_A<-colNames_temp[a]
							msps_B<-colNames_temp[b]

							form_a<-tempA[i]
							form_b<-tempB[k]

							row = cbind(msps_A, msps_B, form_a, form_b, typeFreq_a, typeFreq_b, typeFreq_ab, total_lexemes, prob_a, prob_b, joint_prob_ab, prob_a_given_b)
						
							write.table(row, file="prob_calcs.txt", append=T, quote=F, sep="\t", row.names=F, col.names=F, fileEncoding = "UTF-8")
						
						}	
					}
				}
				cat(".")
			}
		}
		cat("!")
	}
}


#######################

##Add entropy calculations and return them (called by doCalc.fnc())

calcEntropy.fnc = function(typeFreqTable){

#Add column to the data set
	entropy_A<-vector(length=nrow(typeFreqTable))
	entropy_AgB<-vector(length=nrow(typeFreqTable))
	typeFreqTable<-cbind(typeFreqTable, entropy_A, entropy_AgB)

#H(a|A) -- on average, how predictable is a form a belonging to A, without regard for what the 
#particular form a is. General predictability of msps A forms

	msps_A<-unique(as.character(typeFreqTable$msps_A))
	for(i in 1:length(msps_A)){
#Find all the forms with a given MSPS.
		logic_temp<-as.character(typeFreqTable$msps_A)==msps_A[i]
		subset_A<-subset(typeFreqTable, logic_temp)
#Find all the unique forms within that MSPS
		forms_temp<-match(unique(as.character(subset_A$form_a)), subset_A$form_a)
#Hand coding of column number
		prob_temp<-as.numeric(as.character(subset_A[forms_temp,9]))
		temp<-vector(length=length(prob_temp))
		for(k in 1:length(prob_temp)){
			temp[k]<-(prob_temp[k])*(log2(prob_temp[k]))*-1
		}
		entropyA_temp<-sum(temp)
		row_temp<-which(typeFreqTable$msps_A == msps_A[i])
		for(m in 1:length(row_temp)){
			typeFreqTable$entropy_A[row_temp[m]] <- entropyA_temp
		}
		cat("&")
	}

##H(A|B) -- on average, how predictable is a form 'a' belonging to A, given a form 'b' belonging to B,
##without regard for what the particular forms a and b are. General predictability of msps A, given B

	msps_B<-msps_A
	for(i in 1:length(msps_A)){
		logic_temp<-as.character(typeFreqTable$msps_A)==msps_A[i]
		subset_A<-subset(typeFreqTable, logic_temp)
		for(k in 1:length(msps_B)){

#joint distribution of A and B (ordered)
			logic_temp2<-as.character(subset_A$msps_B)==msps_B[k]
			subsetAB<-subset(subset_A, logic_temp2)

#Discard ones where A and B represent the same msps. There should be none of these,
#but just in case...
			if(as.character(msps_A[i]) == as.character(msps_B[k])){
			}
			else{

#Hand coding of column number
#Copy probabilities from table
				prob_aGb<-as.numeric(as.character(subsetAB[,12]))
				prob_ab<-as.numeric(as.character(subsetAB[,11]))
				tempaGb<-vector(length=length(prob_ab))
#Calculate
				for(m in 1:length(tempaGb)){
					tempaGb[m]<-(prob_ab[m])*(log2(prob_aGb[m]))*-1
				}
						
				entropyAgB_temp<-sum(tempaGb)
#Assign to the table
				rowA_temp<-which(typeFreqTable$msps_A == msps_A[i])
				rowB_temp<-which(typeFreqTable$msps_B == msps_B[k])
				row_temp<-intersect(rowA_temp, rowB_temp)
			
				for(m in 1:length(row_temp)){
					typeFreqTable$entropy_AgB[row_temp[m]] <- entropyAgB_temp
				}
			}
		}
		cat("@")
	}
	return(typeFreqTable)
}




