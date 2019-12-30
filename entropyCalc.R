##################################
##################################
##
## INFLECTIONAL NETWORKS: RESOURCES FOR GRAPH-THEORETIC ANALYSIS OF LINGUISTIC MORPHOLOGY
##
## ANDREA D. SIMS
## v. 1.0.1
## version date: 30 December 2019
##
## For how to cite, see the accompanying README.txt
##
##################################
##################################
##
## This script takes as input a tab-delimited file for which:
## each row represents an inflection class,
## - the first column designates the class
## - the second column is the type frequency of the class
## - all subsequent columns (3 to n) represent morphosyntactic property sets (named in any way)
##    and specify the exponent which realizes the property set in that class (language plat)
## 
## It outputs various conditional entropy measures based on pairwise comparison of forms.
## 
#####################################
#####################################

## Read in a TXT file of inflectional formatives and output a TXT file of 
## calculations run over pairs of paradigm cells (A and B, realized by forms a and b): type frequency of a,
## type frequency of b, joint frequency of ab, probability of a, probability of b, joint probability,
## probability of a given b, entropy of the cell A, and conditional entropy of the cell A given cell B.
##
## 'weighted' is boolean. Should the calculations be weighted by type frequency (i.e. by the number of lexical items in a class) 
## or unweighted (i.e. the type frequency column of the output reflects the number of unique classes in 
## which the relevant form occurs, but not the number of words in the class)?
##
## WARNING: Does not check for uniqueness of rows. Non-unique rows will screw up the calculations. (But check of duplicates is 
## done within main.R before calling this function, by checkDuplicates.fnc)

doCalc.fnc = function(forms, language, weighted){

	calcProb.fnc(forms, weighted)
#Read in output of calcProb.fnc()
	pairwiseTable = read.delim("temp.txt", sep="\t", header=T)
	output <- calcEntropy.fnc(pairwiseTable)
	
	if(weighted == T){
		filename = paste0("entropy/", language, "/", language, "_entropy_weighted_Rcalc.txt")
		write.table(output, file=filename, append=F, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")

#Identify unique combinations of A and B
		temp<-aggregate(data=output, FUN=mean, entropy_AgB ~ msps_A+msps_B)
#Calculate mean over unique combos
		mean<- mean(round(temp$entropy_AgB, digits=9))
		cat("Weighted mean is", mean, "\n")
		
	}
	if(weighted == F){
		filename = paste0("entropy/", language, "/", language, "_entropy_nonweighted_Rcalc.txt")
		write.table(output, file=filename, append=F, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
		
		temp<-aggregate(data=output, FUN=mean, entropy_AgB ~ msps_A+msps_B)
		mean<- mean(round(temp$entropy_AgB, digits=9))
		cat("Unweighted mean is", mean, "\n")
	}
	return(output)	
}

######################################################################

## Calculate pairwise probabilities and write out to a file temp.txt. Called by doCalc.fnc().

calcProb.fnc = function(forms, weighted){

# initialize a data frame into which calculations will be dumped
	typeFreqTable<-data.frame(msps_A=character(), msps_B=character(), form_a=character(), form_b=character(), typeFreq_a=numeric(), typeFreq_b=numeric(), typeFreq_ab=numeric(), total_lexemes=numeric(), total_classes=numeric(), prob_a = numeric(), prob_b=numeric(), joint_prob_ab=numeric(), prob_a_given_b=numeric())

# initialize an output file	
write.table(typeFreqTable, file="temp.txt", quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
	
#Count number of lexemes in dataset
	total_lexemes<-sum(as.numeric(forms[,2])) #Type frequency is in column 2
#Count number of classes in dataset
	total_classes<-nrow(forms)
				
# find all unique forms of msps A and msps B
	colNames_temp<-as.character(names(forms))

# going through all columns but the first two (where the class labels and type frequencies are)
	for(a in 3:ncol(forms)){
		tempA<-unique(as.character(forms[ ,a] ))
		for (b in 3:ncol(forms)){
			tempB<-unique(as.character(forms[ ,b]))

#Don't count ones where a and b represent the same msps
			if(a == b){
			}
			else{
#count joint type frequency of all (a,b) where a != b
				for(i in 1:length(tempA)){
					subsetA<-subset(forms, as.character(forms[ ,a])==tempA[i]) #all rows with form a
					for(k in 1:length(tempB)){
						subsetAB<-subset(subsetA, as.character(subsetA[ , b])==tempB[k]) #joint distrib of a and b
#Only count ones where forms a and b co-occur in same lexeme
						if(nrow(subsetAB)>0){
## Version with type frequency weighting
							if(weighted == T){
								typeFreq_ab<-sum(as.numeric(subsetAB[,2])) #Type frequency is in column 2
#count individual type frequencies of forms a and b
								typeFreq_a<-sum(as.numeric(as.character(subsetA[,2]))) #Type frequency is in column 2
								subsetB<-subset(forms, as.character(forms[ ,b])==tempB[k]) #all rows with form b
								typeFreq_b<-sum(as.numeric(as.character(subsetB[,2]))) #Type frequency is in column 2
#Calculate individual form probabilities: p(a|A), p(b|B)						
#p(a|A)
								prob_a<-typeFreq_a/total_lexemes
#p(b|B)
								prob_b<-typeFreq_b/total_lexemes
#p(a,b)
								joint_prob_ab<-typeFreq_ab/total_lexemes
							}
## Version without type frequency weighting
							if(weighted == F){
								typeFreq_ab<-nrow(subsetAB) #Number of inflection classes in which a,b co-occur
#count individual type frequencies of forms a and b
##This is the part that is changed to take out type frequency weighting
								typeFreq_a<-nrow(subsetA)
								subsetB<-subset(forms, as.character(forms[ ,b])==tempB[k]) #all rows with form b
								typeFreq_b<-nrow(subsetB)
#Calculate individual form probabilities: p(a|A), p(b|B)						
#p(a|A)
								prob_a<-typeFreq_a/total_classes
#p(b|B)
								prob_b<-typeFreq_b/total_classes
#p(a,b)
								joint_prob_ab<-typeFreq_ab/total_classes
							}	

#p(a|b)
							prob_a_given_b<-joint_prob_ab/prob_b

#Output to typeFreqTable
							msps_A<-colNames_temp[a]
							msps_B<-colNames_temp[b]

							form_a<-tempA[i]
							form_b<-tempB[k]

							row = cbind(msps_A, msps_B, form_a, form_b, typeFreq_a, typeFreq_b, typeFreq_ab, total_lexemes, total_classes, prob_a, prob_b, joint_prob_ab, prob_a_given_b)
						
							write.table(row, file="temp.txt", append=T, quote=F, sep="\t", row.names=F, col.names=F, fileEncoding = "UTF-8")
						
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

## Add entropy calculations and return them (called by doCalc.fnc())

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
		subset_A<-subset(typeFreqTable, as.character(typeFreqTable$msps_A) == msps_A[i])
#Rows containing the unique forms within that MSPS
		forms_temp<-match(unique(as.character(subset_A$form_a)), subset_A$form_a)
#Hand coding of column number
		prob_temp<-as.numeric(as.character(subset_A[forms_temp,10]))
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

#Copy probabilities from table
				prob_aGb<-as.numeric(as.character(subsetAB$prob_a_given_b))
				prob_ab<-as.numeric(as.character(subsetAB$joint_prob_ab))
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

##########################################################

## Function that reads in files of entropy calculations output by doCalc.fnc and calculates the min, mean and max entropy values

get_cond_ent.fnc = function(filename){
	df <- read.delim(file=filename, sep="\t")
	ent_uncond<-aggregate(data=df, FUN=mean, entropy_A ~ msps_A+msps_B)
	mean_ent_A<-mean(ent_uncond$entropy_A)
	ent_cond<-aggregate(data=df, FUN=mean, entropy_AgB ~ msps_A+msps_B)
	mean_ent_AgB<-mean(ent_cond$entropy_AgB)
	MI<- mean_ent_A - mean_ent_AgB
	if (length(grep("nonweighted",filename)) == 1){
		weighted <- 0
		}
	else{ weighted <- 1
	}

	distribution<-cbind(filename, weighted, mean_ent_A, mean_ent_AgB, MI)
	return(distribution)
}


