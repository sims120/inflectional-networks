##############################################################################
##############################################################################
####
#### script to get entropy, conditional entropy and mutual inforamation from output files
####
#### Jeff Parker, June 2015;  modified Jan 2015 for other languages
####
##############################################################################
##############################################################################

#install.packages("stringr", repos='http://cran.us.r-project.org')
#library("stringr")

#############################
####
####  reading in data from output files, creating data frames to generate graphs from
####
#############################

#list.files(getwd())

## list of filenames; specified manually
#filenames<-c("2_pairwise_numbers/nuer_entropy_nonweighted_Rcalc.txt", "2_pairwise_numbers/nuer_entropy_weighted_Rcalc.txt")


#df <- read.delim(file=filenames[1], sep="\t")
#ent_uncond<-aggregate(data=df, FUN=mean, entropy_A ~ msps_A+msps_B)
#mean_ent_uncond<-mean(ent_uncond$entropy_A)

## function that reads in output files and calculates the min, mean and max entropy values
get_cond_ent.fnc = function(filename){
	df <- read.delim(file=filename, sep="\t")
	ent_uncond<-aggregate(data=df, FUN=mean, entropy_A ~ msps_A+msps_B)
	mean_ent_uncond<-mean(ent_uncond$entropy_A)
	ent_cond<-aggregate(data=df, FUN=mean, entropy_AgB ~ msps_A+msps_B)
	actual_ent<-mean(ent_cond$entropy_AgB)
	MI<- mean_ent_uncond - actual_ent
	if (length(grep("nonweighted",filename)) == 1){
		weighted <- 0
		}
	else{ weighted <- 1
	}

	distribution<-cbind(filename, weighted, mean_ent_uncond, actual_ent, MI)
	return(distribution)
}

## loops through files of simulations and gathers data
#for (i in 1:length(filenames)){
#	temp<-get_cond_ent.fnc(filenames[i])
#	if (i ==1){
#	data<-temp
#	}
#	else{
#	data<-rbind(data, temp)
#	}
#}

#coercing the output to be a data frame
#data<-as.data.frame(data, stringsAsFactors = FALSE)

#write.table(data, file="nuer_entropies.txt", sep="\t", append=F, quote=F, row.names=F, col.names=T)