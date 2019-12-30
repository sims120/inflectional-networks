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
## This script calculates "entropy difference" for a plat -- for each class/node, the difference between the mean conditional entropy 
## of the entire inflection class system (calculated pairwise over paradigm cells) minus the mean conditional entropy of the system with 
## the target class removed, iterating over all classes. It also produces some graphs related to entropy difference.
##
## The functions in this script are called by other scripts, mostly main.R
##
##################################
##################################

library(igraph)
library(ggplot2)
library(ggrepel)
library(effects)
library(gridExtra)

##################################

## See how much each class contributes to the complexity of the system by comparing system complexity of all classes N to
## system complexity of N-1 classes. Interatively over all classes.

#	input = filename for language plat of exponents
#	full = filename for pre-calculated entropy numbers for entire inflectional system -- no dropped classes

complexityIncrease.fnc = function(input, full){
	source("entropyCalc.R")
	weighted = F

	file = read.delim(input, header=T, sep="\t")
	file <- checkDuplicates.fnc(file)
	file[,1] = paste("Class",c(1:nrow(file)), sep="") #Change class id's because they are a pain

#Get file name without extension
		language = unlist(strsplit(input, ".txt"))
		language = unlist(strsplit(language, "language_plats/"))[2]
		
#Run only if the language plat has more than one class.
	if(nrow(file) > 1){
	
		full_summary <- read.delim(full, header=T, sep="\t")
		row <- grep("nonweighted", as.character(full_summary$filename))
		summary <- full_summary[row,]
		full_meanEnt <- as.numeric(as.character(summary$mean_ent_AgB))
	
		for(j in 1:nrow(file)){
#Remove one row (i.e. one inflection class), and calculate entropy over the remainder
			type_freq_dropped = file[j,2] #Type frequencies are in column 2 in input file.
			data = file[-c(j),]
#Output directory
			out_file = paste0("drop", j, "_", language, ".txt")
			out = paste0("subtract_classes/", language, "/", out_file)
#Write out individual inflection class data -- one class removed
			write.table(data, out, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
#Do entropy calculations -- calling script_entropy_calc.R
			output <- doCalc.fnc(data, language, weighted)
#Write out output of doCalc.fnc()
			get_file = paste0("subtract_classes/", language, "/", "drop", j, "_", language, "_entropy.txt")
			write.table(output, get_file, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding="UTF-8")
#Get average numbers (inflectional system entropy) for output of doCalc.fnc() (one dropped class)
			values = get_cond_ent.fnc(get_file)
			values = as.data.frame(values)
			log_type_freq <- log(as.numeric(as.character(type_freq_dropped)))
			ent_difference <- full_meanEnt - as.numeric(as.character(values$mean_ent_AgB))
			values = cbind(type_freq_dropped, log_type_freq, values, ent_difference)
#Write out results to a composite file
			if(j == 1){
				write.table(values, paste("subtract_classes/", language, "/", language, "_dropSummary.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
			}else{
				write.table(values, paste("subtract_classes/", language, "/", language, "_dropSummary.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=F, fileEncoding = "UTF-8", append=T) #append and don't write out column names
			}
		}
	}else{
		cat("Plat ", language, " has only one class. Skipping calculation of entropy difference.\n")
	}
}

########################################

#Graph difference between complexity of full system and complexity with one class removed

#input = filename for table written out as output of complexityIncrease.fnc(), or series of such filenames

complexityIncreaseGraph.fnc = function(input){	

	for(k in 1:length(input)){	

		file_name = unlist(strsplit(input[k], ".txt"))
		file_name = unlist(strsplit(file_name, ".*/"))
		if(length(file_name) > 1){
			file_name = file_name[2]
		}

#Read in "drop class" summary -- summary of conditional entropy of system with one class removed (iterated)
		all_values = read.delim(paste("subtract_classes/", file_name, "/", file_name, "_dropSummary.txt", sep=""), header=T, sep="\t")
		all_values = all_values[order(all_values$log_type_freq),]
		
		all_values = as.data.frame(all_values)
		
		if(length(grep("icelandic", file_name)) > 0){
			trimmed = subset(all_values, all_values$log_type_freq < 4) # Trim high type frequency classes because of data sparsity
		} else if(length(grep("french", file_name)) > 0){
			trimmed = subset(all_values, all_values$log_type_freq < 6) # Trim high type frequency classes because of data sparsity
		} else{
			trimmed = all_values
		}
		linear_fit = lm(trimmed$ent_difference ~ trimmed$log_type_freq)
		pVal = coef(summary(linear_fit))["trimmed$log_type_freq", "Pr(>|t|)"]
		rSquare = summary(linear_fit)$r.squared

#The following for individual graphs of each language

		language <- unlist(strsplit(file_name, "_.*"))

		graph_name = paste0("graphs/",file_name,"/dropClass_", file_name, ".png")
		png(graph_name, width=680, height=680, units="px")
		par(mar = c(6,6,2,2 + 0.1))
			
		plot(all_values$ent_difference ~ all_values$log_type_freq, xlab="Log Type Frequency of Dropped Class", ylab="Entropy Difference", cex=1.25, cex.axis=1.5, cex.lab=2, col="black", pch=21, bg="gray85")
		lines(trimmed$log_type_freq, predict(linear_fit), col="black", lwd=4)
		legend("bottomright", c(paste("p = ", round(pVal, digits = 5), sep=""), paste("R^2 = ", round(rSquare, digits = 3), sep="")), cex=1.5)
		
		dev.off()
	
#		cat("\nDone with", input[k], " entropy diff x type frequency \n")
		
	}
}

###################################

## Create graphs of network properties as predictors of entropy difference. 
## WARNING: This script will run by using a default (full) model, but this is likely to be inappropriate.
## This function needs to be customized in order to fit a model to each data set individually/by hand.

complexityComponentsGraph.fnc = function(file, degree_data){

#Read in complexity increase data
	file_name = unlist(strsplit(file, ".txt"))
	file_name = unlist(strsplit(file_name, "language_plats/"))[2]	
	language = file_name

	all_values2 = read.delim(paste("subtract_classes/", file_name, "/", file_name, "_dropSummary.txt", sep=""), header=T, sep="\t")
	all_values2 = cbind(all_values2, degree_data)
	all_values2 = as.data.frame(all_values2)

# Check that the data sets merged correctly by matching whether the vector of type frequencies from the "drop class" file is the same order 
# as the vector of type frequencies from the "degree data" file.
	for(i in 1:length(all_values2$type_freq_dropped)){
		if(all_values2$type_freq_dropped[i] != all_values2$typeFreq[i]){
			cat("Data sets did not merge correctly!\n")
		}
	}

#Check that the data are paired correctly. Seems to be true on manual check, but this check doesn't work for some reason -- throws false warnings..
#	for(i in 1:nrow(degree_data)){
#		if(degree_data$log_freq[i] != all_values2$log_type_freq[i]){
#			cat("Something wrong here. Type frequencies don't match. row = ", i, ".\n")
#		}
#	}

	all_values_trim = all_values2

## Should nodes with degree 0 be trimmed?
	if(length(grep("russian", file_name)) > 0){
		all_values_trim = subset(all_values2, all_values2$degree > 0)
	}

## Should any other nodes be trimmed?
	if(length(grep("icelandic", file_name)) > 0){
		all_values_trim = subset(all_values_trim, all_values_trim$log_type_freq < 4) # Trim high type frequency classes because of data sparsity
	}
	
	jpeg(paste0("graphs/",file_name,"/complexity_byDegree&SumWeight_panel_", file_name, ".jpeg"), width = 10, height = 3, units="in", res=600)
	par(bg=NA)
	par(mfrow=c(1, 3))
	par(mar=c(3.5,4.5,1,1))
	
	plot(ent_difference ~ degree, data = all_values2, xlab="", ylab="", las=1)
	temp = all_values_trim[order(all_values_trim$degree),]
#	loess_fit <- loess(ent_difference ~ degree, data = temp)
#	lines(temp$degree, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(ent_difference ~ I(degree^2) + degree, data = temp)
#	xs <- seq(25, max(temp$degree))
#	lines(temp$degree, predict(quad_fit), col="red", lwd=2)
	linear_fit <- lm(ent_difference ~ degree, data = temp)
	lines(temp$degree, predict(linear_fit), col="red", lwd=2)
	title(ylab="Entropy Difference", line=3.5, cex.lab=1)
	title(xlab="Node Degree", line=2, cex.lab=1)

	plot(ent_difference ~ resid_weight, data = all_values2, xlab="", ylab="", las=1)
	temp = all_values_trim[order(all_values_trim$resid_weight),]
#	loess_fit <- loess(ent_difference ~ resid_weight, data = temp)
#	lines(temp$resid_weight, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(ent_difference ~ I(resid_weight^2) + resid_weight, data = temp)
#	lines(temp$resid_weight, predict(quad_fit), col="red", lwd=2)
	linear_fit <- lm(ent_difference ~ resid_weight, data = temp)
	lines(temp$resid_weight, predict(linear_fit), col="red", lwd=2)
	title(xlab="Summed Edge Weight (residualized)", line=2, cex.lab=1)

	plot(ent_difference ~ clusterCoef_local, data = all_values2, xlab="", ylab="", las=1)
	temp = all_values_trim[order(all_values_trim$clusterCoef_local),]
#	loess_fit <- loess(ent_difference ~ clusterCoef_local, data = temp)
#	lines(temp$clusterCoef_local, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(ent_difference ~ I(clusterCoef_local^2) + clusterCoef_local, data = temp)
#	lines(temp$clusterCoef_local, predict(quad_fit), col="red", lwd=2)
	linear_fit <- lm(ent_difference ~ clusterCoef_local, data = temp)
	lines(temp$clusterCoef_local, predict(linear_fit), col="red", lwd=2)
	title(xlab="Local Cluster Coefficient", line=2, cex.lab=1)
	
	dev.off()

###

# Make effect plot for entropy difference by degree, local clustering coefficient, and residualized edge weight
	
## WARNING!!!: Starting here, models need to be fit by hand to each data set. If no pre-existing model is found for a given data set, then the
## full model is used as a default. Script will not fail. But use of the full model is unlikely to be appropriate. Warning message is printed.
## CUSTOMIZE SCRIPT HERE FOR EACH DATA SET.

#Full model -- testing interactions
		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*log_type_freq + clusterCoef_local*log_type_freq + resid_weight*log_type_freq + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight + degree*log_type_freq*clusterCoef_local, data = all_values_trim)
	
#Best fit for Russian
	if(length(grep("russian", file)) > 0){
#Untrimmed
#		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)

#Trimmed
		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("chinantec", file)) > 0){
#Best fit for Palantla Chinantec
#		model_diff = lm(ent_difference ~ log_type_freq + degree + resid_weight, data = all_values_trim)
#For graphing:
	model_diff = lm(ent_difference ~ log_type_freq + degree + resid_weight + clusterCoef_local, data = all_values_trim)
	}else if(length(grep("french", file)) > 0){	
#Best fit for French
#		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + degree*log_type_freq + clusterCoef_local*log_type_freq + clusterCoef_local*degree, data = all_values_trim)
#For graphing:
		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*log_type_freq + clusterCoef_local*log_type_freq + clusterCoef_local*degree, data = all_values_trim)
	}else if(length(grep("greek", file)) > 0){
#Best fit for Greek
#		model_diff = lm(ent_difference ~ degree + clusterCoef_local + resid_weight + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
#For graphing:
		model_diff = lm(ent_difference ~ degree + clusterCoef_local + resid_weight + log_type_freq + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("icelandic", file)) > 0){
#Best fit for Icelandic
		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("kadiweu", file)) > 0){
#Best fit for Kadiweu
#		model_diff = lm(ent_difference ~ log_type_freq + clusterCoef_local + resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
#For graphing:
		model_diff = lm(ent_difference ~ log_type_freq + clusterCoef_local + resid_weight + degree + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("nuer", file)) > 0){
#Best fit for Nuer -- CHECK with summed weight -- not residualized on clusterCoef?
		model_diff = lm(ent_difference ~ log_type_freq + resid_weight, data = all_values_trim)
	}else if(length(grep("seri", file)) > 0){
#Best fit for Seri
		model_diff = lm(ent_difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*resid_weight, data = all_values_trim)
	}else if(length(grep("voro", file)) > 0){
#Best fit for Voro
#		model_diff = lm(ent_difference ~  degree + resid_weight, data = all_values_trim)
#For graphing
		model_diff = lm(ent_difference ~  degree + resid_weight + clusterCoef_local + log_type_freq, data = all_values_trim)
	}else{
		cat("Could not identify fitted model for ", file_name, ". Used full regression model. Warning: This may not be appropriate.\n")
	}
####

#Graphing
#For some reason these plots ignore cex and las. Not sure whether that has to do with being an effect plot or being multipanel
	ymax = max(all_values_trim$ent_difference) + 0.0001
	ymin = min(all_values_trim$ent_difference) - 0.0001
	plot3 = plot(effect("resid_weight", model_diff), ylab="Entropy Difference", xlab="Sum Edge Weight (resid.)\n(= # of Cell Overlaps)", main=NULL, xlim=c(-70,45), ylim=c(ymin, ymax))
	plot4 = plot(effect("log_type_freq", model_diff), ylab="Entropy Difference", xlab="Log Node Size\n(= Log Type Frequency)", main=NULL, ylim=c(ymin, ymax))
	
#For Nuer only (because it has fewer main factors):
	if(length(grep("nuer", file)) > 0){
		jpeg(paste0("graphs/", file_name,"/effect_plot_main_", file_name, ".jpg"), width=6, height=3, units="in", res=600)
		grid.arrange(plot3, plot4, ncol=2)
		dev.off()
	}else{
		plot1 =	plot(effect("degree", model_diff), ylab="Entropy Difference", xlab="Degree\n(= # of Class Overlaps)", main=NULL, ylim=c(ymin, ymax), cex.lab=2)
		plot2 = plot(effect("clusterCoef_local", model_diff), ylab="Entropy Difference", xlab="Local Clustering\nCoefficient", main=NULL, xlim=c(0, 1), ylim=c(ymin, ymax))
		jpeg(paste0("graphs/",file_name,"/effect_plot_main_", file_name, ".jpg"), width=8, height=7, units="in", res=600)
		grid.arrange(plot4, plot1, plot3, plot2, ncol=2)
		dev.off()
	}

#This is really ugly, but it seems better not to produce meaningless graphs (which could cause confusion). So the goList is hand-generated to match interactions in the model above
	goList1 <- c("french")
	goList2 <- c("french")
	goList3 <- c("french", "greek")
	goList4 <- c("russian", "greek", "seri")
	goList5 <- c("russian", "greek", "icelandic", "kadiweu")

	if(length(grep(language, goList1)) > 0){
		jpeg(paste0("graphs/",language,"/effect_plot1_", file_name, ".jpg"), width=7.5, height=4, units="in", res=600)
		plot5 = plot(effect("log_type_freq:degree", model_diff), ylab="Contribution to\nSystem Complexity (bits)", xlab="Log Type Frequency", main=NULL)
		grid.arrange(plot5, ncol=1)
		dev.off()
	}
	
	if(length(grep(language, goList2)) > 0){
		jpeg(paste0("graphs/",language,"/effect_plot2_", file_name, ".jpg"), width=7.5, height=4, units="in", res=600)
		plot6 = plot(effect("log_type_freq:clusterCoef_local", model_diff), ylab="Contribution to\nSystem Complexity (bits)", xlab="Local Clustering Coefficient", main=NULL)
		grid.arrange(plot6, ncol=1)
		dev.off()
	}
	
	if(length(grep(language, goList3)) > 0){
		jpeg(paste0("graphs/",language,"/effect_plot3_", file_name, ".jpg"), width=7.5, height=4, units="in", res=600)
		plot7 = plot(effect("degree:clusterCoef_local", model_diff), ylab="Contribution to\nSystem Complexity (bits)", xlab="Local Clustering Coefficient", main=NULL)
		grid.arrange(plot7, ncol=1)
		dev.off()
	}
	
	if(length(grep(language, goList4)) > 0){
		jpeg(paste0("graphs/",language,"/effect_plot4_", file_name, ".jpg"), width=7.5, height=4, units="in", res=600)
		plot8 = plot(effect("degree:resid_weight", model_diff), ylab="Contribution to\nSystem Complexity (bits)", xlab="Degree", main=NULL)
		grid.arrange(plot8, ncol=1)
		dev.off()
	}
	
	if(length(grep(language, goList5)) > 0){
		jpeg(paste0("graphs/",language,"/effect_plot5_", file_name, ".jpg"), width=7.5, height=4, units="in", res=600)
		plot9 = plot(effect("clusterCoef_local:resid_weight", model_diff), ylab="Contribution to\nSystem Complexity (bits)", xlab="Local Clustering Coefficient", main=NULL)
		grid.arrange(plot9, ncol=1)
		dev.off()
	}

	return(model_diff)
	
}

## END CUSTOMIZATION OF SCRIPT
