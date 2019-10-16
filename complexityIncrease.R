##########
#
# Calculate how much any individual inflection class uniquely increases (or decreases) the 
# complexity of an inflection class system
#
# Andrea D. Sims
# v. 6-12-19
#
##########

library(igraph)
library(ggplot2)
library(ggrepel)
library(effects)
library(gridExtra)

##########

#input = filename for raw file of inputs. 
#full = filename for pre-calculated entropy numbers for entire inflectional system -- no dropped classes.

calculate.fnc = function(input, full){
	complexity_data <- complexityIncrease.fnc(input, full) #Run this only if doing for the first time.
	complexityIncreaseGraph.fnc(complexity_data)
}

###########

#See how much each class contributes to the complexity of the system by comparing system complexity of all classes N to
#system complexity of N-1 classes. Interatively over all classes.

#Input is filename for raw matrix of formatives + filename for pre-calculated entropy numbers for entire inflectional system -- no dropped classes

complexityIncrease.fnc = function(input, full){
	source("script_entropy_calc.R")
	source("script_get_values.R")
	weighted = F

	file = read.delim(input, header=T, sep="\t")
	file[,1] = paste("Class",c(1:nrow(file)), sep="") #Change class id's because they are a pain
#Strip extension
	file_name = unlist(strsplit(input, ".txt"))
	for(j in 1:nrow(file)){
#Remove one row (i.e. one inflection class), and calculate entropy over the remainder
		type_freq_dropped = file[j,2] #Type frequencies are in column 2 in input file.
		data = file[-c(j),]
#Output directory
		out_file = paste("drop", j, "_", file_name, sep="")
		out = paste("subtract_classes/", out_file, sep="")
#Write out individual inflection class data -- one class removed
		write.table(data, out, quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
#Do entropy calculations -- calling script_entropy_calc.R
		doCalc.fnc(out, weighted)
#File that is the output of doCalc.fnc()
		get_ent = paste("2_pairwise_numbers/", strsplit(out_file,".txt"), "_entropy_nonweighted_Rcalc.txt", sep="")
#Get average numbers (inflectional system entropy) for output of doCalc.fnc() (one dropped class)
		values = get_cond_ent.fnc(get_ent)
		values = cbind(type_freq_dropped, values)
#Write out results to a composite file
		if(j == 1){
			write.table(values, paste("subtract_classes/", file_name, "_dropSummary.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T, fileEncoding = "UTF-8")
		} else {
			write.table(values, paste("subtract_classes/", file_name, "_dropSummary.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=F, fileEncoding = "UTF-8", append=T) #append and don't write out column names
		}
	}
}

getComplexityIncrease.fnc = function(input, full){
#Calculate complexity increase numbers after generation of drop class data (complexityIncrease.fnc())

	file_name = unlist(strsplit(input, ".txt"))
	file_name = unlist(strsplit(file_name, ".*/"))
	if(length(file_name) > 1){
		file_name = file_name[2]
	}
#Read in composite/summary data					
	all_values = read.delim(paste("subtract_classes/", file_name, "_dropSummary.txt", sep=""), header=T, sep="\t")
	all_values = as.data.frame(all_values)
#Read in pre-calculated numbers for entire inflectional system -- no dropped classes
	full_system = read.delim(full, header=T, sep="\t")		
##BIG FLAW HERE
	full_system = subset(full_system, full_system$weighted == 0)
	if(nrow(full_system) != 1){
		cat("Not extracting actual entropy value correctly.")
	}
	full_ent = full_system$actual_ent
#Calculate difference between entropy of full system and entropy when one class is dropped. Keep track of type frequency of
#the dropped class.
	difference = full_ent - as.numeric(as.character(all_values$actual_ent))
	log_type_freq = log(as.numeric(as.character(all_values$type_freq_dropped)))
	all_values = cbind(all_values, difference, log_type_freq)
	all_values = as.data.frame(all_values)
#Write out results.
	write.table(all_values, paste("entropy/entropy_diff_calcs_", file_name,".txt", sep=""), col.names=T, row.names=F, quote=F, sep="\t")
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

		all_values = read.delim(paste("entropy_diff_calcs_", file_name, ".txt", sep=""), header=T, sep="\t")
		all_values = all_values[order(all_values$log_type_freq),]
		
		all_values = as.data.frame(all_values)
		
		if(length(grep("icelandic", file_name)) > 0){
			trimmed = subset(all_values, all_values$log_type_freq < 4) # Trim high type frequency classes because of data sparsity
		} else if(length(grep("french", file_name)) > 0){
			trimmed = subset(all_values, all_values$log_type_freq < 6) # Trim high type frequency classes because of data sparsity
		} else{
			trimmed = all_values
		}
		linear_fit = lm(trimmed$difference ~ trimmed$log_type_freq)
		pVal = coef(summary(linear_fit))["trimmed$log_type_freq", "Pr(>|t|)"]
		rSquare = summary(linear_fit)$r.squared

#The following for individual graphs of each language

		language <- unlist(strsplit(file_name, "_.*"))

		graph_name = paste0("graphs/",language,"/dropClass_", file_name, ".png")
		png(graph_name, width=680, height=680, units="px")
		par(mar = c(6,6,2,2 + 0.1))
			
		plot(all_values$difference ~ all_values$log_type_freq, xlab="Log Type Frequency of Dropped Class", ylab="Entropy Difference", ylim=c(-0.02, 0.02), cex=1.25, cex.axis=1.5, cex.lab=2, col="black", pch=21, bg="gray85")
		lines(trimmed$log_type_freq, predict(linear_fit), col="black", lwd=4)
		legend("bottomright", c(paste("p = ", round(pVal, digits = 5), sep=""), paste("R^2 = ", round(rSquare, digits = 3), sep="")), cex=1.5)
		
		dev.off()
	
#		cat("\nDone with", input[k], " entropy diff x type frequency \n")
		
	}
}

###################################

complexityComponentsGraph.fnc = function(file, degree_data){

#Read in complexity increase data
	file_name = unlist(strsplit(file, ".txt"))
	file_name = unlist(strsplit(file_name, ".*/"))
#Also identify just language name (may or may not be the same as file_name) in order to define output file later, and which regressions to run
#E.g. delete anything connected by underscore, following naming convention: e.g. russian_classes_79 --> russian
	language = unlist(strsplit(file_name, "_.*"))

	if(length(file_name) > 1){
		file_name = file_name[2]
	}
	all_values2 = read.delim(paste("entropy/entropy_diff_calcs_", file_name, ".txt", sep=""), header=T, sep="\t")
	all_values2 = all_values2[order(all_values2$filename),]
	degree_data = degree_data[order(as.character(degree_data$id)),]
	all_values2 = cbind(all_values2, degree_data)
	all_values2 = as.data.frame(all_values2)

#Check that the data sets merged correctly.
	for(i in 1:length(all_values2$type_freq_dropped)){
		if(all_values2$type_freq_dropped[i] != all_values2$typeFreq[i]){
			cat("Data sets did not merge correctly!\n")
		}
	}
	
#	write.table(all_values2, "temp.txt", col.names=T, row.names=F, quote=F, sep="\t")

#Check that the data are paired correctly. Seems to be true, but this check doesn't work for some reason -- throws false warnings..
#	for(i in 1:nrow(degree_data)){
#		if(degree_data$log_freq[i] != all_values2$log_type_freq[i]){
#			cat("Something wrong here. Type frequencies don't match. row = ", i, ".\n")
#		}
#	}
#	write.table(degree_data, file=paste("degree_data_temp.txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")

	all_values_trim = all_values2

## Should nodes with degree 0 be trimmed?
	if(length(grep("russian", file_name)) > 0){
		all_values_trim = subset(all_values2, all_values2$degree > 0)
	}

## Should any other nodes be trimmed?
	if(length(grep("icelandic", file_name)) > 0){
		all_values_trim = subset(all_values_trim, all_values_trim$log_type_freq < 4) # Trim high type frequency classes because of data sparsity
	}
	
	jpeg(paste0("graphs/",language,"/complexity_byDegree&SumWeight_panel_", file_name, ".jpeg"), width = 10, height = 3, units="in", res=600)
	par(bg=NA)
	par(mfrow=c(1, 3))
	par(mar=c(3.5,4.5,1,1))
	
	plot(difference ~ degree, data = all_values2, xlab="", ylab="", las=1)
	temp = all_values_trim[order(all_values_trim$degree),]
#	loess_fit <- loess(difference ~ degree, data = temp)
#	lines(temp$degree, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(difference ~ I(degree^2) + degree, data = temp)
#	xs <- seq(25, max(temp$degree))
#	lines(temp$degree, predict(quad_fit), col="red", lwd=2)
	linear_fit <- lm(difference ~ degree, data = temp)
	lines(temp$degree, predict(linear_fit), col="red", lwd=2)
	title(ylab="Entropy Difference", line=3.5, cex.lab=1)
	title(xlab="Node Degree", line=2, cex.lab=1)
	
#	all_values4 = all_values2[order(all_values2$summedWeight),]
#	plot(all_values4$difference ~ all_values4$summedWeight, xlab="", ylab="", las=1)
#	temp = all_values4
#	loess_fit <- loess(difference ~ summedWeight, data = temp)
#	lines(temp$summedWeight, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(difference ~ I(summedWeight^2) + summedWeight, data = temp)
#	lines(temp$summedWeight, predict(quad_fit), col="red", lwd=2)
#	linear_fit <- lm(difference ~ summedWeight, data = temp)
#	lines(temp$summedWeight, predict(linear_fit), col="green", lwd=2)
#	title(xlab="Sum Weight of Node's Edges", line=2, cex.lab=1)

	plot(difference ~ resid_weight, data = all_values2, xlab="", ylab="", las=1)
	temp = all_values_trim[order(all_values_trim$resid_weight),]
#	loess_fit <- loess(difference ~ resid_weight, data = temp)
#	lines(temp$resid_weight, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(difference ~ I(resid_weight^2) + resid_weight, data = temp)
#	lines(temp$resid_weight, predict(quad_fit), col="red", lwd=2)
	linear_fit <- lm(difference ~ resid_weight, data = temp)
	lines(temp$resid_weight, predict(linear_fit), col="red", lwd=2)
	title(xlab="Summed Edge Weight (residualized)", line=2, cex.lab=1)

	plot(difference ~ clusterCoef_local, data = all_values2, xlab="", ylab="", las=1)
	temp = all_values_trim[order(all_values_trim$clusterCoef_local),]
#	loess_fit <- loess(difference ~ clusterCoef_local, data = temp)
#	lines(temp$clusterCoef_local, predict(loess_fit), col="blue", lwd=2)
#	quad_fit <- lm(difference ~ I(clusterCoef_local^2) + clusterCoef_local, data = temp)
#	lines(temp$clusterCoef_local, predict(quad_fit), col="red", lwd=2)
	linear_fit <- lm(difference ~ clusterCoef_local, data = temp)
	lines(temp$clusterCoef_local, predict(linear_fit), col="red", lwd=2)
	title(xlab="Local Cluster Coefficient", line=2, cex.lab=1)
	
	dev.off()

###

#This one didn't turn out as interesting (harder to interpret), compared to residualized summed weight

#	tiff("complexity_byDegree&MeanWeight_panel.tiff", width = 7, height = 3, units="in", res=600)
#	par(bg=NA)
#	par(mfrow=c(1, 2))
#	par(mar=c(3.5,4.5,1,1))
#	
#	plot(all_values3$difference ~ all_values3$degree, xlab="", ylab="", las=1)
#	temp = subset(all_values3, all_values3$degree >=25)
#	loess_fit <- loess(difference ~ degree, data = temp)
#	lines(temp$degree, predict(loess_fit), col="red", lwd=2)
#	title(ylab="Entropy Difference", line=3.5, cex.lab=1)
#	title(xlab="Node Degree", line=2, cex.lab=1)
#	
#	all_values5 = all_values2[order(all_values2$meanWeight),]
#	plot(all_values5$difference ~ all_values5$meanWeight, xlab="", ylab="", las=1)
#	temp = subset(all_values5, all_values5$meanWeight >= 2.5)
#	loess_fit <- loess(difference ~ meanWeight, data = temp)
#	lines(temp$meanWeight, predict(loess_fit), col="red", lwd=2)
#	title(xlab="Mean Weight of Node's Edges", line=2, cex.lab=1)
#	dev.off()

###

# Make effect plot for entropy difference by degree, local clustering coefficient, and residualized edge weight
	

#Full model -- testing interactions
#		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*log_type_freq + clusterCoef_local*log_type_freq + resid_weight*log_type_freq + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight + degree*log_type_freq*clusterCoef_local, data = all_values_trim)
	
#Best fit for Russian
	if(length(grep("russian", file)) > 0){
#Untrimmed
#		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)

#Trimmed
		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("chinantec", file)) > 0){
#Best fit for Palantla Chinantec
#		model_diff = lm(difference ~ log_type_freq + degree + resid_weight, data = all_values_trim)
#For graphing:
	model_diff = lm(difference ~ log_type_freq + degree + resid_weight + clusterCoef_local, data = all_values_trim)
	}else if(length(grep("french", file)) > 0){	
#Best fit for French
#		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + degree*log_type_freq + clusterCoef_local*log_type_freq + clusterCoef_local*degree, data = all_values_trim)
#For graphing:
		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*log_type_freq + clusterCoef_local*log_type_freq + clusterCoef_local*degree, data = all_values_trim)
	}else if(length(grep("greek", file)) > 0){
#Best fit for Greek
#		model_diff = lm(difference ~ degree + clusterCoef_local + resid_weight + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
#For graphing:
		model_diff = lm(difference ~ degree + clusterCoef_local + resid_weight + log_type_freq + clusterCoef_local*degree + degree*resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("icelandic", file)) > 0){
#Best fit for Icelandic
		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("kadiweu", file)) > 0){
#Best fit for Kadiweu
#		model_diff = lm(difference ~ log_type_freq + clusterCoef_local + resid_weight + clusterCoef_local*resid_weight, data = all_values_trim)
#For graphing:
		model_diff = lm(difference ~ log_type_freq + clusterCoef_local + resid_weight + degree + clusterCoef_local*resid_weight, data = all_values_trim)
	}else if(length(grep("nuer", file)) > 0){
#Best fit for Nuer -- CHECK with summed weight -- not residualized on clusterCoef?
		model_diff = lm(difference ~ log_type_freq + resid_weight, data = all_values_trim)
	}else if(length(grep("seri", file)) > 0){
#Best fit for Seri
		model_diff = lm(difference ~ log_type_freq + degree + clusterCoef_local + resid_weight + degree*resid_weight, data = all_values_trim)
	}else if(length(grep("voro", file)) > 0){
#Best fit for Voro
#		model_diff = lm(difference ~  degree + resid_weight, data = all_values_trim)
#For graphing
		model_diff = lm(difference ~  degree + resid_weight + clusterCoef_local + log_type_freq, data = all_values_trim)
	}else{
		cat("Could not identify language to produce regression model.\n")
	}
####

#Graphing
#For some reason these plots ignore cex and las. Not sure whether that has to do with being an effect plot or being multipanel
	ymax = max(all_values_trim$difference) + 0.0001
	ymin = min(all_values_trim$difference) - 0.0001
	plot3 = plot(effect("resid_weight", model_diff), ylab="Entropy Difference", xlab="Sum Edge Weight (resid.)\n(= # of Cell Overlaps)", main=NULL, xlim=c(-70,45), ylim=c(ymin, ymax))
	plot4 = plot(effect("log_type_freq", model_diff), ylab="Entropy Difference", xlab="Log Node Size\n(= Log Type Frequency)", main=NULL, ylim=c(ymin, ymax))
	
#For Nuer only (because it has fewer main factors):
	if(length(grep("nuer", file)) > 0){
		jpeg(paste0("graphs/",language,"/effect_plot_main_", file_name, ".jpg"), width=6, height=3, units="in", res=600)
		grid.arrange(plot3, plot4, ncol=2)
		dev.off()
	}else{
		plot1 =	plot(effect("degree", model_diff), ylab="Entropy Difference", xlab="Degree\n(= # of Class Overlaps)", main=NULL, ylim=c(ymin, ymax), cex.lab=2)
		plot2 = plot(effect("clusterCoef_local", model_diff), ylab="Entropy Difference", xlab="Local Clustering\nCoefficient", main=NULL, xlim=c(0, 1), ylim=c(ymin, ymax))
		jpeg(paste0("graphs/",language,"/effect_plot_main_", file_name, ".jpg"), width=8, height=7, units="in", res=600)
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

