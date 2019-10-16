##################################
#
# calculate and graph measures for inflection class typology
#
# Andrea D. Sims
# v. 8-5-19
# 
##################################

## Create input data structure for toy IC system (not called).

toy.fnc = function(){
	X = c("a", "b", "c")
	Y = c("d", "d", "e")
	Z = c("f", "f", "f")
	id = c("I", "II", "III")
	typeFreq = c(1,1,1)

	toy_data <- cbind(id, typeFreq, X, Y, Z)
	toy_data = as.data.frame(toy_data)

## Plot undirected network of inflection class overlap for toy IC system

	jpeg("toy_graph.jpg", width=2, height=2, units="in")
	par(mar=c(1,1,1,1))
	par(bg=NA)
	adjMatrix2 <- classMatrix.fnc(toy_data)
	plot <- graph_from_adjacency_matrix(adjMatrix2, mode="upper", weight=T)
	layout = layout.fruchterman.reingold(plot)
	plot(plot, layout=layout, edge.width=(E(plot)$weight)^1.5, vertex.color = "white", vertex.label.color = "black", edge.color = "black", vertex.size = 40, edge.label=c("d,f", "f", "f"), edge.label.x=c(-0.1,-0.45,0.8), edge.label.y = c(0.95, -0.35, -0.05), edge.label.color = "black")
	dev.off()
}

##################################

## Calculate graphing parameters for an undirected graph of class overlaps (node = inflection class) using igraph package
## Called by graphing function 

parameters.fnc = function(data3, adjMatrix){

	size = log(as.numeric(as.character(data3$typeFreq))+1) * 1.5	#Plot size of nodes -- weighted by log type frequency of class

# Calculate graph containing all nodes/edges
	plot5 <- graph_from_adjacency_matrix(adjMatrix, mode="upper", weight=T)

# All edges to be plotted
	half_cells = round((ncol(data3)-2)/2)	# Calculate how many cells = half of them. (-2 in equation because first two columns of file aren't cells)

# Color of edges
	color = vector(length = length(E(plot5)$weight))
	for(i in 1:length(color)){
		if(E(plot5)$weight[i] >= 0){
			color[i] = "green"
		}
		if(E(plot5)$weight[i] == (half_cells)){
			color[i] = "grey"
		} 
		if(E(plot5)$weight[i] >= (half_cells+1)){
			color[i] = "black"
		}
	}
	par = vector("list", 2)
	par[[1]] = size
	par[[2]] = color
	return(par)
}

###############

## Plot network for inflection class system according to various parameter settings

# modules (numeric) = vector containing group that each node belongs to, to be used for coloring. Grouping can be either by module or by betweenness centrality
#If moduleType (character) = "modules", modules can be any number of modules
#If moduleType = "betweenness", expects number of modules to be 9
#If moduleType = anything else, all nodes are colored red

#colorByModule (boolean) = should nodes be colored according to their groups (modules)?

#trim (boolean) = should weak edges in the graph be trimmed, or not?

plotNetwork.fnc = function(language, plat, adjMatrix, modules, colorByModule, moduleType, trim){

#Calculate node size and edge colors
	par_out <- parameters.fnc(plat, adjMatrix)
	size = par_out[[1]]
	color = par_out[[2]]
		
	plot <- graph_from_adjacency_matrix(adjMatrix, mode="upper", weight=T)
##Trim all but the strongest edges -- used to redraw strongest edges as last step, for better visualization. Also for layout
	plot.copy.2 <- delete.edges(plot, which(color != "black"))

	layout = layout.fruchterman.reingold(plot.copy.2)
	
	numModules <- length(unique(modules))

#Calculate color of nodes by groups
	if(colorByModule == T){
		if(moduleType == "modules"){
			
			if(numModules >= 3){
				nodeColors <- brewer.pal(numModules, "Spectral")
			}else{
				nodeColors <- brewer.pal(numModules, "Set3")
			}
			for(r in 1:length(modules)){
				temp <- modules[r]
				V(plot)$nodeColor[r] <- nodeColors[temp+1] #+1 b/c modules start at value 0, but indexed positions start at 1
			}

#Color nodes from light to dark red, based on betweenness centrality (dark = higher)
		}else if(moduleType == "betweenness"){
			nodeColors <- brewer.pal(numModules, "Reds")
			for(r in 1:length(modules)){
				temp <- modules[r]
				V(plot)$nodeColor[r] <- nodeColors[temp+1] #+1 b/c modules start at value 0, but indexed positions start at 1
			}
		}
	}else{
#All nodes red
		V(plot)$nodeColor <- "red" 
	}

#Name output files
	if(trim == T){
		file_out <- paste0("graphs/", language, "/", language, "_network_trimmed_", moduleType,".jpg")
	}else{
		file_out <- paste0("graphs/", language, "/", language, "_network_full_", moduleType,".jpg")
	}
#General graph parameters
	jpeg(file_out, width=4, height=4, units="in", res=600)
	par(mar=c(1,1,1,1))
	par(bg=NA)			

#Graph depending on whether weak nodes should be trimmed, and 	
	if(trim==T){
#Plot graph with weak edges trimmed
		color = par_out[[2]] #Edge color calculated by parameters.fnc, for trimmed graph only

##Trim weak edges
		plot.copy <- delete.edges(plot, which(color == "green"))
		color2 = color[which(color != "green")]
##Trim all but the strongest edges -- used to redraw strongest edges as last step, for better visualization
		plot.copy.2 <- delete.edges(plot, which(color != "black"))
		color3 = color[which(color == "black")]

		plot(plot, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.width = 0.05, edge.color = "gray80")	
		plot(plot.copy, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.color = color2, add=T)
		plot(plot.copy.2, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.color = color3, add=T)
	}else{
#Plot graph with all edges, with edge color and thickness according to weight	
		maxWeight<- as.numeric(as.character(max(E(plot)$weight)))
		palette <- colorRampPalette(c("gray99","gray0"))
		edgeColors <- palette(maxWeight)
	
#Edge colors
		E(plot)$edgeColor <- edgeColors[E(plot)$weight]

#Rescaling to facilitate line width by edge wedge
		E(plot)$rescale = rescale(E(plot)$weight, to = c(0,1))

		half_cells = round(maxWeight/2)
		half_plus_one_cells = round(maxWeight/2)+1
		quarter_cells = round(maxWeight) - round(maxWeight/4)
##Progressively trim weak edges -- used to redraw strongest edges as last step, for better visualization
		plot.2 <- delete.edges(plot, which(E(plot)$weight < half_cells))
		E(plot.2)$edgeColor = E(plot)$edgeColor[which(E(plot)$weight >= half_cells)]
		plot.3 <- delete.edges(plot, which(E(plot)$weight < half_plus_one_cells))
		E(plot.3)$edgeColor = E(plot)$edgeColor[which(E(plot)$weight >= half_plus_one_cells)]
		plot.4 <- delete.edges(plot, which(E(plot)$weight < quarter_cells))
		E(plot.4)$edgeColor = E(plot)$edgeColor[which(E(plot)$weight >= quarter_cells)]
		plot.5 <- delete.edges(plot, which(E(plot)$weight < maxWeight))
		E(plot.5)$edgeColor = E(plot)$edgeColor[which(E(plot)$weight >= maxWeight)]


		plot(plot, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.width = (E(plot)$rescale^1.5)*2, edge.color = E(plot)$edgeColor)	
##Replot only black lines (and nodes) so that they are on top of grey ones
		plot(plot.2, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.width = (E(plot.2)$rescale^1.5)*2, edge.color = E(plot.2)$edgeColor, add=T)
		plot(plot.3, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.width = (E(plot.3)$rescale^1.5)*2, edge.color = E(plot.3)$edgeColor, add=T)
		plot(plot.4, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.width = (E(plot.4)$rescale^1.5)*2, edge.color = E(plot.4)$edgeColor, add=T)
		plot(plot.5, layout=layout, vertex.size = size, vertex.color = V(plot)$nodeColor, vertex.label=NA, edge.width = (E(plot.5)$rescale^1.5)*2, edge.color = E(plot.5)$edgeColor, add=T)
	}
	dev.off()
}


##############

##Calculate and graph degree distribution of networks

degreeDistrib.fnc = function(language, plat, adjMatrix, networkProperties){

##Full version -- all edges included
	
#How many nodes have a given degree k?
	countNodes <- table(networkProperties$degree)
#In the input file, each node is a row, so nrow(networkProperties) = number of nodes in the network
	distrib <- as.vector(countNodes)/nrow(networkProperties)
	degrees <- as.numeric(names(countNodes)) 

	ymax = max(distrib) + max(distrib)*0.05
	
	file_out <- paste0("graphs/",language,"/degreeDistrib")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	plot(distrib ~ degrees, xlab="Degree (k)", ylab="p(k)", las=1, ylim=c(0,ymax)) #Use same y-axis scale as in trimmed version, to facilitate comparison
	
	
##Fit a Laplace distribution for languages which this visually is well suited to (minimal skew but high kurtosis): Russian (trimmed + untrimmed)
	goList <- c("russian")
	if(length(grep(language, goList)) > 0){
		m = median(as.numeric(as.character(networkProperties$degree)))
		t = sd(as.numeric(as.character(networkProperties$degree)))
		networkProperties = networkProperties[order(as.numeric(as.character(networkProperties$degree))),]
		fit <- dlaplace(as.numeric(as.character(networkProperties$degree)), mean = m, sd = t)
		lines(as.numeric(as.character(networkProperties$degree)),fit, col="green", lwd=2)
	}	
## Goodness of fit for laplace fit:	
#	laplace.test(networkProperties$degree)
##	for Russian, D = 0.8438159 (= Kolmogorov-Smirnov statistic; significance at 0.05 is appx. 0.920 for a sample size of 75, and 1.074 for significance at 0.01 for same sample size, based on tables in Puig and Stephens 2000. So, reject null hypothesis of a laplace distribution)
	
#Fit a quadratic for languages which this visually seems well suited to: Seri (untrimmed), Chinantec (untrimmed)
#	goList <- c("seri", "chinantec")
#	if(length(grep(language, goList)) > 0){
##Add a red density line
#		lines(density(networkProperties$degree), col="red", lwd=2)
#Add quadratic line
#		quad_fit <- lm(distrib ~ degrees + I(degrees^2))
#		lines(degrees, predict(quad_fit), col="green", lwd=2)
#	}
#Add a vertical line at mean degree
	abline(v=mean(as.numeric(as.character(networkProperties$degree))), col="blue", lwd=2)
	dev.off()

#Histogram	
	file_out <- paste0("graphs/",language,"/degreeDistrib_hist")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	hist(as.numeric(as.character(networkProperties$degree)), breaks=16, xlab="Degree (k)", ylab="Number of Nodes", las=1, main="") #Use same y-axis scale as in trimmed version, to facilitate comparison
	dev.off()


##Trimmed version based on the half cells criterion

	half_cells = round((ncol(plat)-2)/2)	# Calculate how many cells = half of them. (-2 in equation because first two columns of file aren't cells)

	#Trim the adjacency matrix
	adjMatrix2 = adjMatrix
	for(k in 1:nrow(adjMatrix)){
		for(m in 1:ncol(adjMatrix)){
			if(adjMatrix[k,m] < half_cells){
				adjMatrix2[k,m] = 0
			}
		}
	}
	cat("degreeDistrib.fnc: Recalculating network after trimming weak edges.\n")
	networkProperties2 <- calculateDegree.fnc(adjMatrix2, plat)
	countNodes2 <- table(networkProperties2$degree)
#In the input file, each node is a row, so nrow(networkProperties) = number of nodes in the network
	distrib2 <- as.vector(countNodes2)/nrow(networkProperties2)
	degrees2 <- as.numeric(names(countNodes2)) 
	
	ymax = max(distrib2) + max(distrib2)*0.05
	
	file_out <- paste0("graphs/",language,"/degreeDistrib_trimmed")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	plot(distrib2 ~ degrees2, xlab="Degree (k)\n(weak edges trimmed)", ylab="p(k)", las=1, ylim=c(0,ymax))
#Add a red density line
#	lines(density(networkProperties2$degree), col="red", lwd=2)

##Fit a normal distribution for languages which this visually seems well suited to: Voro (trimmed), Nuer (trimmed), Kadiweu (trimmed + untrimmed), Greek (trimmed), Chinantec (trimmed)
	goList <- c("voro", "nuer", "kadiweu", "greek", "chinantec")
	if(length(grep(language, goList)) > 0){
		xfit <- seq(min(as.numeric(as.character(networkProperties2$degree))), max(as.numeric(as.character(networkProperties2$degree))))
		yfit <- dnorm(xfit, mean=mean(as.numeric(as.character(networkProperties2$degree))),sd=sd(as.numeric(as.character(networkProperties2$degree))))
		lines(xfit, yfit, col="green", lwd=2)
	}
##Fit a Laplace distribution for languages which this visually is well suited to: Russian (trimmed)
	goList <- c("russian")
	if(length(grep(language, goList)) > 0){
		m = median(as.numeric(as.character(networkProperties2$degree)))
		t = sd(as.numeric(as.character(networkProperties2$degree)))
		networkProperties2 = networkProperties2[order(as.numeric(as.character(networkProperties2$degree))),]
		fit <- dlaplace(as.numeric(as.character(networkProperties2$degree)), mean = m, sd = t)
		lines(as.numeric(as.character(networkProperties2$degree)),fit, col="green", lwd=2)
	}

## Goodness of fit for laplace fit:	
#	laplace.test(as.numeric(as.character(networkProperties2$degree)))
##	for Russian, D = 1.181342 (= Kolmogorov-Smirnov statistic; significance at 0.05 is appx. 0.920 for a sample size of 75, and 1.074 for significance at 0.01 for same sample size, based on tables in Puig and Stephens 2000. So, reject null hypothesis of a laplace distribution)

#Fit a power law distribution for languages which this visually seems well suited to: Icelandic (trimmed), French (trimmed)
	goList <- c("french", "icelandic")
	if(length(grep(language, goList)) > 0){

#exponent = shape, minimum value = location
## This draws a Pareto (power law) curve, but doesn't actual estimate alpha (shape) from the data directly. Need to work on this.
		dpareto=function(x, shape=1, location=0.1) shape * location^shape / x^(shape + 1)
		plot(function(x) {dpareto(x, shape = 0.5)}, min(degrees2), max(degrees2), col="green", lwd=2, add=T)
	}

#Add a vertical line at mean degree
	abline(v=mean(as.numeric(as.character(networkProperties2$degree))), col="blue", lwd=2)
	dev.off()

###A version that takes account of edge weighting -- degree multiplied by mean edge weight to get total number of cell overlaps with other classes
	temp <- as.numeric(as.character(networkProperties$degree))*as.numeric(as.character(networkProperties$meanWeight))
	countNodes <- table(temp)
#In the input file, each node is a row, so nrow(networkProperties) = number of nodes in the network
	distrib3 <- as.vector(countNodes)/nrow(networkProperties)
	degrees3 <- as.numeric(names(countNodes)) 

	ymax = max(distrib3) + max(distrib3)*0.05
	
	file_out <- paste0("graphs/",language,"/degreeDistrib_weighted")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	plot(distrib3 ~ degrees3, xlab="Degree (k)\n(multiplied by edge weight)", ylab="p(k)", las=1, ylim=c(0.01,ymax)) #Use same y-axis scale as in trimmed version, to facilitate comparison

#Add a vertical line at mean degree
	abline(v=mean(temp), col="blue", lwd=2)
	dev.off()
	
###A version that calculates cumulative probability (unweighted, all edges)
	countNodes <- table(as.numeric(as.character(networkProperties$degree)))
#In the input file, each node is a row, so nrow(networkProperties) = number of nodes in the network
	distrib <- as.vector(countNodes)/nrow(networkProperties)
	cumulative <- vector(length = length(distrib))
	for(m in 1:length(cumulative)){
		if(m == 1){
			cumulative[m] <- distrib[m]
		}else{
			cumulative[m] <- sum(cumulative[m-1],distrib[m])
		}
	}
	
	degrees <- as.numeric(names(countNodes)) 

	ymax = max(cumulative) + max(cumulative)*0.05
	
	file_out <- paste0("graphs/",language,"/degreeDistrib_cumulative")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	plot(cumulative ~ degrees, xlab="Degree (k)", ylab="Cumulative probability", las=1, ylim=c(0,ymax)) #Use same y-axis scale as in trimmed version, to facilitate comparison
	dev.off()
	
###Another cumulative probability (*weighted*, all edges)
	countNodes <- as.vector(table(as.numeric(as.character(networkProperties$summedWeight)))) * as.numeric(labels(countNodes)[[1]])
#In the input file, each node is a row, so nrow(networkProperties) = number of nodes in the network

	distrib <- countNodes/sum(as.numeric(as.character(networkProperties$summedWeight)))
	cumulative <- vector(length = length(countNodes))
	for(m in 1:length(countNodes)){
		if(m == 1){
			cumulative[m] <- distrib[m]
		}else{
			cumulative[m] <- sum(cumulative[m-1],distrib[m])
		}
	}
	
	
	degreeWeight <- as.numeric(names(table(as.numeric(as.character(networkProperties$summedWeight))))) 

	ymax = max(cumulative) + max(cumulative)*0.05
	
	file_out <- paste0("graphs/",language,"/degreeDistrib_cumulative_weighted")
	jpeg(paste0(file_out,".jpg"), width=4.5, height=4.5, units="in", res=600)
	par(bg=NA)
	plot(cumulative ~ degreeWeight, xlab="Degree (k) x Edge Weight", ylab="Cumulative probability", las=1, ylim=c(0,ymax)) #Use same y-axis scale as in trimmed version, to facilitate comparison
	dev.off()

}

####################

# Plot scatters of degree x edge weight for a network 
# as a two-panel graph, using the half-nodes measure as the minimal requirement for defining an edge
# (same as undirected graph above)

networkPropertyScatters.fnc = function(language, plat, degree_data){
	
#	networkProperties <- read.delim(paste0(language,"_graphProperties.txt"), header=T, sep="\t")
#	plat <- read.delim(file_in, header=T, sep="\t")

	log_freq = log(as.numeric(as.character(plat$typeFreq)))
	degree_data <- as.data.frame(cbind(degree_data, log_freq))
	degree_data3 = degree_data[order(degree_data$log_freq),]

	file_out <- paste0("graphs/",language,"/degree_weight_byFreq_panel")
	jpeg(paste0(file_out,".jpg"), width=6, height=4, units="in", res=600)
	par(bg=NA)
	par(mfrow=c(1, 2))

	ymax = max(as.numeric(as.character(degree_data3$degree))) + (max(as.numeric(as.character(degree_data3$degree)))*0.05)
	plot(as.numeric(as.character(degree)) ~ as.numeric(as.character(log_freq)), degree_data3, xlab="Log Node Size", ylab="Node Degree", ylim=c(0,ymax), las=1)
	#loess_fit <- loess(degree ~ log_freq, degree_data3)
	#lines(degree_data3$log_freq, predict(loess_fit), col="blue", lwd=2)
	lm_fit <- lm(as.numeric(as.character(degree)) ~ as.numeric(as.character(log_freq)), degree_data3)
	lines(as.numeric(as.character(degree_data3$log_freq)), predict(lm_fit), col="red", lwd=2)

	plot(as.numeric(as.character(meanWeight)) ~ as.numeric(as.character(log_freq)), degree_data3, xlab="Log Node Size", ylab="Mean Edge Weight", las=1)
	#loess_fit <- loess(meanWeight ~ log_freq, degree_data3)
	#lines(degree_data3$log_freq, predict(loess_fit), col="blue", lwd=2)
	lm_fit <- lm(as.numeric(as.character(meanWeight)) ~ as.numeric(as.character(log_freq)), degree_data3)
	lines(as.numeric(as.character(degree_data3$log_freq)), predict(lm_fit), col="red", lwd=2)

	dev.off()

## Plot degree against weight -- for Russian this produced a quadratic fit.
	degree_data4 = degree_data[order(as.numeric(as.character(degree_data$degree))),]
	ymax = max(as.numeric(as.character(degree_data$meanWeight))) + max(as.numeric(as.character(degree_data$meanWeight)))*0.05
	
	file_out <- paste0("graphs/",language,"/degree_by_weight")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	plot(as.numeric(as.character(meanWeight)) ~ as.numeric(as.character(degree)), degree_data4, ylim=c(0,ymax), xlab="Node Degree", ylab="Mean Edge Weight", las=1)

##The following only for Russian. Raises questions about whether other languages should be trimmed...

	goList <- c("russian")
	if(length(grep(language, goList)) > 0){
		degree_data4 <- degree_data4[which(as.numeric(as.character(degree_data4$degree)) > 1),]

	}
	goList <- c("russian", "chinantec")
	if(length(grep(language, goList)) > 0){
		quad_fit <- lm(as.numeric(as.character(meanWeight)) ~ as.numeric(as.character(degree)) + I(as.numeric(as.character(degree))^2), data = degree_data4)
		lines(as.numeric(as.character(degree_data4$degree)), predict(quad_fit), col="red", lwd=2)
	}
	goList <- c("french", "greek")
	if(length(grep(language, goList)) > 0){
		linear_fit <- lm(as.numeric(as.character(meanWeight)) ~ as.numeric(as.character(degree)), data = degree_data4)
		lines(as.numeric(as.character(degree_data4$degree)), predict(linear_fit), col="red", lwd=2)
	}
	
	dev.off()
}

###########

##Identify modules using simulated annealing (rnetcarto package) and calculate modularity of the network and related scores
## See Guimera et al. 2007 (Classes of complex networks deined by role-to-role connectivity profiles) for discussion, and inter alia (11, 12)

##The output column "connectivity" is used to differentiate hubs from non-hubs based on module-internal connectivity (~connectivity in my network). 
## Connectivity is a z-score that (seems to) calculates within-module connectivity, as discussed in Guimera et al. 2007 and
##earlier work. They define >=2.5 as a hub and <2.5 as a non-hub. (For Russian, everything is a non-hub.)

##The output column "participation" is used to subdivide hubs and non-hubs into subtypes based on connectivity across modules.
##Participation co-efficient as defined in Guimera et al. 2007. 
##For non-hubs:
##P <= 0.05 "ultra-peripheral"
## 0.05 < p <= 0.62 "peripheral"
## 0.62 < P <= 0.8 "satellite" (called "connector" in the output)
## 0.8 < P "kinless"

##Total modularity of the graph under the proposed partition is also calculated, as the second element in the list.

identifyModules.fnc = function(language, plat, adjMatrix){
	modules <- netcarto(adjMatrix)
	file_out <- paste0("network_properties/modules_", language, ".txt")
	write.table(modules[[1]], file_out, quote=F, sep="\t", col.names=T, row.names=F)

	color = T	#Should nodes in the network be colored according to their module?
	moduleType = "modules"  #What kind of grouping of nodes is being passed?

#Order the data frame
	modules[[1]] <- modules[[1]][order(modules[[1]]$name),]
	modules <- modules[[1]]$module
	trim=T
	plotNetwork.fnc(language, plat, adjMatrix, modules, color, moduleType, trim)
	trim=F
	plotNetwork.fnc(language, plat, adjMatrix, modules, color, moduleType, trim)
}

#############

##Plot distributions of shortest path length node-by-node, clustering coefficient, and betweenness centrality, node-by-node
##Plot network graph with nodes colored according to quartile on betweenness centrality.

plotShortestPath.fnc = function(language, networkProperties){

##Weight is usually *cost* -- e.g., length of a route. So it is normally thought of as a correlate to distance. But in my graphs it is a correlate to *closeness*. In calculating shortest path length, using weighting leads the system to look for paths from one node to another that bias towards going through minimally similar nodes, when in fact it makes more sense for the path to go through maximally similar nodes. So the "inverted" version of edge weights (edge weights reversed) is better.)

##Standard approach -- goal is to minimize edge weight (odd assumption in present context -- produces wonky results)
#	file_out <- paste0("graphs/",language,"/pathLengthDistrib")
#	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
#	par(bg=NA)
#	hist(networkProperties$meanShortestPath_local, xlab="Mean Shortest Path Length\n(Edge Weighting)", main="")
#	dev.off()

#Version without any edge weight in shortest path length calculation.	
	file_out <- paste0("graphs/",language,"/pathLengthDistrib_noWeight")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	hist(as.numeric(as.character(networkProperties$meanShortestPath_local_noWeight)), xlab="Mean Shortest Path Length\n(No Edge Weighting)", main="")
	dev.off()

#Version with edge weights reversed, to better fit how algorithm traverses the graph and calculates shortest path.	
	file_out <- 	file_out <- paste0("graphs/",language,"/pathLengthDistrib_reverseWeight")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	hist(as.numeric(as.character(networkProperties$meanShortestPath_local_reverseWeight)), xlab="Mean Shortest Path Length\n(Reverse Edge Weighting)", main="")
	dev.off()
	
#Clustering coefficient distribution
	file_out <- paste0("graphs/",language,"/clusterCoefDistrib")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	hist(as.numeric(as.character(networkProperties$clusterCoef_local)), xlab="Local Clustering Coefficient", main="")
	dev.off()
	
#Betweenness centrality distribution
	file_out <- paste0("graphs/",language,"/betweennessCentralityDistrib")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	hist(as.numeric(as.character(networkProperties$betweennessCentrality)), xlab="Betweenness Centrality", main="")
	dev.off()
	
#Scatterplot of relationship between type frequency (node size) and betweenness centrality
	file_out <- paste0("graphs/", language, "/betweennessByTypeFreq")
	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
	par(bg=NA)
	plot(as.numeric(as.character(networkProperties$betweennessCentrality)) ~ log(as.numeric(as.character(networkProperties$typeFreq))), xlab="Log Type Frequency of Class\n(= Node Size)", ylab="Node Betweenness Centrality")
	fit <- lm(as.numeric(as.character(networkProperties$betweennessCentrality)) ~ log(as.numeric(as.character(networkProperties$typeFreq))))
#Fit not significant for Russian
#	lines(log(networkProperties$typeFreq), predict(fit), col="green", lwd=2)
	dev.off()
	
#For Russian (87 classes), 8-4-19:
#Call:
#lm(formula = networkProperties$betweennessCentrality ~ log(networkProperties$typeFreq))
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-25.766 -14.088  -1.709   9.536  76.209 
#
#Coefficients:
#                                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                      25.7661     3.0884   8.343 1.15e-12 ***
#log(networkProperties$typeFreq)  -0.9204     0.7900  -1.165    0.247    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 18.69 on 85 degrees of freedom
#Multiple R-squared:  0.01572,	Adjusted R-squared:  0.004138 
#F-statistic: 1.357 on 1 and 85 DF,  p-value: 0.2473
	
#	temp <- ggplot(networkProperties, aes(x=log(typeFreq), y=betweennessCentrality)) + geom_point() + geom_smooth(method=lm, color="green") + xlab("Log Type Frequency of Class\n(= Node Size)") + ylab("Node Betweenness Centrality")
#	file_out <- paste0("graphs/", language, "/betweennessByTypeFreq")
#	jpeg(paste0(file_out,".jpg"), width=4, height=4, units="in", res=600)
#	grid.arrange(temp, ncol=1)
#	dev.off()

}

#############

#Network graph with nodes colored according to betweenness centrality, with scalar color (darker = more centrality)

plotNetworkBetweenness.fnc = function(language, plat, adjMatrix, networkProperties){
	
	breaks = max(as.numeric(as.character(networkProperties$betweennessCentrality)), na.rm=T)/9

	groups <- vector(length=nrow(networkProperties))
	for(i in 1:nrow(networkProperties)){
		if(is.na(as.numeric(as.character(networkProperties$betweennessCentrality[i])))){
			groups[i] <- 0
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks){
			groups[i] <- 0
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*2){
			groups[i] <- 1
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*3){
			groups[i] <- 2
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*4){
			groups[i] <- 3
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*5){
			groups[i] <- 4
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*6){
			groups[i] <- 5
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*7){
			groups[i] <- 6
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*8){
			groups[i] <- 7
		}else if(as.numeric(as.character(networkProperties$betweennessCentrality[i])) <= breaks*9){
			groups[i] <- 8
		}
	}
	colors <- T #Should nodes be colored according to their group?
	moduleType <- "betweenness" #What kind of grouping of nodes is being passed? (Affects color scheme used for nodes and output file name.)
	trim=T
	plotNetwork.fnc(language, plat, adjMatrix, groups, colors, moduleType, trim)
	trim=F
	plotNetwork.fnc(language, plat, adjMatrix, groups, colors, moduleType, trim)
}

#############

# Plot complexity increase against type frequency, degree, edge weight, clustering coefficient for a language.

complexityScatters.fnc = function(file, entropies_file, adjMatrix, degree_data){
	
	data <- read.delim(file, header=T, sep="\t")

#First file is file name for raw file of inputs (not dropped class -- used to identify dropped class file). Second file is name for pre-calculated entropy numbers for entire inflectional system -- no dropped classes.

#Calculate and then plot difference between complexity of full system and complexity of system with one class removed

	getComplexityIncrease.fnc(file, entropies_file)

#Plot results of calculations. 
	complexityIncreaseGraph.fnc(file)

#Calculate
#	adjMatrix2 <- classMatrix.fnc(data)
#	degree_data <- calculateDegree.fnc(adjMatrix, plat)
	id <- 1:nrow(data) #Initialize id for classes
	typeFreq <- data[,2] #Column of type frequency (node size) for inflection classes
	degree_data <- cbind(id, typeFreq, degree_data)
	filename = unlist(strsplit(file, ".txt"))
	filename = unlist(strsplit(filename, ".*/"))
	if(length(filename) > 1){
		filename = filename[2]
	}
#	write.table(degree_data, file=paste("degree_data_", filename, ".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
#Add to entropy difference calculations and plot results
	model_out = complexityComponentsGraph.fnc(file, degree_data)
	summary(model_out)	
#Checking confidence intervals on regression slopes
	confint(model_out, parm="log_type_freq")
	confint(model_out, parm="degree")
	confint(model_out, parm="clusterCoef_local")
	confint(model_out, parm="resid_weight")
	confint(model_out, parm="log_type_freq:degree")
	confint(model_out, parm="log_type_freq:clusterCoef_local")
	confint(model_out, parm="degree:clusterCoef_local")
	confint(model_out, parm="degree:resid_weight")
	confint(model_out, parm="clusterCoef_local:resid_weight")

#Checking collinearity (variance inflation)
# This doesn't work because not fixef reported... Why???!
#	source("vif-mer.R")
#	vif.mer(model_out)
	
}

########################

#Graphing t-values for main effect node-level factors (node size, degree, summed edge weight (residualized), local clustering coefficient)

graphTVals.fnc = function(){

##Where does this data set come from?
	tVals = read.delim("factor_tValues.txt", header=T, sep="\t")


##Plot raw t-values

#Set x- and y-limits
	xmax = max(tVals$nodeSize)
	xmin = min(tVals$nodeSize)
	if(xmax^2 > xmin^2){
		xmin = xmax*-1
	}else{
		xmax = xmin*-1
	}

	ymax = max(tVals$clusterCoefLocal)
	ymin = min(tVals$clusterCoefLocal)
	if(ymax^2 > ymin^2){
		ymin = ymax*-1
	}else{
		ymax = ymin*-1
	}

	plot1 <- ggplot(tVals, aes(x = nodeSize, y = clusterCoefLocal)) + geom_point(color = "gray") + geom_text_repel(aes(x = nodeSize, y = clusterCoefLocal, color = "red", label = language), size = 5) + theme_bw() + theme_minimal() + theme(legend.position="none") + xlim(xmin,xmax) + ylim(ymin,ymax) + xlab("Node Size (t-values)") + ylab("Clustering (t-values)") + geom_vline(xintercept=0) + geom_hline(yintercept=0)

#Set x- and y-limits
	xmax = max(tVals$edgeWeight)
	xmin = min(tVals$edgeWeight)
	if(xmax^2 > xmin^2){
		xmin = xmax*-1
	}else{
		xmax = xmin*-1
	}

	ymax = max(tVals$clusterCoefLocal)
	ymin = min(tVals$clusterCoefLocal)
	if(ymax^2 > ymin^2){
		ymin = ymax*-1
	}else{
		ymax = ymin*-1
	}

	plot2 <- ggplot(tVals, aes(x = edgeWeight, y = clusterCoefLocal)) + geom_point(color = "gray") + geom_text_repel(aes(x = edgeWeight, y = clusterCoefLocal, color = "red", label = language), size = 5) + theme_bw() + theme_minimal() + theme(legend.position="none") + xlim(xmin,xmax) + ylim(ymin,ymax) + xlab("Edge Weight (t-values)") + ylab("Clustering (t-values)") + geom_vline(xintercept=0) + geom_hline(yintercept=0)

	xmax = max(tVals$degree)
	xmin = min(tVals$degree)
	if(xmax^2 > xmin^2){
		xmin = xmax*-1
	}else{
		xmax = xmin*-1
	}

	ymax = max(tVals$clusterCoefLocal)
	ymin = min(tVals$clusterCoefLocal)
	if(ymax^2 > ymin^2){
		ymin = ymax*-1
	}else{
		ymax = ymin*-1
	}

	plot3 <- ggplot(tVals, aes(x = degree, y = clusterCoefLocal)) + geom_point(color = "gray") + geom_text_repel(aes(x = degree, y = clusterCoefLocal, color = "red", label = language), size = 5) + theme_bw() + theme_minimal() + theme(legend.position="none") + xlim(xmin,xmax) + ylim(ymin,ymax) + xlab("Degree (t-values)") + ylab("Clustering (t-values)") + geom_vline(xintercept=0) + geom_hline(yintercept=0)

	xmax = max(tVals$nodeSize)
	xmin = min(tVals$nodeSize)
	if(xmax^2 > xmin^2){
		xmin = xmax*-1
	}else{
		xmax = xmin*-1
	}

	ymax = max(tVals$degree)
	ymin = min(tVals$degree)
	if(ymax^2 > ymin^2){
		ymin = ymax*-1
	}else{
		ymax = ymin*-1
	}

	plot4 <- ggplot(tVals, aes(x = nodeSize, y = degree)) + geom_point(color = "gray") + geom_text_repel(aes(x = nodeSize, y = degree, color = "red", label = language), size = 5) + theme_bw() + theme_minimal() + theme(legend.position="none") + xlim(-7,7) + ylim(ymin,ymax) + xlab("Degree (t-values)") + ylab("Node Size (t-values)") + geom_vline(xintercept=0) + geom_hline(yintercept=0)

	xmax = max(tVals$degree)
	xmin = min(tVals$degree)
	if(xmax^2 > xmin^2){
		xmin = xmax*-1
	}else{
		xmax = xmin*-1
	}

	ymax = max(tVals$edgeWeight)
	ymin = min(tVals$edgeWeight)
	if(ymax^2 > ymin^2){
		ymin = ymax*-1
	}else{
		ymax = ymin*-1
	}

	plot5 <- ggplot(tVals, aes(x = degree, y = edgeWeight)) + geom_point(color = "gray") + geom_text_repel(aes(x = degree, y = edgeWeight, color = "red", label = language), size = 5) + theme_bw() + theme_minimal() + theme(legend.position="none") + xlim(xmin,xmax) + ylim(ymin,ymax) + xlab("Degree (t-values)") + ylab("Edge Weight (t-values)") + geom_vline(xintercept=0) + geom_hline(yintercept=0)


	plot1_lim <- plot1 <- ggplot(tVals, aes(x = nodeSize, y = clusterCoefLocal)) + geom_point(color = "gray") + geom_text_repel(aes(x = nodeSize, y = clusterCoefLocal, color = "red", label = language), size = 5) + theme_bw() + theme_minimal() + theme(legend.position="none") + xlim(xmin,xmax) + ylim(ymin,ymax) + xlab("Node Size (t-values)") + ylab("Clustering (t-values)") + geom_vline(xintercept=0) + geom_hline(yintercept=0)

	jpeg("graphs/tVals_1.jpg", width=4, height=4, units="in", res=600)
	par(bg=NA)
	grid.arrange(plot1, ncol=1)
	dev.off()
	
	jpeg("graphs/tVals_2.jpg", width=4, height=4, units="in", res=600)
	par(bg=NA)
	grid.arrange(plot2, ncol=1)
	dev.off()
	
	jpeg("graphs/tVals_3.jpg", width=4, height=4, units="in", res=600)
	par(bg=NA)
	grid.arrange(plot3, ncol=1)
	dev.off()
	
	jpeg("graphs/tVals_4.jpg", width=4, height=4, units="in", res=600)
	par(bg=NA)
	grid.arrange(plot4, ncol=1)
	dev.off()
	
	jpeg("graphs/tVals_5.jpg", width=4, height=4, units="in", res=600)
	par(bg=NA)
	grid.arrange(plot5, ncol=1)
	dev.off()

	jpeg("graphs/tVals_panel.jpg", width=8, height=4, units="in", res=600)
	par(bg=NA)
	grid.arrange(plot1_lim, plot5, ncol=2)
	dev.off()
}

###########

graphMonteCarlo.fnc = function(filenames){

#Scatter plot of each language's network's global clustering coefficient as a function of its mean path length
#Also a simulated version of each language that randomizes how exponents are assigned to classes

	#Extract language names from file name
	languages = filenames
	languages = gsub(".txt", "", languages)
	languages = gsub(".*/", "", languages)
	languages = gsub("_.*", "", languages)
	languages[which(languages == "palantla")] = "Chinantec"
	languages[which(languages == "nuer2")] = "Nuer"
	#Capitalize language names
	languages = paste(toupper(substr(languages,1,1)), substr(languages, 2, nchar(languages)), sep="")


#Calculate global clustering coefficient and mean shortest path length for each language
	values = data.frame(lang = vector(), meanPathLength = vector(), clusterCoef = vector())
	simValues = data.frame(lang = vector(), meanPathLength = vector(), clusterCoef = vector())
#	set.seed(500) #For rnorm() below
	for(i in 1:length(filenames)){
		temp_plat <- read.delim(filenames[i], header=T, sep="\t")
#Check for duplicate rows and fix as needed.
		temp_plat <- checkDuplicates.fnc(temp_plat)
#Check for non-unique inflection class labels and assign new ones if needed
		temp_plat <- checkLabels.fnc(temp_plat)
		
		adjMatrix2 <- classMatrix.fnc(temp_plat)
#		write.table(adjMatrix2, "adjMatrix2_temp.txt", sep="\t", quote=F)
		graph <- graph_from_adjacency_matrix(adjMatrix2, mode="upper", weight=T)
		meanPathLength = mean_distance(graph, directed=F, unconnected=T) #According to documention edge weights aren't included in calculations -- weighting not implemented for mean_distance(). 
		clusterCoef = transitivity(graph, type = "global", isolates="NaN")
		
#Version with reverse edge weight (reversed b/c weighted distance in network is correlate to closeness, but weighted path length treats it as a correlate to distance/cost and tries to minimize it.) reverseWeight.fnc() returns output of graph_from_adjacency_matrix
		graph_reverseWeight <- reverseWeight.fnc(adjMatrix2, temp_plat)
		shortestPath_weighted = distances(graph_reverseWeight, v=V(graph_reverseWeight), to = V(graph_reverseWeight)) #Edge weight used if attribute
		shortestPath_weighted <- shortestPath_weighted[which(is.finite(shortestPath_weighted))] #For disconnected graphs, where distance produces Inf
		shortestPath_weighted <- shortestPath_weighted[which(shortestPath_weighted > 0)] #Exclude node-to-same-node from mean calculation
		meanPathLength_weighted = mean(shortestPath_weighted, na.rm=T)
		clusterCoef_weighted <- transitivity(graph_reverseWeight, type = "global", isolates="NaN") #Edge weighted used if attribute

		lang = languages[i]
		
		values = rbind(values, cbind(lang, meanPathLength, clusterCoef, meanPathLength_weighted, clusterCoef_weighted))
		
		adjMatrixBinary = adjMatrix2
		for(j in 1:nrow(adjMatrixBinary)){
			for(k in 1:ncol(adjMatrixBinary)){
				if(adjMatrixBinary[j,k] > 0){
					adjMatrixBinary[j,k] <- 1
				}else{
					adjMatrixBinary[j,k] <- 0
				}
			}
		}

#Calculate values for Monte Carlo simulation of system (unweighted)		
		simNodes = nrow(adjMatrixBinary)
#Does this need to be divided in 2 because of symmetry of adjMatrix?
		simProb = sum(adjMatrixBinary[,])/(simNodes*(simNodes-1)) #Denominator is total number of possible class overlaps
		
		simGraph <- sample_gnp(simNodes, simProb, directed=F, loops=F)
		
		meanPathLength = mean_distance(simGraph, directed=F, unconnected=T)
		clusterCoef = transitivity(simGraph, type = "global", isolates="NaN")		

#Calculate values for Monte Carlo simulation of system (weighted). Using same graph but assigning edge weights to it. Actual weights from real language are kept (and thus their distrib is the same), but they are assigned at random to edges.
		meanWeight = mean(E(graph)$weight) #For existing edges in real language, what is mean edge weight

		nSamples <- length(E(simGraph)) #Sample number must be same as number of edges in simGraph
		maxWeight <- ncol(temp_plat)-3 #Minus 2 because first two columns in the plat aren't cells, and one more because a cell can't overlap with itself
		origWeight <- c(1:maxWeight)
		invertWeight <- c(maxWeight:1)
		reverseWeight <- vector(length=length(E(simGraph)))

		weightedDistrib1 <- sample(E(graph)$weight, nSamples, replace=T)
#Reverse weight
		for(y in 1:length(weightedDistrib1)){
			match <- which(origWeight == weightedDistrib1[y])
			reverseWeight[y] <- invertWeight[match]
		}
		simGraph1 <- simGraph
		E(simGraph1)$weight <- reverseWeight
		
		shortestPath_weighted = distances(simGraph1, v=V(simGraph1), to = V(simGraph1)) #Calculates using 'weight' attribute if it exists
		shortestPath_weighted <- shortestPath_weighted[which(is.finite(shortestPath_weighted))]
		shortestPath_weighted <- shortestPath_weighted[which(shortestPath_weighted > 0)] #Exclude node-to-same-node from calc
		
		meanPathLength_weighted = mean(shortestPath_weighted, na.rm=T)
		clusterCoef_weighted = transitivity(simGraph1, type = "global", isolates="NaN") #Uses 'weight' attribute if it exists
		
		lang = paste("Sim", languages[i], sep="")
			
		simValues = rbind(simValues, cbind(lang, meanPathLength, clusterCoef, meanPathLength_weighted, clusterCoef_weighted))
	}
	
	xmax = max(as.numeric(as.character(values$meanPathLength)))+0.05

	simulated = vector(length = nrow(simValues))
	simulated[] = 1
	simValues = cbind(simulated, simValues)
	simulated[] = 0
	values = cbind(simulated, values)
	data = rbind(values, simValues)
	write.table(data, "network_properties/monteCarlo_data.txt", quote=F, sep="\t", col.names=T, row.names=F)

#Graph	
	xmax = max(as.numeric(as.character(data$meanPathLength)))
	
	jpeg("graphs/pathLength_and_clusterCoef.jpg", width=6, height=6, units="in", res=600)
	par(bg=NA)

#The following produces overlapping labels
#	plot(y = as.numeric(as.character(values$clusterCoef)), x = as.numeric(as.character(values$meanPathLength)), main="", type="n", xlab="Mean Shortest Path Length", ylab="Global Clustering Coefficient", ylim=c(0,1), xlim=c(1,xmax))
#	text(x = as.numeric(as.character(simValues$meanPathLength)), y = as.numeric(as.character(simValues$clusterCoef)), simValues$lang, col="darkgray")
#	text(x = as.numeric(as.character(values$meanPathLength)), y = as.numeric(as.character(values$clusterCoef)), values$lang, col="blue")

#Fixing overlapping labels problem
	ggplot(data, aes(x = as.numeric(as.character(meanPathLength)), y = as.numeric(as.character(clusterCoef)))) + geom_point(color = "gray") + geom_text_repel(aes(x = as.numeric(as.character(meanPathLength)), y = as.numeric(as.character(clusterCoef)), group = factor(simulated), color = factor(simulated), label = lang), size = 5) + xlim(1,xmax) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("red", "black"), breaks=c("0", "1")) + theme(legend.position="none") + xlab("Mean Shortest Path Length\n(unweighted)") + ylab("Global Clustering Coefficient")
	dev.off()

#Another version with larger axes/labels, blue, lower res for digital display
	jpeg("graphs/pathLength_and_clusterCoef_blue.jpg", width=6, height=6, units="in", res=300)
	par(bg=NA)
	ggplot(data, aes(x = as.numeric(as.character(meanPathLength)), y = as.numeric(as.character(clusterCoef)))) + geom_point(aes(group = factor(simulated), color = factor(simulated)), size=3) + geom_text_repel(aes(x = as.numeric(as.character(meanPathLength)), y = as.numeric(as.character(clusterCoef)), group = factor(simulated), color = factor(simulated), label = lang), size = 6.5) + xlim(1,xmax) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("blue", "black"), breaks=c("0", "1")) + theme(legend.position="none") + xlab("Mean Shortest Path Length\n(unweighted)") + ylab("Global Clustering Coefficient") + theme(axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"))
	dev.off()
	
#Another version with larger axes/labels, red, lower res for digital display
	jpeg("graphs/pathLength_and_clusterCoef_red.jpg", width=6, height=6, units="in", res=300)
	par(bg=NA)
	ggplot(data, aes(x = as.numeric(as.character(meanPathLength)), y = as.numeric(as.character(clusterCoef)))) + geom_point(aes(group = factor(simulated), color = factor(simulated)), size=3) + geom_text_repel(aes(x = as.numeric(as.character(meanPathLength)), y = as.numeric(as.character(clusterCoef)), group = factor(simulated), color = factor(simulated), label = lang), size = 6.5) + xlim(1,xmax) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("red", "black"), breaks=c("0", "1")) + theme(legend.position="none") + xlab("Mean Shortest Path Length\n(unweighted)") + ylab("Global Clustering Coefficient") + theme(axis.text=element_text(size=14), axis.title=element_text(size=16, face="bold"))
	dev.off()
	
#Version with (reverse) edge weight used for shortest path length -- distribution of edge weights in simLanguage (probabilistically) the same as in real language
	jpeg("graphs/pathLength_and_clusterCoef_weighted.jpg", width=6, height=6, units="in", res=600)
	par(bg=NA)
	ggplot(data, aes(x = as.numeric(as.character(meanPathLength_weighted)), y = as.numeric(as.character(clusterCoef)))) + geom_point(color = "gray") + geom_text_repel(aes(x = as.numeric(as.character(meanPathLength_weighted)), y = as.numeric(as.character(clusterCoef)), group = factor(simulated), color = factor(simulated), label = lang), size = 5) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("red", "black"), breaks=c("0", "1")) + theme(legend.position="none") + xlab("Mean Shortest Path Length\n(weighted)") + ylab("Global Clustering Coefficient")
	dev.off()
	
#Version with (reverse) edge weight used for shortest path length -- mean path length in log space to space out data points
	jpeg("graphs/pathLength_and_clusterCoef_weighted_log.jpg", width=6, height=6, units="in", res=600)
	par(bg=NA)
	ggplot(data, aes(x = log(as.numeric(as.character(meanPathLength_weighted))), y = as.numeric(as.character(clusterCoef)))) + geom_point(color = "gray") + geom_text_repel(aes(x = log(as.numeric(as.character(meanPathLength_weighted))), y = as.numeric(as.character(clusterCoef)), group = factor(simulated), color = factor(simulated), label = lang), size = 5) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("red", "black"), breaks=c("0", "1")) + theme(legend.position="none") + xlab("Log of Mean Shortest Path Length\n(weighted)") + ylab("Global Clustering Coefficient")
	dev.off()

}

############


#Principle components analysis of t-values

#pca <- prcomp(tVals[,2:5], center=T, scale.=T)
#print(pca)
##                        PC1         PC2        PC3         PC4
##nodeSize          0.5871975 -0.02196286 -0.2971318 -0.75261507
##degree            0.5745102 -0.03844458 -0.5001071  0.64680209
##clusterCoefLocal -0.3104406 -0.86980144 -0.3774578 -0.06780617
##edgeWeight        0.4782926 -0.49141134  0.7205076  0.10305307
