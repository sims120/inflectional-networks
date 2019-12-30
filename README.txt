INFLECTIONAL NETWORKS: RESOURCES FOR GRAPH-THEORETIC ANALYSIS OF LINGUISTIC MORPHOLOGY

Author: Andrea D. Sims (The Ohio State University, Department of Linguistics, Department of Slavic & East European Languages & Cultures)
Contributor: Jeff Parker (Brigham Young University, Department of Linguistics & English Language)

Version: 1.0.1
Version date: 30 December 2019

####################

INTRODUCTION

This release contains R code that implements graph-theoretic and information-theoretic analyses of inflectional data. The core idea underlying the analyses is that inflection class systems can be conceptualized as networks, with the classes as nodes. Edges connect the nodes if they share exponents. Conceptualized thusly, the ways in which inflection classes are related – the extent to which they share inflectional exponents, and the distribution of exponents across classes more generally – can be quantified in terms of network properties.

First and foremost, the R code in this release implements network visualizations and a range of standard graph-theoretic calculations over inflection class networks: degree, edge weight, clustering coefficient, connected components, shortest path length, betweenness centrality, modules, node role. Mostly the calculations are for individual networks, but the code implements some minimal comparison based on global clustering coefficient and mean shortest path length, comparing the individual networks to each other, and also to randomized versions of the same class systems (based on Monte Carlo simulations).

Additionally, the code implements probability and entropy calculations for individual paradigm cells, and mean values for the inflectional system as a whole. Conditional entropy is calculated pairwise over individual cells in the style of Ackerman et al. (2009, “Parts and wholes: Implicative patterns in inflectional paradigms”) and Ackerman and Malouf (2013, “Morphological organization: The Low Conditional Entropy Conjecture”). 

Finally, it calculates "entropy difference" -- the amount of average conditional entropy in a network minus the amount of average conditional entropy with one node (class) removed -- and graphs correlations between network properties and entropy difference, to facilitate evaluation of how network properties contribute to the 'complexity' of the inflection class system.

Ultimately, the scripts in this release are designed to facilitate quantitative, cross-linguistic comparison of inflectional systems based on their systemic internal organization.

####################

GETTING STARTED

This release was built using R 3.4.1 (R Core Team 2017) in MacOS. 

It consists of the following six scripts:
main.R
adjacencyMatrix.R
entropyCalc.R
complexityIncrease.R
calcDegree.R
graphs.R

PREREQUISITES

It requires the following R packages to be installed:
igraph			built under version 1.1.1 (Csardi and Nepusz 2006)
effects			built under version 3.1-2
gridExtra		built under version 2.2.1
ggplot2			built under version 2.2.1
ggrepel			built under version 0.6.5
rnetcarto		built under version 0.2.4. For calculating network modularity (Guimerà and Amaral 2005a,b)
RColorBrewer		built under version 1.1-2
scales			built under version 0.4.1
jmuOutlier  		built under version 2.1. Not required for basic functionality -- only if fitting a laplace distribution
lawstat			built under version 3.3. Not required for basic functionality -- only if fitting a laplace distribution


It requires input data to be in the folder (relative to the local working directory that contains the scripts):

./language_plats

For format of the input data files, see below.

RUNNING

The main function calls are in main.R. When the working directory is the same as the directory containing the scripts, 

source("main.R")
main.fnc(all_files)

will run all major functions -- calculations and graphs -- for all files in language_plats/ . An unlimited number of language plats is accommodated, although the usual caveats about run time apply.

However, it is also possible to feed individual file names to main.fnc() and also to run individual subfunctionalities of the code only. Most notably:

network.fnc() calculates various graph-theoretic measures of each inflectional network, makes network and other graphs, and writes out summary files of network properties. It only calculates measures for each language individually.

comparison.fnc() graphs some limited comparisons across all plats -- e.g. plotting mean shortest path length against global clustering coefficient -- and adds Monte Carlo simulations of each language.

entropy.fnc() makes probability/entropy calculations: probability of exponents (unconditioned and conditioned on one other exponent), unconditioned entropy of a cell, entropy of a cell conditioned on one other cell (in the style of Ackerman et al. 2009 and Ackerman and Malouf 2013). It also calculates "entropy difference" -- the difference between the mean conditional entropy of paradigm cells in a full inflectional system and the mean conditional entropy of paradigm cells of the same system with one class removed. It iterates through all classes in the plat and outputs a summary file for each modified IC system.

All three of these functions take a vector of file names as an argument (it assumes the listed files are in language_plats/ -- the directory is not passed as part of the file name).

####################

INPUT FILE FORMAT

A 'language plat' is .txt file containing a tab-delimited table of inflectional exponents. Except for the header row, each row in the input file represents an inflection class. Except for the first two columns, each column represents a paradigm cell (i.e. set of inflectional values). Each [row,column] combination is filled with an inflectional exponent for that paradigm cell in that inflection class. The form should NOT be a full word-form. This script does not do segmentation into stems and affixes (or themes and distinguishers in the terminology of Stump and Finkel 2013). The assumption is that appropriate segmentation has already been done and the input language plat reflects the result.

The first column of the input file must contain a label for the inflection class. (These can be non-unique, even all "NA". Uniqueness is checked within the script and random numbers are assigned if non-unique labels are found.) The second column of the input file must contain type frequency counts for the class (i.e. the number of lexemes in that class). (If type frequency is irrelevant to the analysis or unknown, these can all be 1.) The script assumes that paradigm cells begin in column 3. Any number of columns and rows is accommodated. In general, column names are not important (as long as they exist). The script checks for uniqueness of exponents (columns 3 to ncol) and combines any rows that are found to be identical. It also prints to the console the row numbers that have been combined.

An example of the input format:

class	Freq	A	B	C	D	E	F		#header row. A-F are paradigm cells (morphosyntactic property sets)
1	279	b	m	y	vv	xx	iii		#each subsequent row is an inflection class
2	169	i	t	ff	qq	eee	ppp
3	67	d	o	aa	ll	zz	kkk
4	53	e	p	bb	mm	aaa	lll
5	41	d	p	aa	mm	zz	lll
6	38	l	x	ii	uu	hhh	ttt
7	34	i	u	ff	rr	eee	qqq
8	29	l	w	ii	tt	hhh	sss

Example input files are made available with this release, for demo purposes. 

This input format was chosen because it is similar to the language plat input required for Raphael Finkel and Gregory Stump's Principal Parts Analyzer (https://www.cs.uky.edu/~raphael/linguistics/analyze.html) (Stump and Finkel 2013). The main difference is that a column of frequency counts is required here.

It is also the format that is input to / output by the iterated Bayesian learning model for analogical learning built by Rob Reynolds, Andrea D. Sims, and Jeff Parker (Parker et al. forthcoming). The code for that model is available at: https://github.com/reynoldsnlp/bayes-morph-soc-net . (The demo files available with this release are actual output files from the Bayesian model.)

####################

OUTPUT FILES

The code outputs the following ("language" is replaced with input file name):

./network_properties/networkProperties_language.txt: calculations of node-level properties: number of components in a graph and which component a node belongs to, node degree, mean weight of a node's edges, sum weight of all of node's edges, node (local) clustering coefficient, mean shortest path length between node and all other nodes (weighted by type frequency and unweighted). resid_weight is edge weight residualized on degree and clustering coefficient. This number is useful for regression, since edge weight is heavily correlated with degree and clustering coefficient. The file also includes class labels and type frequency counts from input language plat.

./network_properties/modules_language.txt: Calculations of graph modularity and role-to-role connectivity profiles of nodes: how many modules there are in the graph (based on simulated annealing), which module a node belongs to, its connectivity score, participation score, and role. These numbers are calculated using the netcarto package (Guimera and Amaral 2005a, b); see https://amaral.northwestern.edu/resources/software/netcarto for more details.

./entropy/language/language_entropy_weighted_Rcalc.txt: Calculations of probability (for individual inflected forms, and joint prob for pairs of forms) and entropy (for paradigm cells and pairs of paradigm cells). Probabilities/entropy are calculated based on type-frequency-weighted forms. 

./entropy/language/language_entropy_nonweighted_Rcalc.txt: Same as above, except calculations are unweighted -- probability is based on the number of *classes* a form appears in.

./entropy/language/language_entropy_summary.txt: A summary file for the input data set as a whole, with mean entropy and mutual information (MI) values, calculated over the values for the individual cells. Calculated from weighted and nonweighted files above.

./subtract_classes/language/drop1_language.txt: Contains the original input language plat, but with one class removed ("drop class" file). The name of the file corresponds to the row that was removed from the original input. There is one such file for each row of the original file.

./subtract_classes/language/drop1_language_entropy.txt: Entropy calculations based on the language plat with one class removed (corresponding to above).

./subtract_classes/language/language_dropSummary.txt: A summary file giving mean values for each of the drop class versions of the inflectional systems, plus a calculation of the difference between the cell-pairwise conditional entropy of the full system and the conditional entropy of the system with one class removed (drop class system).

graphs/language/: the code outputs 22 total graphs for each input file. 
- A range of network graphs reflecting: 
-- either all edges visualized ("full") or with weak edges trimmed. An edge is 'weak' if edge weight represents half, or fewer, of all paradigm cells. (For example, in a paradigm with 12 cells, maximum possible edge weight is 11, and the edge is considered weak if an edge connecting two classes has a weight of 6 or less.)
-- nodes colored uniformly ("noGroups"), according to the betweenness centrality of nodes ("betweenness", darker = more central), or according to which module they fall into ("modules").

- A variety of graphs showing the distributions and correlations of network properties

- An effect panel plot of log type frequency (log node size), degree, residualized edge weight, and local clustering coefficient as main effect predictors of entropy difference.

####################

MISCELLANEOUS NOTES

This code is written in R, which means that it can be CPU-bound. In general, the scripts run well on a local machine, but large, dense networks can take a long time to run. The hold up is at two places in the code: generating network graphs, and calculating entropy difference. Run time is primarily a function of the number of nodes in a network (which affects esp. entropy difference calculations) and edges in a network (which affects esp. how long it takes to generate network graphs). (Among the original datasets, Seri takes a particularly long time to generate network graphs, and Icelandic takes a particularly long time to make entropy difference calculations.) Maybe someday I'll get around to fixing this...

Since entropy difference does not need to be calculated each time a small change is made, the code checks whether "drop class" (e.g. subtract_classes/language/drop1_...) files already exist for the dataset. If yes, calculation of entropy difference is skipped. Note that this will cause problems (only) if the content of a language plat is changed without renaming. In this case, to run new entropy difference calculations, the easiest thing is to simply delete or rename the directory with the previous calculations.

####################

LICENSE

This release is licensed under GNU General Public License 3.0.

####################

CITATION

If you use or modify this code, please cite it as:

Sims, Andrea D., with contribution by Jeff Parker. 2019. Inflectional networks: Resources for graph-theoretic analysis of linguistic morphology, v. 1.0.0 [24 December 2019]. DOI: 10.5281/zenodo.3594436 

and

Sims, Andrea D. 2020. Inflectional networks: Graph-theoretic tools for inflectional typology. Proceedings of the Society for Computation in Linguistics 3: Article 10, 88-98. https://scholarworks.umass.edu/scil/vol3/iss1/10

####################

CONTACT

This code is offered without any support or warranty. However, comments and suggestions are welcomed. Email Andrea Sims: sims.120@osu.edu

####################

ACKNOWLEDGMENTS

A big 'thank you' to Jeff Parker, who in addition to being a good collaborator in general, contributed pieces of the scripts, provided the Russian data, and arranged access to some of the other datasets that are made available with this code (and used for Sims and Parker (2016) and Sims (2020)). Thanks to Matthew Baerman for providing the Kadiwéu, Nuer, Seri, and Võro data. Thanks to Raphael Finkel and Greg Stump for the French and Icelandic data.

####################

REFERENCES

Ackerman, Farrell, James P. Blevins, and Robert Malouf. 2009. Parts and wholes: Implicative patterns in inflectional paradigms. In Analogy in grammar: Form and acquisition, ed. by James P. Blevins and Juliette Blevins, 54-82. Oxford: Oxford University Press.

Ackerman, Farrell and Robert Malouf. 2013. Morphological organization: The Low Conditional Entropy Conjecture. Language 89(3): 429-464.

Csardi, G. and T. Nepusz. 2006. The igraph software package for complex network research. InterJournal, Complex Systems 1695. http://igraph.org.

Guimerà, Roger & Luís A.N. Amaral. 2005a. Functional cartography of complex metabolic networks, Nature 433: 895-900.

Guimerà, Roger & Luís A.N. Amaral. 2005b. Cartography of complex networks: modules and universal roles, J. Stat. Mech.-Theory Exp., art. no. P02001.

Parker, Jeff, Robert Reynolds, and Andrea D. Sims. forthcoming. The role of language-specific network properties in the emergence of inflectional irregularity. In Morphological typology and linguistic cognition, ed. by Andrea D. Sims, Adam Ussishkin, Jeff Parker, and Samantha Wray. Cambridge: Cambridge University Press.

R Core Team. 2017. R: A language and environment for statistical computing: R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org

Sims, Andrea D. and Jeff Parker. 2016. How inflection class systems work: On the informativity of implicative structure. Word Structure 9(2): 215-239.

Stump, Gregory T. and Raphael A. Finkel. 2013. Morphological typology: From word to paradigm. Cambridge: Cambridge University Press.
