### R code from vignette source './SigMod_manual.Rnw'

###### CREATE A FOLDER UNDER SigMod V4 to save output files
#####mkdir Sigmod_Output
###### cd SigMod_Output
##### R
###################################################
### code chunk number 1: SigMod_manual.Rnw:41-42
###################################################
options(width = 60)


###################################################
### code chunk number 2: SigMod_manual.Rnw:51-55
###################################################
# load the igraph package required by SigMod
require( igraph )
if(packageVersion("igraph") < "1.0.0") {
    stop("Need to install igrah version 1.0.0 or above")
}


###################################################
### code chunk number 3: SigMod_manual.Rnw:58-66
###################################################
# load SigMod functions
# setwd('path/to/SigMod')
source('../R/construct_scored_network.R')
source('../R/construct_string_net.R')
source('../R/convert_pvalues_2_scores.R')
source('../R/remove_isolated_nodes.R')
source('../R/SigMod_core.R')
source('../R/selection_path.R')
source('../R/selection_path_fast.R')
source('../R/SigMod_bisection.R')
source('../R/SigMod_v00.R')
source('../R/determine_selection_fast.R')
load('../data/OR_genes/OR_genes.Rdata')


###################################################
### code chunk number 4: SigMod_manual.Rnw:91-96
###################################################
# read gene score data into R as a data frame
gene_scores_file = '../data/gene_p_score/gene_scores.tab'
gene_scores = read.table(gene_scores_file, header = TRUE)
## Convert all gene names to upper case to prevent case mismatch:
gene_scores$gene = toupper(gene_scores$gene)
head(gene_scores)


###################################################
### code chunk number 5: SigMod_manual.Rnw:100-109
###################################################
# read gene network data into R
network_data_file = '../data/network/PINA_net.tab'
network_data = read.table(network_data_file,
                          header = TRUE, fill=TRUE,
                          stringsAsFactors = FALSE, 
                          sep='\t', quote='"')

## Convert all gene names to upper case to prevent case mismatch:
network_data[,3] = toupper(network_data[,3])
network_data[,4] = toupper(network_data[,4])


###################################################
### code chunk number 6: SigMod_manual.Rnw:113-117
###################################################
# construct a node-scored gene network
net = construct_scored_net(network_data, interaction_indices=c(3,4),
                           weight_index=NULL, gene_scores )
## an overview of the constructed net:
net 

## SAVE THE SCORED NETWORK included in the cosntruct_scored net function 
## will be saved in SigMod_Output: file name "constructed_scored net"

###################################################
### code chunk number 7: SigMod_manual.Rnw:121-127
###################################################
# If gene p-values are provided, convert them to scores using the p2score function
gene_ps_file = '../data/gene_p_score/gene_ps.tab'
gene_ps = read.table(gene_ps_file, header = TRUE)
head(gene_ps)
## convert gene p-values to scores using the p2score function:
gene_scores = p2score( gene_ps )
head(gene_scores)


###################################################
### code chunk number 8: SigMod_manual.Rnw:196-199 (eval = FALSE)
###################################################
# run the SigMod_bisection function to identify module(s)
sink("res_info_SigMod_bisection.txt")
res_info = SigMod_bisection( net,
                              lambda_min=0, lambda_max=1,
                              nmax=300, maxjump=10)
sink()							  
## SAVE THE OUTPUTS at each step in the bisecton function: information ( lamba, size of S, jump..) appearing on the scrren will be saved in a log file
# will be saved in SigMod_Output: file names "SigMod_computation_details.R;selected_genes.tab;selected_genes_next.tab;selected_module;selected_module_next"
### IT WILL BE NECESSARY TO ANNOTATE the list of genes in selected_genes.tab with crom, position (strat,end), p values, best SNP..."
###################################################
### code chunk number 9: SigMod_manual.Rnw:202-203
###################################################
load('./res_info.Rdata')
## LOAD IN DIR WHERE THE DETAILED OUTPUT OF THE SigMod_bisection

###################################################
### code chunk number 10: SigMod_manual.Rnw:206-207
###################################################
# Display the information of res_info at the each step
names( res_info )


###################################################
### code chunk number 11: SigMod_manual.Rnw:211-212
###################################################
# Display the information of optimal module(s):
res_info$opt_module


###################################################
### code chunk number 12: SigMod_manual.Rnw:216-222
###################################################
# Retrieve genes in the selected module
selected_module = res_info$opt_module[['selected']]
selected_genes = V(selected_module)$name
## order the genes by gene score:
selected_genes = selected_genes[order( V(selected_module)$weight, 
                                       decreasing=TRUE ) ]
head(selected_genes)


###################################################
### code chunk number 13: SigMod_manual.Rnw:228-231
###################################################
# plot the selected module / do not show gene names
set.seed(1)
par(mar=c(0,0,0,0))
plot( selected_module, vertex.size = 5, vertex.label=NA )


###################################################
### code chunk number 14: SigMod_manual.Rnw:237-240
###################################################
# plot the selected module / show gene names on the scrren (files saved in bissection function)
set.seed(1)
par(mar=c(0,0,0,0))
plot( selected_module, vertex.size = 5, vertex.label.cex=0.5 )





#  for STRING
###################################################
### code chunk number 15 SigMod_manual.Rnw:304-314 (eval = FALSE)
###################################################
# Construct a scored-network using the construct_string_net function
# If a gene score file is provided
require(STRINGdb)
gene_scores_file = '../data/gene_p_score/gene_scores.tab'
gene_scores = read.table(gene_scores_file, header = TRUE)
string_net = construct_string_net(
  gene_scores = gene_scores,
  evidence2exclude = c('textmining_transferred', 'textmining'),
  genes2exclude = "OR_genes",
  confidence=0.7 )
## an overview of the constructed string_net:
string_net
sink("sting_net_mattrix.txt")
string_net[]
sink()

## Search module(s) on the scored-string-network
sink("res_info_SigMod_bisection_string.txt")
res_info = SigMod_bisection( string_net,
                              lambda_min=0, lambda_max=1,
                              nmax=400, maxjump=10)
sink()
  

