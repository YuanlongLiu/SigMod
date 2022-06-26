## The following code performs SigMod analysis using the STRING network
## Created by YL, 20-12-2017


## Please first create a working folder under SigMod_v4, then change the directory to that folder
## Then run the following script to conduct the analysis 

## Load the igraph package
require( igraph )
if(packageVersion("igraph") < "1.0.0") {
    stop("Need to install igrah version 1.0.0 or above")
}

## Load the SigMod functions
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
source('../R/p_association_cir.R')
source('../R/p_enrichment_of_score.R')
load('../data/OR_genes/OR_genes.Rdata')

## Construct a scored-string-network
require(STRINGdb)
gene_ps_file = '../data/gene_p_score/gene_ps.tab'
gene_ps = read.table(gene_ps_file, header = TRUE)
## construct the string-based scored-network:
string_net = construct_string_net(
  gene_ps = gene_ps,
  evidence2exclude = c('textmining_transferred', 'textmining'),
  genes2exclude = "OR_genes",
  confidence=0.7 )

sink("sting_net_matrix.txt")
string_net[]
sink()
## Search module(s) on the scored-string-network
# try nmax= 300 350 400 450 500
# and try max jump = 5 10 15 20
sink("res_info_SigMod_bisection_string.txt")
res_info = SigMod_bisection( string_net,
                              lambda_min=0, lambda_max=1,
                              nmax=300, maxjump=10)
sink()

