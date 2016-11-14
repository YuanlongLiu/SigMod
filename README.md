SigMod: an exact and efficient method to identify a strongly interconnected disease-associated module in gene connectome
==================================

**Motivation:** Apart from single marker-based tests classically used in genome-wide association studies (GWAS), network-assisted analysis has become a promising approach to identify the joint effect of multiple genetic factors on disease. To date, most network-assisted methods aim at finding genetic units connected in a background network, whatever their density or strength of connections. This can hamper the findings as sparse connections are non-robust against noise from either the GWAS results or the network resource.

**Results:** We developed SigMod, a novel and efficient method that integrates GWAS results with gene network to identify a strongly interconnected gene module enriched in high association signals. SigMod identifies the gene module by solving the optimization problem: 

* `maxmize: obj(u) = c'u + \lambda u'Au + \eta |u|,`

where `u` is a vector of binary variables to be optimized. `u=1` means a gene is selected and `u=0` otherwise. `c` is a vector of gene scores, can be computed from the GWAS results using gene-based methods. `A` is the adjacency matrix of the gene network. `\lambda` and `\eta` are tuning parameters. The meaning of each term in `obj(u)` is:

* `c'u`: the combined effect of all selected genes
* ` u'Au`: the overall connection strength among selected genes
* `|u|`: the number of selected genes

Therefore the first term `c'u` encourages high score genes to be selected; the second term `\lambda u'Au` encourages the strong interconnection among selected genes; the third term `\eta |u|` limits the number of genes to be selected. We showed that the maximization of `obj(u)` can be solved exactly and efficiently through min-cut algorithms. Compared to existing methods, SigMod has several desirable properties: (i) edge weights quantifying evidence of connections between genes are taken into account; (ii) the selection path corresponds to variation of sparsity parameter can be computed rapidly; (iii) the identified gene module is strongly interconnected, hence includes genes of highly functional relevance; and (iv) it is robust against noise from either the GWAS results or the network resource. We applied SigMod to both simulated and real data. It outperforms state-of-the-art network-assisted methods in identifying disease-association genes.


## Usage
Usage and examples will be added soon

* Author: Yuanlong LIU
* Affiliation: French National Institute of Health and Medical Research, Unit 946,  Paris, France
* Mail: yuanlong.liu@inserm.fr
