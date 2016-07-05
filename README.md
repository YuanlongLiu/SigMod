SigMod: an exact and efficient method to identify a strongly interconnected disease-associated module in gene connectome

Motivation: Apart from single marker-based tests classically used in genome-wide association studies (GWAS), network-assisted analysis has become a promising approach to identify the joint effect of multiple genetic factors on disease. To date, most network-assisted methods aim at finding genetic units connected in a background network, whatever their density or strength of connections. This can hamper the findings as sparse connections are non-robust against noise from either the GWAS results or the network resource.

Results: SigMod is a novel and efficient method integrating GWAS results and gene network to identify a strongly interconnected gene module enriched in high association signals. It is formulated as a binary quadratic optimization problem, which can be solved exactly through min-cut algorithms. Compared to existing methods, SigMod has several desirable properties: (i) edge weights quantifying evidence of connections between genes are taken into account (ii) the selection path corresponds to variation of sparsity parameter can be computed rapidly (iii) the identified gene module is strongly interconnected, hence includes genes of highly functional relevance and (iv) it is robust against noise from either the GWAS results or the network resource. We applied SigMod to both simulated and real data. It outperforms state-of-the-art network-assisted methods in identifying disease-association genes.

SigMod tries to select gene module by solving the optimization function: obj(u) = c'u + \lambda*u'*A*u + \eta*|u|, where u is a vector of binary variables to be optimized. u=1 means a gene is selected and u=0 otherwise. c is a vector of gene scores, can be computed from the GWAS results using gene-based methods. 'A' is the adjacence matrix of the  gene network. \lambda and \eta are tuning parameters. The first term c'u of obj(u) encourages high score genes to be selected; the second term 


More information wil be added soon.
Contact: yuanlong.liu@inserm.fr
