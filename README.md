SigMod: an exact and efficient method to identify a strongly interconnected disease-associated module in gene connectome
==================================

Motivation: Apart from single marker-based tests classically used in genome-wide association studies (GWAS), network-assisted analysis has become a promising approach to identify the joint effect of multiple genetic factors on disease. To date, most network-assisted methods aim at finding genetic units connected in a background network, whatever their density or strength of connections. This can hamper the findings as sparse connections are non-robust against noise from either the GWAS results or the network resource.

Results: We developed SigMod, a novel and efficient method that integrates GWAS results with gene network to identify a strongly interconnected gene module enriched in high association signals. SigMod identifies the gene module by solving the optimization problem: 

* `maxmize: obj(u) = c'u + \lambda u'Au + \eta |u|,`

where `u` is a vector of binary variables to be optimized. `u=1` means a gene is selected and `u=0` otherwise. `c` is a vector of gene scores, can be computed from the GWAS results using gene-based methods. `A` is the adjacency matrix of the gene network. `\lambda` and `\eta` are tuning parameters. The meaning of each term in `obj(u)` is:

* `c'u`: the combined effect of all selected genes
* ` u'Au`: the overall connection strength among selected genes
* `|u|`: the number of selected genes

Therefore the first term `c'u` encourages high score genes to be selected; the second term `\lambda u'Au` encourages the strong interconnection among selected genes; the third term `\eta |u|` limits the number of genes to be selected. We showed that the maximization of `obj(u)` can be solved exactly and efficiently through min-cut algorithms. Compared to existing methods, SigMod has several desirable properties: (i) edge weights quantifying evidence of connections between genes are taken into account; (ii) the selection path corresponds to variation of sparsity parameter can be computed rapidly; (iii) the identified gene module is strongly interconnected, hence includes genes of highly functional relevance; and (iv) it is robust against noise from either the GWAS results or the network resource. We applied SigMod to both simulated and real data. It outperforms state-of-the-art network-assisted methods in identifying disease-association genes.

==================================

MID for measuring statistical dependence between two random variables.


Summary
-------

An estimation algorithm for MID (Mutual Information Dimension), which measures statistical dependence between two random variables.
This algorithm has the following advantages:

* **Nonlinear dependences** (and also linear dependences) can be measured,
* **Scalable**; the average-case time complexity is O(nlogn), where *n* is the number of data points, and
* **Parameter-free**.

Please see the following article for detailed information and refer it in your published research:

* Sugiyama, M., Borgwardt, K. M.: **Measuring Statistical Dependence via the Mutual Information Dimension**,
	*Proceedings of the 23rd International Joint Conference on Artificial Intelligence (IJCAI 2013)*, to appear.


Installation
------------

The code consists of only one C file "MID.c".
Thus you can use by compiling it, for example, type into your terminal:

	$ gcc -O3 MID.c -o MID


Usage
-----

To calculate MID between two variables, type:

	$ ./MID <input_file>
	
`<input_file>` is a comma-separated text file with two columns without row and column names.
Columns correspond to respective variables.
For example,	

	0.921,0.930
	0.491,0.492
	0.990,0.993
	0.775,0.777
	...
	0.577,0.561

The followings are shown at standard output.

* **dimX** (the information dimension of the first variable)
* **dimY** (the information dimension of the second variable)
* **dimXY** (the information dimension of X and Y)
* **MID** (equivalent to dimX + dimY - dimXY)

Example
-------

	$ gcc -O3 MID.c -o MID
	$ ./MID ./sampledata/linear.csv
	idim_x:  0.994690
	idim_y:  0.994690
	idim_xy: 0.994690
	MID:     0.994690
	$ ./MID ./sampledata/noise.csv
	idim_x:  0.995130
	idim_y:  0.996233
	idim_xy: 1.755107
	MID:     0.236256


Information
-----------

* Author: Mahito Sugiyama
* Affiliation: Machine Learning & Computational Biology Research Group, MPIs Tübingen, Germany
* Mail: mahito.sugiyama@tuebingen.mpg.de
