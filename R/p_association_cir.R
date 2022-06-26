	##  this function evaluates whether the the selected gene module is significantly associated with a disease status using circular genomic permutation of SNP p-values 
	##  Developed by Yuanlong LIU, contact: yuanlong.liu@inserm.fr
	##  Latest update: 20-12-2017
	
	##	If you this method, please cite in your publication: Liu, Y., et al. "Network-assisted analysis of GWAS data identifies a functionally-relevant gene module for childhood-onset asthma." Scientific Reports 7.1 (2017): 938.
	
	##	this sub-function constructs the 'circle' to perform fastCGP, which includes information of ordered snp p-values, and their associated genes 
	fastCGP_cir_info <- function( snp2gene, snp_chr_pos_p )
	{

		genes = unique( snp2gene$gene )
		mapped_snp  = unique( snp2gene$SNP ) ##snps in genes
		all_snps = unique( snp_chr_pos_p$SNP ) ##the complete snp list
		if( any(duplicated( all_snps )) ) stop('\nSome SNPs appear more than once in your snp_chr_pos_p file\nPlease remove duplicates first')
		
		mapped_snp_nin_all_snps = setdiff( mapped_snp, all_snps )
		diff_len = length(mapped_snp_nin_all_snps) ##SNPs mapped to gene but do not have a position in snp_chr_pos_p file
		if( diff_len > 0 ) stop('\n', len, ' SNPs are mapped to genes in your snp2gene file \nBut their chromosomal or p-value information are missed in your snp_chr_pos_p file')

		genes_count = length(genes)
		snps_count = length(all_snps)
		mapped_snp_count = length(mapped_snp)
		# cat('\nNumber of genes in your snp2gene file:', genes_count, '\n')
		# cat('Number of SNPs in your snp_chr_pos_p file:', snps_count, '\n')
		# cat('Number of SNPs mapped to genes:', mapped_snp_count, '\n')
		
		snps_p = snp_chr_pos_p$p
		names( snps_p ) = snp_chr_pos_p$SNP ##named p values for complete SNP list
		
		snp_prk = rank( snps_p, ties.method='random' ) ##p value rank of snps

		snps_pos_table = 1:snps_count
		names(snps_pos_table) = snp_chr_pos_p[ order(snp_chr_pos_p$chr, snp_chr_pos_p$pos), ]$SNP ##place snps on the table according to their position on the chromosome
		
		table_rk = snp_prk[ names(snps_pos_table) ] ##rank of the snps on the table; rank(c(1,2)): 1 2
		
		gene_snps_poses = tapply( snps_pos_table[ snp2gene$SNP ], snp2gene$gene, function(v) {return(v)} ) ##position of snps on table of gene
		
		if(any(is.na(unlist(gene_snps_poses, recursive=TRUE)))) {cat( 'NA information from SNP table again!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' ); return(0)}
	
		lens = table( snp2gene$gene )[1:genes_count] ##number of snps in each gene
		gene_minrk = as.list(tapply( snp_prk[ snp2gene$SNP ], snp2gene$gene, min ))##best rank of each gene

		snp2gene_p = snp2gene; snp2gene_p$p = snps_p[ snp2gene$SNP ]
		gene_minP = tapply(snp2gene_p$p, snp2gene_p$gene, min) ##best rank of each gene
		
		return( list(gene_snps_poses=gene_snps_poses, gene_lens = lens, gene_minrk=gene_minrk, gene_minP=gene_minP, table_rk=table_rk) )	
	}

	## the function that computes gene p-value using fastCGP_v1 (the same)
	fastCGP_slim <- function( snp2gene, snp_chr_pos_p, genes2compute ) ## this function is the same with the fastCGP_lib
	{
		
		genes = unique( snp2gene$gene )
		unmapped_genes = setdiff( genes2compute, genes )
		if( length( unmapped_genes ) > 0 ) cat( sep='', '\n', length( unmapped_genes ), ' genes in your genes2compute file do not have SNP_to_gene mapping information\nThey are removed for downstream computation\n' )
		
		cir_info = fastCGP_cir_info( snp2gene, snp_chr_pos_p )

		sorted_genes = names(sort(unlist(cir_info$gene_minrk), decreasing=TRUE)) ## sort genes by their p-value rank. More significant genes will be ranked at the tail
		genes2compute = intersect( sorted_genes, genes2compute ) ##ordered according to sorted_genes
		
		x = cir_info$gene_minrk[genes2compute]
		if(!all(x==cummin(x))) stop('\nThe genes to be computed need to be ranked by minP') ##check whether genes are ranked by minP; this is required for the purpose of increasing efficiency
		table_len = length( cir_info$table_rk ) ##total number of SNPs on the circle
		pos2compare = cbind(rk=unname(cir_info$table_rk), pos=1:table_len)
		sig_pos_ini = pos2compare ## the rank at each position
		gene_lens = cir_info$gene_lens[genes2compute]
		
		genes_p =list()
		i=1
		
		# cat('\n\nNumber of genes already computed:\n')
		##compute corrected gene-level p-value gene by gene
		for(gene in genes2compute)
		{
			gene_len = gene_lens[gene]
			
			if( gene_len == 1 ) { genes_p[[ gene ]] = cir_info$gene_minP[ gene ]; next } ##if the gene only contains one SNP, the p-value will be that SNP p-value
			
			thresh = cir_info$gene_minrk[ gene ] ##Test statistics // cir_info$gene_minrk['GSDMB']: 1
			sig_pos = sig_pos_ini[which(sig_pos_ini[,'rk'] <= thresh), 'pos' ] ## CONFRIMED
			if(length(sig_pos) ==1) { genes_p[[gene]]= ( gene_len + 1 )/( table_len + 1 ); next}  ##only extreme cases, that a gene contains the best snp in all_snps 
			gaps = c( diff( sig_pos ) - 1, (table_len - sig_pos[length(sig_pos)]) + (sig_pos[1] - 1) ) ##CONFRIMED
			
			big_gaps = gaps[gaps >= gene_len]
			gap_width = big_gaps - gene_len + 1 ##CONFRIMED
			genes_p[[gene]]  = 1 - sum(gap_width) / ( table_len + 1)
			
			if( i %% 100 ==0) cat(i, '|')
			sig_pos_ini = pos2compare[ sig_pos, ] ##because genes are ranked according to minP
			i=i+1
		}

	
		genes_p_frame = data.frame( gene=names( genes_p ), p=unname(unlist( genes_p )) )
		# write.table(genes_p_frame, file='./computed_gene_p_values.tab', row.names=FALSE, sep='\t', quote=FALSE)
		
		return( genes_p_frame )
	}
	
	## Just to have a look how snp_chr_pos_p looks like // YL 19-12-2017
	# > head(snp_chr_pos_p)
			 # SNP chr       pos        p
	# 1  rs1000000  12 126890980 0.394150
	# 2 rs10000010   4  21618674 0.821190
	# 3 rs10000012   4   1357325 0.623990
	# 4 rs10000013   4  37225069 0.081821
	# 5 rs10000017   4  84778125 0.211210
	# 6  rs1000002   3 183635768 0.305950

	## the function that generate the permuted snp_chr_pos_p information
	permutation <- function( snp_chr_pos_p, shift ) ## this function tries to to permute the SNP p-values for one time
	{
		tmp = snp_chr_pos_p[ order(snp_chr_pos_p$chr, snp_chr_pos_p$pos), ]
		len = nrow(tmp)
		if(shift > len) stop('The position to be shifted in CGP should be less than the total number of SNPs analyzed\n')
		indices = 1:len
		shifted_index = c(tail(indices, shift), head(indices, length(indices) - shift))
		tmp$p = tmp$p[shifted_index]
		return( tmp )
	}
	
	## the function that converts gene p-values to scores
	p2score <- function( gene_ps )
	{
		gene_ps$score = qnorm( 1 - gene_ps$p ) ## convert p-values to upper quantile of normal distribution
		
		neg_inf_index = which( is.infinite( gene_ps$score ) & ( 0 > gene_ps$score ) ) ##	p-values converted to -inf
		pos_inf_index = which( is.infinite( gene_ps$score ) & ( 0 < gene_ps$score ) ) ## p-values converted to +inf
		
		neg_inf_index_len = length(neg_inf_index)
		pos_inf_index_len = length(pos_inf_index)

		if( neg_inf_index_len > 0 ) warning('Warning:\n The p-value of ', neg_inf_index_len, ' genes are converted to a score of negative infinity. We reset their scores as -9 for computation consideration.\n'  )
		if( pos_inf_index_len > 0 ) warning('Warning:\n The p-value of ', pos_inf_index_len, ' genes are converted to a score of positive infinity. We reset their scores as 9 for computation consideration.\n'  )
		
		gene_ps$score[ neg_inf_index ] = -9
		gene_ps$score[ pos_inf_index ] = +9
		# write.table( gene_ps, file='./gene_scores.tsv', sep='\t', quote=FALSE, row.names=FALSE )
		gene_scores = gene_ps
		return( gene_scores )
	}

	## the main function to compute the p-value
	p_association_cir <- function(snp2gene_file, snp_chr_pos_p_file, module_genes_file, nperm=1E5)
	{
		time_0 = proc.time()
		##read files:
		snp2gene = read.table( snp2gene_file, header=TRUE, stringsAsFactors=FALSE, comment.char='' )
		snp_chr_pos_p = read.table( snp_chr_pos_p_file, header=TRUE, stringsAsFactors=FALSE, comment.char='' )
		genes2compute = read.table( module_genes_file, header=FALSE, stringsAsFactors=FALSE, comment.char='' )
		genes2compute = unique( toupper( genes2compute$gene ) )
		
		cat('\nNote: all gene names will be converted to uppercase\n')
		##lowercase to uppercase
		snp2gene$gene = toupper( snp2gene$gene )
		if( !all( genes2compute %in% snp2gene$gene ) ) stop('The function cannot be executed since some genes in your module_genes_file do not exist in the snp2gene_file')

		snp_count = nrow( snp_chr_pos_p )
		set.seed(1)
		shifts = sample(2:(snp_count - 1), size=nperm, replace = FALSE, prob = NULL)
		shifts = c(0, shifts)
		permed_module_scores = numeric()
		i = 0
		cat('\n\n\n\n@@@@-------------------------------- COMPUTATION BEGINS --------------------------------@@@@\n')

		cat(c('n_perm', 'module_score'), '\n', file='./p_association_cir.log', append=FALSE, sep='\t')
		for( shift in shifts )
		{
			cat(i, ' ')
			permed_snp_chr_pos_p = permutation( snp_chr_pos_p, shift )
			permed_gene_ps = fastCGP_slim( snp2gene, permed_snp_chr_pos_p, genes2compute )
			permed_gene_scores = p2score( permed_gene_ps )
			permed_module_scores[i+1] = sum( permed_gene_scores$score )
			cat(c(i, permed_module_scores[i+1]), '\n', file='./p_association_cir.log', append=TRUE, sep='\t')
			i = i + 1
		}
		# cat('\n\n\n\n@@@@--------------------------------- COMPUTATION ENDS ---------------------------------@@@@\n')
		
		delta_time = proc.time() - time_0
		minutes = format( round( delta_time[3] / 60, 1), nsmall=1)
		cat('\n\nTotal computational time:', minutes, 'minutes', '\n')
		p_association_cir = sum(permed_module_scores[1] <= permed_module_scores) / nperm
		cat('p_association_cir', '\n', p_association_cir, file='./p_association_cir_result.txt', append=FALSE, sep='')

		return( permed_module_scores )
	}
	





