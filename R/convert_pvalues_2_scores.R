	##	this function used the bisection algorithm to find the optimal module
	
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
		write.table( gene_ps, file='./gene_scores.tsv', sep='\t', quote=FALSE, row.names=FALSE )
		gene_scores = gene_ps
		return( gene_scores )
	}

	