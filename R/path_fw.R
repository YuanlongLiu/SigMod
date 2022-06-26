	##	this function used the bisection algorithm to find the optimal module
	##	updated 02-01-2018
	
	path_fw <- function( net, lambda, edge_limit )
	{
		genes = V(net)$name
		gene_weights = V(net)$weight
		names(gene_weights) = genes
		cands = genes
		regul_path = character()
		
		first = genes[ which.max( gene_weights ) ] ## the first gene to select
		cands = setdiff( cands, first )
		deg_ini = rep( 0, length(genes) )
		names(deg_ini) = genes
		next_in = first
		regul_path = c( regul_path, next_in )
		regul_path_subnet = induced.subgraph( net, regul_path )
		
		while( ( ecount(regul_path_subnet) < edge_limit ) & ( length(regul_path) < vcount( net ) ) )
		{
			previous_in = next_in
			degs_suf = deg_ini[ cands ] + ( cands %in% ( ego( net, 1, previous_in, mindist=1 )[[1]] )$name )
			
			next_in = cands[ which.max( 2*lambda*degs_suf + gene_weights[cands] )]
			regul_path = c( regul_path, next_in )
			cands = setdiff( cands, next_in )
			deg_ini = degs_suf
			regul_path_subnet = induced.subgraph( net, regul_path )
			# cat( ecount(regul_path_subnet), '\n' )
		}
		return( regul_path )
	}
	
	
