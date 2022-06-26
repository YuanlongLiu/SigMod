	##	this function used the bisection algorithm to find the optimal module
	##	updated 02-01-2018
	
	path_fw_weighted <- function( net, lambda, edge_limit ) ## edge weight added 23_03_2016
	{
		strength_net = strength( net )
		genes = V(net)$name
		gene_len = length(genes)

		gene_weights = V(net)$weight
		names(gene_weights) = genes
		cands = genes
		regul_path = character()

		# net = remove.vertex.attribute(net, 'name')
		
		first = genes[ which.max( gene_weights ) ] ## the first gene to select
		cands = setdiff( cands, first )
		deg_ini = rep( 0, gene_len )
		names(deg_ini) = genes
		next_in = first
		regul_path = c( regul_path, next_in )

		regul_path_subnet = induced.subgraph( net, regul_path )
		
		while( ( ecount(regul_path_subnet) < edge_limit ) & ( length(regul_path) < vcount( net ) ) )
		{
			previous_in = next_in

			pair = character(2*length(cands))
			pair[ (1:length(pair)) %% 2 == 1  ] = cands
			pair[ (1:length(pair)) %% 2 == 0  ] = previous_in
			e_id <- get.edge.ids(net, pair)
			
			temp = numeric( length(cands) )
			temp[ e_id!=0 ] <- edge_attr(net, "weight", e_id)	## this takes time
			degs_suf = deg_ini[ cands ] + temp
		
			next_in = cands[ which.max( 2*lambda*degs_suf + gene_weights[cands] )]
			regul_path = c( regul_path, next_in )
			cands = setdiff( cands, next_in )
			deg_ini = degs_suf
			regul_path_subnet = induced.subgraph( net, regul_path )

			# cat(ecount(regul_path_subnet), '\n')
		}
		return( regul_path )
	}	
	
