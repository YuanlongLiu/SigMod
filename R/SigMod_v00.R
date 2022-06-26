	# source('../R/construct_scored_network.R')
	# source('../R/construct_string_net.R')
	# source('../R/convert_pvalues_2_scores.R')
	# source('../R/remove_isolated_nodes.R')
	# source('../R/SigMod_core.R')
	# source('../R/selection_path.R')
	# source('../R/selection_path_fast.R')
	# source('../R/SigMod_bisection.R')
	# load('../data/OR_genes/OR_genes.Rdata')
	
	# source_root = '/ext/egea/Yuanlong Liu/SigMod_lib'
	
	# source(file.path(source_root,'header.R'))
		
	# while( substr(as.character(Sys.time()),12,13)!='20' ) {}
	
	compute_selection_path <- function( net, lambda, size_min=0, size_diff=1, silent=1, cp=1000, nmax=500 )
	{
		change_points = SigMod_core(net, lambda, eta=1E7, silent)$change_points
		cp_len = length( change_points )
		cp = min(cp, cp_len)
				
		while(1){
			cp_info = SigMod_core(net, lambda, eta = change_points[cp], silent)
			num = length(cp_info$selected)
			if(num >= nmax) break	## num is big enough
			if(num <= nmax) cp = min(cp + 100, cp_len)
			cat('cp too small. Increase cp to: ', cp, '\n' )
		}
		
		eta_min = change_points[cp]
		tmp = selection_path_fast( net, lambda, eta_max=1E7, eta_min, size_diff, size_min, selection_max = nmax, silent=silent )
		tmp = tmp[order( sapply( tmp, length ) )]
		keep_index = min( which((sapply(tmp, length) >= nmax)) )
		tmp = tmp[1:keep_index]	
	}
	
	

	
	compute_selection_paths <- function( net, lambda_min=0.001, lambda_max=0.1, k=1000, size_min=0, size_diff=1, silent=1, cp=1000, nmax=1000 )
	{
		cat('Begin to compute the selection paths\n')
		cat('For', k, 'lambda values in [', lambda_min, ',', lambda_max, ']\n'  )
		lambdas = seq(lambda_min, lambda_max, length=k)
		paths = list()
		for( i in 1:k )
		{
			lambda = lambdas[k]
			path_k = compute_selection_path( net=net, lambda=lambda, size_min=size_min, size_diff=size_diff, silent=silent, cp=cp, nmax=nmax )
			paths[[i]] = path_k
			cat(i, '\n')
		}
	}
	
	# tmp = compute_selection_paths( net, lambda_min=0.001, lambda_max=0.1, k=1000, size_diff=1, silent=1, cp=1000, nmax=500 )
	
	