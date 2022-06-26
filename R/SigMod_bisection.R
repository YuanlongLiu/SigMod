	##	this function used the bisection algorithm to find the optimal module
	##	updated 02-01-2018
	
	SigMod_bisection <- function( net, lambda_min=0, lambda_max=1, eta_max=1E7, size_diff=1, size_min=0, nmax=300, maxjump=30, epsilon=1E-5, cp=1000, silent=1 )
	{
		if( max( size_diff, size_min, nmax, maxjump ) > vcount(net) ) stop(paste('Non of the parameters size_diff, size_min, nmax, maxjump should not be greater than the total number of nodes in the whole network:', vcount(net)))
		# edge_limit = min(1e5, ecount(net)) ## 1e5 to 5e5
		edge_limit = min(200000, ecount(net)) ## 1e5 to 5e5

		net_full = net
		
		lambdas = numeric()	##stores the lambda values that have been computed
		paths = list()	##stores the selection path of each lambda
		lambda_mid = lambda_max
		time_0 = proc.time()
		i = 1
		total_count = ceiling( log2((lambda_max - lambda_min)/epsilon) )
		while( i < 100 )
		{
			cat('\n\n\n\n@@@@-------------------------------- Compute next selection path --------------------------------@@@@\n')
			cat('lambda values already computed:', format(lambdas, scientific=TRUE, digits=4))
			cat('\nNumber of lambda values left to compute:', total_count - i + 1 )
			
			time_1 = proc.time()
			
			lambda = lambda_mid
			cat('\nCurrent lambda:', lambda, '\n')
			
			stepwise_flag = ecount( net_full ) > edge_limit
			stepwise_flag_warning = "There are more than 200000 edges in your network. To increase computation speed, a stepwise forward approximation algorithm was first called to preselect a subnetwork that can contribute most to the objective function f(u) = z'u + lambda*u'Au - eta*|u|_0. This algorithm starts from an empty subnetwork and adds node by node that can lead to the maximum gain of f(u). The iteration terminates until the number of edges in the resulted subnetwork reaches edge_limit. Then the exact algorithm is called to compute the selection path over this preselected subnetwork. This approximation+exact strategy was tested to obtain similar results as using the exact algorithm solely\n"
			if(stepwise_flag)
			{
				if(i==1) cat(stepwise_flag_warning)
				if( length(unique(E(net_full)$weight)) != 1 ) nodes_pool = path_fw_weighted( net_full, lambda, edge_limit = edge_limit )
				if( length(unique(E(net_full)$weight)) == 1 ) nodes_pool = path_fw( net_full, lambda, edge_limit = edge_limit )

				if( length( nodes_pool ) < nmax ) stop('The number of genes preselected by path_fw_weighted is too small. Please check the path_fw_weighted algorithm')
				net_pool = induced.subgraph( net_full, nodes_pool)
			}
			else
			{
				net_pool = net_full
			}	
	
	
			change_points = SigMod_core(net_pool, lambda, eta=1E7, silent)$change_points
			cp_len = length( change_points )
			cp = min(cp, cp_len)
			
			while(1){
				cp_info = SigMod_core(net_pool, lambda, eta = change_points[cp], silent)
				num = length(cp_info$selected)
				if(num >= nmax) break	## num is big enough
				if(num <= nmax) cp = min(cp + 100, cp_len)
				cat('cp too small. Increase cp to: ', cp, '\n' )
			}
			
			eta_min = change_points[cp]
			tmp = selection_path_fast( net=net_pool, lambda=lambda, eta_max=1E7, eta_min=eta_min, selection_diff=size_diff, selection_min=size_min, selection_max = nmax, silent )
			tmp = tmp[order( sapply( tmp, length ) )]
			keep_index = min( which((sapply(tmp, length) >= nmax)) )
			tmp = tmp[1:keep_index]
			
			paths[[i]] = tmp
			lambdas[i] = lambda
			
			maxjump_obs = max( diff(sapply( tmp, length )) )
			cat('Max size jump in this path:', maxjump_obs, '\n')
			flag = maxjump_obs > maxjump
			if( (i==1) && (flag == 0)) stop(paste('\nThe observed \'maxjump\' is:', maxjump_obs, '\nPlease decrease the \'maxjump\' parameter to at least:', maxjump_obs, '\nOr increase the \'lambda_max\' parameter'))
			if( flag == 1 ) lambda_max = lambda_mid
			else lambda_min = lambda_mid
			
			if( abs( lambda_min - lambda_max ) < epsilon ) break
			lambda_mid = (lambda_min + lambda_max) / 2
			i = i + 1
			
			delta_time_1 = proc.time() - time_1
			minutes = format( round( delta_time_1[3] / 60, 1), nsmall=1)
			cat('Spend time (computing current path):', minutes, 'minutes', '\n')
			
			delta_time = proc.time() - time_0
			minutes = format( round( delta_time[3] / 60, 1), nsmall=1)
			cat('Spend time (computing all paths):', minutes, 'minutes', '\n')			
		}

		cat('\n\n\n\n@@@@----------------------------------- FINISH ALL COMPUTATION ----------------------------------@@@@\n')
		
		delta_time = proc.time() - time_0
		minutes = format( round( delta_time[3] / 60, 1), nsmall=1)
		cat('Spend time (computing all paths):', minutes, 'minutes', '\n')
		
		opt_index = which(lambdas == lambda_max) ## maybe to change as opt_index = which(lambdas == lambda_min)
		# opt_index = which(lambdas == lambda_min) ## changed by Yuanlong LIU, 07-12-2017

		opt_path = paths[[opt_index]] ## the optimal path
		if( length( tail( opt_path, 1 )[[1]] ) == nmax )
		{
			opt_mod_raw = tail( opt_path, 1 )
			names( opt_mod_raw ) = 'opt_module'
		}
		
		if( length( tail( opt_path, 1 )[[1]] ) != nmax ) ## current and previous
		{
			opt_mod_raw = tail( opt_path, 2 )
			names( opt_mod_raw ) = c('selected', 'selected_next')
		}		
		
		
		opt_module = lapply( opt_mod_raw, function(v) remove_isolated( induced.subgraph( net_full, v ) ) )
		info = list( opt_module=opt_module, lambdas=lambdas, paths=paths, opt_index=opt_index, net=net_full )
		save(info, file='./SigMod_computation_details.Rdata')
		obtain_results( info, net_full )
		return( info )
	}
	

	
