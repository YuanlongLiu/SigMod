	##	this function used the bisection algorithm to find the optimal module
	
	SigMod_bisection <- function( net, lambda_min=0, lambda_max=1, eta_max=1E7, size_diff=1, size_min=0, nmax=300, maxjump=30, epsilon=1E-5, cp=1000, silent=1 )
	{
		if( max( size_diff, size_min, nmax, maxjump ) > vcount(net) ) stop(paste('Non of the parameters size_diff, size_min, nmax, maxjump should not be greater than the total number of nodes in the whole network:', vcount(net)))

		lambdas = numeric()	##stores the lambda values that have been computed
		reses = list()	##stores the selection path of each lambda
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
			tmp = selection_path_fast( net=net, lambda=lambda, eta_max=1E7, eta_min=eta_min, selection_diff=size_diff, selection_min=size_min, selection_max = nmax, silent )
			tmp = tmp[order( sapply( tmp, length ) )]
			keep_index = min( which((sapply(tmp, length) >= nmax)) )
			tmp = tmp[1:keep_index]
			
			reses[[i]] = tmp
			lambdas[i] = lambda
			
			maxjump_obs = max( diff(sapply( tmp, length )) )
			cat('Max size jump in this path:', maxjump_obs, '\n')
			flag = maxjump_obs > maxjump
			if( (i==1) && (flag == 0)) stop(paste('\nThe observed \'maxjump\' is:', maxjump_obs, '\nPlease decrease the \'maxjump\' parameter to at least:', maxjump_obs, '\nOr maybe the \'lambda_max\' parameter you set is too small'))
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

		opt_path = reses[[opt_index]] ## the optimal path
		if( length( tail( opt_path, 1 )[[1]] ) == nmax )
		{
			opt_mod_raw = tail( opt_path, 1 )
			names( opt_mod_raw ) = 'opt_module'
		}
		
		if( length( tail( opt_path, 1 )[[1]] ) != nmax )
		{
			opt_mod_raw = tail( opt_path, 2 )
			names( opt_mod_raw ) = c('selected', 'selected_next')
		}		
		
		
		opt_module = lapply( opt_mod_raw, function(v) remove_isolated( induced.subgraph( net, v ) ) )
		info = list( opt_module=opt_module, lambdas=lambdas, paths=reses, opt_index=opt_index, net=net )
		save(info, file='./SigMod_computation_details.Rdata')
		obtain_results( info, net )
		return( info )
	}
	
	obtain_results <- function( SigMod_computation_details, net )
	{
		
		res_info = SigMod_computation_details
		flag_next = length(res_info$opt_module) > 1
		flag_p = 'p' %in% vertex_attr_names(net) ## whether the network has a p-value attributes

		# Retrieve genes in the selected module
		selected_module = res_info$opt_module[[1]]
		selected_genes = V(selected_module)$name
		## order the genes by gene score:
		selected_genes = selected_genes[order( V(selected_module)$weight, 
											   decreasing=TRUE ) ]

		if(flag_p)
		{
			selected_genes_ps = V(net)[selected_genes]$p
			selected_genes_info = data.frame(gene = selected_genes, p = selected_genes_ps)
		}
	
		if(!flag_p)
		{
			selected_genes_info = data.frame(gene = selected_genes)
		}
		write.table(selected_genes_info, file='./selected_genes.tab', quote=FALSE, row.names=FALSE, col.names=TRUE)
		
		# plot the selected module
		main = paste('Number of nodes:', vcount(selected_module), '\n', 'Number of interactions:', ecount(selected_module))
		pdf('./selected_module.pdf')
		par(mar=c(0,0,2,0))
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label=NA, main=main)
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label.cex=0.5, main=main)
		dev.off()
		
		if( !flag_next ) return()
		
		
		selected_module = res_info$opt_module[[2]]
		selected_genes = V(selected_module)$name
		## order the genes by gene score:
		selected_genes = selected_genes[order( V(selected_module)$weight, 
											   decreasing=TRUE ) ]

		if(flag_p)
		{
			selected_genes_ps = V(net)[selected_genes]$p
			selected_genes_info = data.frame(gene = selected_genes, p = selected_genes_ps)
		}
	
		if(!flag_p)
		{
			selected_genes_info = data.frame(gene = selected_genes)
		}
		write.table(selected_genes_info, file='./selected_genes_next.tab', quote=FALSE, row.names=FALSE, col.names=TRUE)
		
		# plot the selected module
		main = paste('Number of nodes:', vcount(selected_module), '\n', 'Number of interactions:', ecount(selected_module))
		pdf('./selected_module_next.pdf')
		par(mar=c(0,0,2,0))
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label=NA, main=main)
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label.cex=0.5, main=main)
		dev.off()
		cat('\nResults saved under:\n', getwd(), '\n')
	}
	
