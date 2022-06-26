	##this function computes the selection path for a fixed \lambda, i.e., all distinct selections by varying \eta from eta_min to eta_max
	##for details please refer to Supplementary Notes Section 3.1 'Computing the selection path at any given \lambda value' of the SigMod paper
	##eta_max is set as very big (1E7) by default, which leads to the empty selection	
	##last modified: 03-11-2015

	##compute k(\lambda) for each cut
	##this function calculates the capacity line
	capac_info <- function(net, lambda, eta, silent=0)
	{
		st_info = SigMod_core(net, lambda, eta, silent) 
		change_points = st_info$change_points
		
		cut_info = st_info$st_cut
		st_net = st_info$st_net
		selected = st_info$selected

		slope = length( selected )
		##correct the capacity line to make it convex
		{
			d_s = degree( st_net, v='s' ) ##number of nodes connected with 's'
			delta = -sum( rep(eta, d_s) - change_points[1:d_s] )
			cut_val = cut_info$value - delta
		}

		intercept = cut_val - slope*eta
		capac_line = list(ax=slope, b=intercept)
		return( list(capac_line=capac_line, selected=selected) ) ##ax + b
	}	

	##reg_path computes the regularization path/selection path
	##selection_diff: resolution of the path. Stop computing if two consecutive selections differ no more than "selection_diff" nodes
	##selection_min: do not compute the selections if its nodes is less than selection_min. A proper specification of selection_min can save computational time
	##selection_max: do not compute the selections if its nodes is more than selection_max. A proper specification of selection_max can save computational time
	reg_path <- function( net, lambda, eta_max, eta_min, capac_infos, selection_diff=1, selection_min, selection_max, silent )
	{	
		selected = list()
		selected[[1]] = capac_infos[[1]]$selected
		selected[[2]] = capac_infos[[2]]$selected
		
		##stop if the difference between selection is smaller than selection_diff
		# cat('\n', selection_diff, '\n')
		if( abs( length(selected[[1]]) - length(selected[[2]]) ) <= selection_diff ) return( list( selected[[1]], selected[[2]] ) )
		if( max( length(selected[[1]]) , length(selected[[2]]) ) <= selection_min  ) return( list( selected[[1]], selected[[2]] ) )
		if( min( length(selected[[1]]) , length(selected[[2]]) ) >= selection_max  ) return( list( selected[[1]] ) )

		##divide and conquer
		# cat('Compute selections with size falling between:', length(selected[[1]]), '||', length(selected[[2]]), '\n')
		cat('(', length(selected[[1]]), ',', length(selected[[2]]), ') ', sep='')

		etas = c( eta_max, eta_min )
		
		cut_lines = list()
		cut_lines[[1]] = capac_infos[[1]]$capac_line
		cut_lines[[2]] = capac_infos[[2]]$capac_line
		
		##else sparsity_mid equals the intersection of two cut_val lines
		eta_mid = ( cut_lines[[2]]$b - cut_lines[[1]]$b ) / ( -cut_lines[[2]]$ax + cut_lines[[1]]$ax )
		##no more selections between
		criterion_b = any( abs(c(eta_mid - etas[1], eta_mid - etas[2])) < 1E-7 )
		if( criterion_b )
		{ #cat(min(abs(c(eta_mid - etas[1], eta_mid - etas[2]))), '\n'); 
		return( list( selected[[1]], selected[[2]] ) ) }
		
		##error if the intersection falls out of the original interval
		non_concave = ( eta_mid > max( etas[1], etas[2] ) ) | ( eta_mid < min( etas[1], etas[2] ) )
		if( is.na(non_concave) ) { cat('incorrect slope of capac_line\n'); return(-1) }
		# if( non_concave ) {cat('capac_line neither concave nor convex\n'); return(-1)}
		if( non_concave ) {cat(eta_min, eta_mid, eta_max); return(-1)}
		
		capac_info_mid = capac_info( net, lambda, eta=eta_mid, silent )
		
		# net_head = induced.subgraph( net, capac_info_mid$selected ):: capac_info not updated
		reg_paths_head  = reg_path( net, lambda, etas[1], eta_mid, list(capac_infos[[1]], capac_info_mid), selection_diff, selection_min, selection_max, silent )
		reg_paths_tail  = reg_path( net, lambda, eta_mid, etas[2], list(capac_info_mid, capac_infos[[2]]), selection_diff, selection_min, selection_max, silent )
	
		return( c( reg_paths_head, reg_paths_tail ) )
	}

	selection_path <- function( net, lambda, eta_max, eta_min, selection_diff, selection_min, selection_max, silent )
	{
		attribute_names = vertex_attr_names(net)
		for(name in setdiff( attribute_names, c('name', 'weight') )) net = remove.vertex.attribute(net, name)
		
		capac_infos = list()
		capac_infos[[1]] = capac_info(net=net, lambda, eta=eta_max, silent )
		capac_infos[[2]] = capac_info(net=net, lambda, eta=eta_min, silent )
		
		cat('\nBEGIN-------------------------------- \n')
		cat('Maximum number of nodes selected:', length(capac_infos[[2]]$selected), '\n')
		cat('Computing modules with size falling between:\n')
		result = reg_path( net, lambda, eta_max, eta_min, capac_infos, selection_diff, selection_min, selection_max, silent )
		
		cat('\n----------------------------------END\n')
		result = result[!duplicated(sapply(result,length))]
		
		result = result[sapply(result,length) >= selection_min]
		# result = result[sapply(result,length) <= selection_max]  ##	commented by Yuanlong LIU, 05-11-2017. This line removes too much of the modules
		return( result[ order(sapply( result, length )) ] )
	}
	
	
	