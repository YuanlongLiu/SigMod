	##	The SigMod_core function computes the nodes selection by maximizing the objective function: c'f + lambda*f'Af - eta*f'f
	##	Two parameters: interconnectivity parameter: lambda; sparsity parameter: eta
	##	Library dependency: igraph 1.0.0 or above http://igraph.org/r/
	##	Last updated: 03-07-2017
	##	Maintained by Yuanlong LIU, contact: yuanlong.liu@inserm.fr

	SigMod_core <- function(net, lambda, eta, silent=0) 
	{
		require( igraph )
		if(packageVersion("igraph") < "1.0.0") {
			stop("\nNeed to install igraph version 1.0.0 or above\nLink to the igraph package: http://igraph.org/r/")
		}

		v_names = V(net)$name
		if( any(c( 's', 't' ) %in% v_names) ) stop('\nVertex names should be different from \'s\', \'t\' \nThey are reserved names used for constructing the augmented graph\n')
		##nodes should have the 'weight' attribute
		if( !'weight' %in% vertex_attr_names( net ) ) stop('\nVertices should have a \'weight\' attribute\n')
		
		##edges should have the 'weight' attribute
		if( !'weight' %in% edge_attr_names( net ) ) {
		stop('\nEdges should have a \'weight\' attribute\nPlease use \'E(net)$weight = values\' to specify edge weights. Edge weights can be set as 1 for unweighted graph\n') }
		
		if( is_directed( net ) ) {
		stop('\nYour input network is directed. Please convert it an undirected network use the as.undirected function from the igraph package\n') }		
		
		v_weight = V(net)$weight

		d =  strength( net ) ## L = D - W for weighted graph. Modified 29-05-2015 

		##construct the augmented graph
		{
			net = net + 's' + 't'
			net <- add.edges(net, as.vector(rbind('s', v_names)))
			net <- add.edges(net, as.vector(rbind('t', v_names)))

			W = as_adjacency_matrix(net, attr='weight') ##slow
			A = lambda * W
			cf = v_weight + lambda*d - eta ##	favors nodes with high degree
			change_points = sort( pmax(v_weight + lambda*d, 0),  decreasing=TRUE )
			change_points = unname(change_points[change_points > 0])
			# if( eta > max( change_points ) ) { cat('\n Parameter inappropriate. Eta should be in:', range(change_points), '\n'); return(-1) }
				
			A['s', v_names] =  cf*(cf >= 0) ##keep edges on 's' side
			A['t', v_names] = -cf*(cf < 0) ##keep edges on 't' side 
			A[,'s'] = A['s',]
			A[,'t'] = A['t',]
			st_net = graph_from_adjacency_matrix( A, mode="undirected", weighted=TRUE )
		}
		
		if(!silent) cat( '\nmin_cut starts::' )	
		st_cut = graph.maxflow( st_net, source='s', target='t', capacity=E(st_net)$weight )	## find the min-cut
		if(!silent) cat( 'min_cut ends.\n' )	
		selected = sort( setdiff( V(st_net)[st_cut$partition1]$name, 's' ) ) 	## the selected nodes
		
		C = -sum( cf*(cf>=0) )	## a constant derived from the objective function
		##	returns the min-cut details
		st_cut_info = list( st_net = st_net, st_cut=st_cut, selected=selected, change_points=change_points, obj=-(st_cut$value + C) )
		return( st_cut_info )
	}




