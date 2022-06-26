	## Yuanlong LIU
	## 08-12-2017


	construct_scored_net <- function( gene_network_data, interaction_indices, weight_index=NULL, gene_scores )
	{
		if(!is.data.frame( gene_network_data )) { stop('\ngene_network should be a data frame') }
		if(!is.data.frame( gene_scores ))  { stop('\ngene_scores should be a data frame') }
		if( !all(c('gene', 'score') %in% names( gene_scores ) ) ) { stop('\nColumn names of gene_scores should contain \'gene\' and \'score\'') }
		if( any( duplicated(colnames(gene_network_data)) ) ) { stop('\nColumn names of gene_scores should be different') }
		
		ncols = ncol( gene_network_data )
		weight_real = which( colnames( gene_network_data ) == 'weight' )

		if( is.null(weight_index) )
		{
			if( length(weight_real) !=0 )
			{ 
				stop('\nColumn ', weight_index, ' of the gene_network_data has a name \'weight\'\nBy default, this column will be converted to the interaction/edge weight\nIf you do not want to set this column as the interaction/edge weight, you need to change its name as a value other than \'weight\'. To set this column as the interaction/edge weight, you can assign its index to the \'weight_col\' parameter\n')
			}
		}
		
		if( !is.null(weight_index) )
		{
			if( length(weight_real) !=0 )
			if( weight_index != weight_real  )
			{ 
				stop('\nColumn ', weight_real, ' of the gene_network_data has a name \'weight\'\nBy default, the values of this column will be converted to the interaction/edge weight\nBut as you specified a different column as interaction/edge weight, you need to change the name of this column as a value other than \'weight\'')
			}
			names( gene_network_data )[ weight_index ] = 'weight'			
		}

		
		main_indices = c(interaction_indices, weight_index)
		gene_network_data = gene_network_data[, c( main_indices, setdiff( 1:ncols, main_indices ) )]

		raw_net = graph_from_data_frame( gene_network_data, directed = FALSE, vertices = NULL )
		if( !is_simple( raw_net ) ) 
		{
			## the following comments should be removed
			## cat('\nNote:\nThe network you provided is not a simple network\nA simple network is a network which do not contain loop and multiple edges\nWe convert it to a simple network using the igraph::simplify function\n')
			net = simplify( raw_net )			
		}
		
		if( !'weight' %in% vertex_attr_names( net ) )
		{
			## the following comments should be removed
			# cat('The network provided by the user does not have interaction/edge weight information. We set all interaction/edge weights as 1\n')
			E( net )$weight = 1
		}
		
		V(net)$weight = gene_scores[match(V(net)$name, gene_scores[, 'gene']), 'score']
		net = induced.subgraph(net, !is.na(V(net)$weight))
		net = induced.subgraph( net, degree(net) >=1 )

		logfile = './log.txt'
		cat(paste(rep("#", 100), collapse=''), file=logfile, append=TRUE)
		cat('\nA scored-network has been successfully constructed\n', file=logfile, append=TRUE)
		cat('This scored-network consists', vcount(net), 'nodes and ', ecount(net), 'interactions\n', file=logfile, append=TRUE)
		cat('This scored-network is saved as an igraph object named \"net\", and in \'./constructed_scored_net.Rdata\'\n', file=logfile, append=TRUE)
		cat(paste(rep("#", 100), collapse=''), file=logfile, append=TRUE)
	
		save(net, file='./constructed_scored_net.Rdata')
		return( net )
	}

	
	




