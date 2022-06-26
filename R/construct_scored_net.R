	##	this function constructs a scored gene network that will be used for searching gene modules
	##	Yuanlong LIU
	##	08-12-2017
	##	updated 05-01-2018

	construct_scored_net <- function( gene_network_data, interaction_indices, weight_index=NULL, gene_ps, gene_scores, genes2exclude = NULL )
	{

		if( !missing( gene_ps ) & !missing( gene_scores ) ) { stop('\nPlease provide either a gene p-value data or a gene score data, not both') }

		if( !missing( gene_ps ) )
		{
			cat('\nYou have provided a gene_ps data. Gene p-values will be first converted to scores by calling the \"p2score\" function\n')
			if( !all(c('gene', 'p') %in% names( gene_ps ) ) ) { stop('\nColumn names of gene_ps should contain \'gene\' and \'p\'') }
			gene_scores = p2score( gene_ps )
		}

		if( !missing( gene_scores ) )
		if( !all(c('gene', 'score') %in% names( gene_scores ) ) ) { stop('\nColumn names of gene_scores should contain \'gene\' and \'score\'') }

		gene_scores$gene = toupper( gene_scores$gene )
		gene_network_data[, interaction_indices[1]] = toupper( gene_network_data[, interaction_indices[1]] )
		gene_network_data[, interaction_indices[2]] = toupper( gene_network_data[, interaction_indices[2]] )
		
		if( !is.null( genes2exclude ) )
		{
			gene2exclude = toupper( genes2exclude )
			gene_scores = subset( gene_scores, !(gene %in% gene2exclude) )
		}
		
		
		ncols = ncol( gene_network_data )
		weight_real = which( colnames( gene_network_data ) == 'weight' )

		if( is.null(weight_index) )
		{
			if( length(weight_real) !=0 )
			{ 
				stop('\nThe name of the ', weight_real, 'th column of the gene_network_data is \'weight\'. By default (of the igraph package), this column will be converted to the interaction/edge weight. If you do not want this column to be converted to the interaction/edge weight, please change the column name as a value other than \'weight\'. If you want to convert this column as the interaction/edge weight, please assign the column index to the \'weight_index\' argument\n')
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
			n_is_numeric = any(!is.numeric(gene_network_data$weight) )
			is_na = any(is.na(gene_network_data$weight) )
			is_neg = any(gene_network_data$weight < 0 )
			
			if( n_is_numeric ) stop('Please make sure all interaction weights should be in a numeric format')
			if( is_na ) stop('Please make sure all interaction weights should not be a missing value')
			if( is_neg ) stop('Please make sure all interaction weights are equal to or greater than 0')
		}

		
		main_indices = c(interaction_indices, weight_index)
		gene_network_data = gene_network_data[, c( main_indices, setdiff( 1:ncols, main_indices ) )]

		raw_net = graph_from_data_frame( gene_network_data, directed = FALSE, vertices = NULL )
		if( !is_simple( raw_net ) ) 
		{
			## the following comments should be removed
			## cat('\nNote:\nThe network you provided is not a simple network\nA simple network is a network which do not contain loop and multiple edges\nWe convert it to a simple network using the igraph::simplify function\n')
			raw_net = simplify( raw_net )			
		}
		
		net = raw_net
		if( !'weight' %in% edge_attr_names( net ) )
		{
			## the following comments should be removed
			# cat('The network provided by the user does not have interaction/edge weight information. We set all interaction/edge weights as 1\n')
			E( net )$weight = 1
		}

		V(net)$weight = gene_scores[match(V(net)$name, gene_scores[, 'gene']), 'score']
		net = induced.subgraph(net, !is.na(V(net)$weight))
		net = induced.subgraph( net, degree(net) >=1 )
		if( 'p' %in% colnames(gene_scores) ) V(net)$p = gene_scores$p[match( V(net)$name, gene_scores$gene )]

		
		logfile = './constructed_scored_net.log'
		cat(paste(rep("#", 100), collapse=''), file=logfile, append=FALSE)
		cat('\n\nA scored-network has been successfully constructed\n', file=logfile, append=TRUE)
		cat('This scored-network consists', vcount(net), 'nodes and ', ecount(net), 'interactions\n', file=logfile, append=TRUE)
		cat('This scored-network is saved as an igraph object named \"net\" in \'constructed_scored_net.Rdata\'\n\n', file=logfile, append=TRUE)
		cat(paste(rep("#", 100), collapse=''), file=logfile, append=TRUE)
	
		save(net, file='constructed_scored_net.Rdata')
		return( net )
	}

	
	




