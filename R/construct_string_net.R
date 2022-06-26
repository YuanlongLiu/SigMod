	## Yuanlong LIU
	## 09-12-2016
	## 07-05-2017
	## 04-07-2017
	## 13-07-2017	

	remove_isolated <- function( input_net )
	{
		output_net = induced.subgraph( input_net, which( degree( input_net ) > 0 ) )
		return( output_net )
	}

	construct_string_net <- function( gene_ps, gene_scores, evidence2exclude = c('textmining_transferred', 'textmining'), genes2exclude = NULL, confidence=0.7 )
	{
		require( 'STRINGdb' )
		if( !missing( gene_ps ) & !missing( gene_scores ) ) { stop('\nPlease provide either a gene p-value file or a gene score file, not both') }
		if( !missing( gene_ps ) )
		{
			cat('\nYou have provided a gene_ps.\nGene p-values will be first converted to scores by calling the p2score function\n')
			
			gene_scores = p2score( gene_ps )
		}
		if( !missing( gene_scores ) ) 	gene_scores = gene_scores
		if( (confidence > 1) | (confidence <= 0) ) 	{ stop('\nconfidence should fall in (0, 1]') }
		
		cat('\n\nConstructing the network takes around 2 minutes (depends on your Internet connection)\n')
		
		if( !is.null( genes2exclude ) )
		{
			if( genes2exclude == 'OR_genes' )
			{
				if( !exists('OR_gene') ) stop( '\nYou have specified to exclude the genes in the OR gene family for downstream analysis. Please load ./OR_genes.Rdata into your workspace\n' )
				cat('\n\ngenes in the OR gene family are excluded for the downstream analysis\n\n')
				gene2exclude = OR_gene
			}
			
			if( genes2exclude != 'OR_genes' )
			{
				tmp = read.table( genes2exclude, stringsAsFactors = FALSE, header=TRUE )
				gene2exclude = toupper( unlist(tmp) )
			}

			gene_scores = subset( gene_scores, !toupper(gene) %in% gene2exclude )

		}
		
		## STRINGdb is an R function available which allows to access STRING online
		string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
		gene_mapped = string_db$map( gene_scores[, 'gene', drop=FALSE], "gene", removeUnmappedRows = TRUE, quiet=FALSE )
		cat('\nThis resulted', dim(gene_mapped)[1], 'genes mapped to STRING\n' )
		
		## convert gene names to upper case
		gene_mapped$gene = toupper( gene_mapped$gene )
		net = string_db$get_subnetwork(gene_mapped$STRING_id)

		e_attrs = as_data_frame(net, what = c("edges")) ## from the igraph package
		cat('\nSources of interaction evidence in the original STRING network:\n')
		evidences = setdiff( colnames(e_attrs), c('combined_score', 'from', 'to' ) )
		print(evidences)	
		
		evidence_keep = setdiff(evidences, evidence2exclude)
		cat('\nSources of interaction evidence kept for further analysis:\n')
		print(evidence_keep)
		
		combined_score_new = apply(e_attrs[, evidence_keep]/1000, 1, function(v) { 1-prod(1-v) } )
		E(net)$combined_score_new = combined_score_new
		E(net)$weight = E(net)$combined_score_new

		cor_combined_score_old_new = cor(combined_score_new, e_attrs[,'combined_score']/1000)
		cat('\nCorrelation between the original combined edge scores and the new edge scores \ncomputed by excluding evidences from', evidence2exclude, ':\n', cor_combined_score_old_new, '\n' )

		V(net)$ens = V(net)$name
		V(net)$name = gene_mapped$gene[match( V(net)$ens, gene_mapped$STRING_id )]
		V(net)$weight = gene_scores$score[match( V(net)$name, toupper(gene_scores$gene) )]
		
		net_hic = delete_edges(net, which(E(net)$weight < confidence))
		net_hic = remove_isolated( simplify( net_hic ) )
		
		# cat('\nNumber of STRING edges in the mapped network is', ecount(net))
		percent = format(100*ecount(net_hic)/ecount(net), scientific=FALSE, digits=2)
		# cat('\nNumber of STRING edges in the mapped network with confidence greater than', confidence, 'is', ecount(net_hic), '(', percent, '% )', '\n')
		cat('\n\nSTRING net successfully constructed\n')
		
		string_net = net_hic
		save(string_net, file='./constructed_scored_string_net.Rdata')
		return( string_net )
	}

	
	




