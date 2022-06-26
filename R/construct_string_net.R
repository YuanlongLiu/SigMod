	## Yuanlong LIU
	## 09-12-2016
	## 07-05-2017
	## 04-07-2017
	## 13-07-2017	


	construct_string_net <- function( gene_ps, gene_scores, evidence2exclude = c('textmining_transferred', 'textmining'), genes2exclude = NULL, confidence=0.7 )
	{
		if( (confidence > 1) | (confidence < 0) ) 	{ stop('\nconfidence should fall in [0, 1]') }
		if( !missing( gene_ps ) & !missing( gene_scores ) ) { stop('\nPlease provide either a gene p-value data or a gene score data, not both') }

		require( 'STRINGdb' )
		if(!"STRINGdb" %in% (.packages())) {
			stop("\nNeed to install STRINGdb package\nLink to the STRINGdb package: https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html")
		}
		
		if( !missing( gene_ps ) )
		{
			cat('\nYou have provided a gene_ps data. Gene p-values will be first converted to scores by calling the \"p2score\" function\n')
			gene_scores = p2score( gene_ps )
		}

		if( !missing( gene_scores ) ) 	gene_scores = gene_scores
		gene_scores$gene = toupper( gene_scores$gene )

		cat('\nAll gene names have been converted to uppercases\n')
		cat('\n\nConstructing the network takes around 2 minutes (depends on your Internet connection)\n')
		
		if( !is.null( genes2exclude ) )
		{
			if( genes2exclude == 'OR_genes' )
			{
				cat('\n\ngenes in the OR gene family are excluded for the downstream analysis\n\n')
				OR_gene = OR_gene_family_info()
				gene2exclude = toupper( OR_gene )
			}
			
			if( genes2exclude != 'OR_genes' )
			{
				gene2exclude = toupper( gene2exclude )
			}

			gene_scores = subset( gene_scores, !(gene %in% gene2exclude) )

		}
		
		## STRINGdb is an R function available which allows to access STRING online
		string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" )
		gene_mapped = string_db$map( gene_scores[, 'gene', drop=FALSE], "gene", removeUnmappedRows = TRUE, quiet=FALSE )
		cat('\nThis resulted', dim(gene_mapped)[1], 'genes mapped to STRING\n' )
		
		net = string_db$get_subnetwork(gene_mapped$STRING_id)

		e_attrs = as_data_frame(net, what = c("edges")) ## from the igraph package
		cat('\nSources of interaction evidence in the original STRING network:\n')
		evidences = setdiff( colnames(e_attrs), c('combined_score', 'from', 'to' ) )
		print(evidences)	
		
		evidence_keep = setdiff(evidences, evidence2exclude)
		cat('\nSources of interaction evidence kept for further analysis:\n')
		print(evidence_keep)
		
		combined_score = e_attrs[,'combined_score']/1000	##	the original combined score
		if(length( evidence_keep ) != length( evidences )) combined_score_new = apply(e_attrs[, evidence_keep]/1000, 1, function(v) { 1-prod(1-v) } )
		if(length( evidence_keep ) == length( evidences )) combined_score_new = combined_score

		E(net)$combined_score_new = combined_score_new
		E(net)$weight = E(net)$combined_score_new

		cor_combined_score_old_new = cor(combined_score_new, combined_score)
		cat('\nCorrelation between the original combined edge scores and the new edge scores computed by excluding evidences from', evidence2exclude, ':\n', cor_combined_score_old_new, '\n' )

		V(net)$ens = V(net)$name
		V(net)$name = gene_mapped$gene[match( V(net)$ens, gene_mapped$STRING_id )]
		V(net)$weight = gene_scores$score[match( V(net)$name, gene_scores$gene )]
		if( 'p' %in% colnames(gene_scores) ) V(net)$p = gene_scores$p[match( V(net)$name, gene_scores$gene )]
		
		net_hic = delete_edges(net, which(E(net)$weight < confidence))
		net_hic = remove_isolated( simplify( net_hic ) )
		
		# cat('\nNumber of STRING edges in the mapped network is', ecount(net))
		percent = format( 100*ecount(net_hic)/ecount(net), scientific=FALSE, digits=2 )
		# cat('\nNumber of STRING edges in the mapped network with confidence greater than', confidence, 'is', ecount(net_hic), '(', percent, '% )', '\n')
		cat('\n\nSTRING net successfully constructed\n')
		
		string_net = net_hic

		logfile = './constructed_string_net.log'
		cat(paste(rep("#", 100), collapse=''), file=logfile, append=FALSE)
		cat('\n\nA scored-network has been successfully constructed\n', file=logfile, append=TRUE)
		cat('This scored-network consists', vcount(string_net), 'nodes and ', ecount(string_net), 'interactions\n', file=logfile, append=TRUE)
		cat('This scored-network is saved as an igraph object named \"string_net\" in \'constructed_scored_string_net.Rdata\'\n\n', file=logfile, append=TRUE)
		cat(paste(rep("#", 100), collapse=''), file=logfile, append=TRUE)
		
		save(string_net, file='./constructed_scored_string_net.Rdata')
		return( string_net )
	}






