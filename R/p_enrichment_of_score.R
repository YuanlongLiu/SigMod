p_enrichment_score <- function(scored_net, module_genes_file, nperm=1E5)
{
	module_genes_frame = read.table( module_genes_file, header=TRUE, stringsAsFactors=FALSE, comment.char='' )
	module_genes = module_genes_frame$gene
	module_len = length( module_genes )
	scores = V(scored_net)$weight
	mod_score = sum( V(scored_net)[module_genes]$weight )
	perm_scores = sapply(1:nperm, function(v) sum(sample( scores, module_len )))
	p = (1 + sum(mod_score <= perm_scores)) / nperm
	cat('p_enrichment_score', '\n', p, file='./p_enrichment_score_result.txt', append=FALSE, sep='')	
}

