gene_rank <- function( path, selected )
{
	counts = sort(table(unlist(path)))
	gene_rank = rank( counts[ selected ] )
	return( gene_rank )
}


