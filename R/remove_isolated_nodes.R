remove_isolated <- function( input_net )
{
	output_net = induced.subgraph( input_net, which( degree( input_net ) > 0 ) )
	return( output_net )
}

