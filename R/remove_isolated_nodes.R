	##	this function removes isolated nodes from a network
	##	updated 05-01-2018
	##	Yuanlong LIU
	remove_isolated <- function( input_net )
	{
		output_net = induced.subgraph( input_net, which( degree( input_net ) > 0 ) )
		return( output_net )
	}
