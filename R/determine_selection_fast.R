
	local_opt <- function( path, maxjump)
	{
		delta_select_size = diff(sapply( path, length ))
		maxjump_obs = max( delta_select_size )
		cat('Max size jump in this path:', maxjump_obs, '\n')
		flag = maxjump_obs > maxjump
		if(flag == 0)
		{
			opt_mod = list()
			return( opt_mod )
		}
		
		interest_indices = which( delta_select_size >= maxjump )
		interest_index = min(interest_indices)
		
		opt_mod = path[[ interest_index ]]
		return( opt_mod )	
	}
	max_comp <- function(net) { comps=decompose( net ); return( comps[[which.max( sapply( comps, vcount ) )]] ) }
	
	opt_mod <- function( net, paths, maxjump)
	{
		local_opts = lapply( paths, function(v) { local_opt(v, maxjump=maxjump) } )
		local_opts_len = sapply(local_opts, length)
		local_opts  = local_opts[which( local_opts_len ) > 0]

		local_opts = lapply( local_opts, function(v) { g=induced.subgraph(net, v); remove_isolated(g) } )
		return( local_opts ) ## this should be modified as the one has the maximal score
	}
	
	remove_isolated <- function( input_net )
	{
		output_net = induced.subgraph( input_net, which( degree( input_net ) > 0 ) )
		return( output_net )
	}

	
	