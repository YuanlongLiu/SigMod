	##this function computes the selection path for a fixed \lambda, i.e., all distinct selection by varying \eta from eta_min to eta_max
	##for details please refer to Supplementary Notes Section 3.1 'Computing the selection path at any given \lambda value' of the SigMod paper
	##eta_max is set as very big (1E7) by default, which corresponds to the empty selection	
	##last modified: 03-Jan-2017

	##this is a more efficient version of the selection_path function
	##for details please refer to Section 2.3.1 'Computing the selection path at any given \lambda value' of the SigMod paper
	selection_path_fast <- function( net, lambda, eta_max=1E7, eta_min, selection_diff=1, selection_min=0, selection_max=Inf, silent )
	{
		if( eta_max <= eta_min ) { stop( '\neta_max should be greater than eta_min' ) }
		st_info = SigMod_core(net, lambda, eta=eta_min, silent) 
		selection = st_info$selected
		if( length(selection) == 0 ) return(selection)
		
		net_reduced = induced.subgraph( net, selection )
		# if( ecount( net_reduced ) < 10000 ) { cat('Loose network selected, quit\n'); return( selection ) }
		selections = selection_path( net_reduced, lambda, eta_max=eta_max, eta_min=eta_min, selection_diff, selection_min, selection_max, silent )
		return( selections )
	}
	
	