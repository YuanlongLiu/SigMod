	##	this function plots the histogram of the p-values of selected genes against the background of all gene p-values in the input network

	hist_p_values <- function( whole_ps, selected_ps, output='./hist_selected_gene_p_values.pdf' )
	{
		pdf( output, width=10, height=5 )
		par( mar=c( 4,5,2,1 ) )
		col_whole = 'lightgrey'
		col_selected = rgb(1,0,0,1)
		xlim = c(0, 1)
		p1 <- hist( whole_ps, breaks=seq( xlim[1], xlim[2], length=100+1), col=col_whole, xlim=xlim, main='Histogram of gene p-values', xlab='Gene p-value' ) 
		axis(1, at=0.05, labels=0.05, las=2, cex=.1)
		p2 <- hist( selected_ps, col=col_selected, xlim=xlim, breaks=seq( xlim[1], xlim[2], length=100+1), add=T )
		legend("topright", c("whole", "selected"), col=c(col_whole, col_selected), lwd=5, bty = "n")
		dev.off()
	}
	
	hist_scores <- function( whole_scores, selected_scores, output='./hist_selected_gene_p_values.pdf' )
	{
		pdf(output, width=10, height=5)
		par(mar=c(4,5,2,1))
		xlim = c(-9, 9)
		p1 <- hist(whole_scores, breaks=seq( xlim[1], xlim[2], length=100+1), col='lightgrey', xlim=xlim, main='Grey: histogram of p-values for all genes in the input network\nRed: histogram of p-values for selected genes', xlab='Gene p-value') 
		axis(1, at=0.05, labels=0.05, las=2, cex=.1)
		p2 <- hist(selected_scores, col=rgb(1,0,0,1), xlim=c(-9,9), breaks=seq( xlim[1], xlim[2], length=100+1), add=T)
		dev.off()
	}	
	selected_ps = 1-pnorm(V(res_info$opt_module[[1]])$weight)
	whole_ps = 1-pnorm(V(net)$weight)
	hist_p_values( whole_ps, selected_ps )
	
	

	

