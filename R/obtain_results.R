	##	this function retrives the final results from the computed selection information
	##	updated 05-01-2018
	##	Yuanlong LIU
	
	obtain_results <- function( computation_details, net )
	{
		
		res_info = computation_details
		flag_next = length(res_info$opt_module) > 1
		flag_p = 'p' %in% vertex_attr_names(net) ## whether the network has a p-value attributes

		
		
		# Retrieve genes in the selected module
		selected_module = res_info$opt_module[[1]]
		selected_module_mat = data.frame(as.matrix(selected_module[]))
		write.table(selected_module_mat, file = "selected_module_interactions.tab", col.names=NA, sep="\t")
		
		selected_genes = V(selected_module)$name
		## order the genes by gene score:
		selected_genes = selected_genes[order( V(selected_module)$weight, decreasing=TRUE ) ]

		if(flag_p) selected_genes_info = data.frame(gene = selected_genes, p = V(net)[selected_genes]$p, score = V(net)[selected_genes]$weight)
		if(!flag_p) selected_genes_info = data.frame(gene = selected_genes, score = V(net)[selected_genes]$weight)

		write.table(selected_genes_info, file='./selected_genes.tab', quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
		
		# plot the selected module
		main = paste('Number of nodes:', vcount(selected_module), '\n', 'Number of interactions:', ecount(selected_module))
		pdf('selected_module.pdf')
		par(mar=c(0,0,2,0))
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label=NA, main=main)
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label.cex=0.5, main=main)
		dev.off()
		
		if(flag_p) hist_selected_p_values( whole_ps=V(net)$p, selected_ps=V(selected_module)$p, output='hist_selected_genes_p_values.pdf' )
		
		if( !flag_next ) return()

		#################	if there are more than one module    #################	

		selected_module = res_info$opt_module[[2]]
		selected_module_mat = data.frame(as.matrix(selected_module[]))
		write.table(selected_module_mat, file = "selected_module_next_interactions.tab", col.names=NA, sep="\t")

		
		selected_genes = V(selected_module)$name
		## order the genes by gene score:
		selected_genes = selected_genes[order( V(selected_module)$weight, decreasing=TRUE ) ]

		if(flag_p) selected_genes_info = data.frame(gene = selected_genes, p = V(net)[selected_genes]$p, score = V(net)[selected_genes]$weight)
		if(!flag_p) selected_genes_info = data.frame(gene = selected_genes, score = V(net)[selected_genes]$weight)
	
		write.table(selected_genes_info, file='./selected_genes_next.tab', quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
		
		# plot the selected module
		main = paste('Number of nodes:', vcount(selected_module), '\n', 'Number of interactions:', ecount(selected_module))
		pdf('./selected_module_next.pdf')
		par(mar=c(0,0,2,0))
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label=NA, main=main)
		set.seed(1)
		plot( selected_module, vertex.size = 5, vertex.label.cex=0.5, main=main)
		dev.off()
		
		if(flag_p) hist_selected_p_values( whole_ps=V(net)$p, selected_ps=V(selected_module)$p, output='hist_selected_genes_next_p_values.pdf' )
		cat('\nResults saved under:\n', getwd(), '\n')
	}
	