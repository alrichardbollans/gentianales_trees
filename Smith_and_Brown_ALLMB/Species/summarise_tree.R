library(here)
source(here('prepared_data.R'))
source(here('helpful_phyl_methods.R'))

species_heatmap_plot(labelled_tree, data_with_tree_labels[c('label', 'APM.Activity')], 'Activity', 'species_actvity.svg', tipnames=TRUE )
species_heatmap_plot(labelled_tree, data_with_tree_labels[c('label', 'APM.Activity')], 'Activity', 'species_actvity.jpg', tipnames=TRUE,plotwidth = 11, plotheight=10,variable_offset=30,labelsize = 1.5 )
species_heatmap_plot(labelled_tree, data_with_tree_labels[c('label', 'APM.Activity')], 'Activity', 'species_actvity_no_tips.svg', tipnames=FALSE)

num_sp_in_original_tree = length(prepared_smb_species_tree$tip.label)
num_sp_in_labelled_tree = length(labelled_tree$tip.label)

fileConn<-file(file.path('outputs','tree summaries','summary.txt'))
writeLines(c(paste('num_sp_in_original_tree: ', num_sp_in_original_tree, sep=''),paste('num_sp_in_labelled_tree: ', num_sp_in_labelled_tree, sep = '')), 
           fileConn)
close(fileConn)
