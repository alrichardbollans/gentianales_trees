library(here)
library(dplyr)

## Process
# - Extract gentianales clade
# - Standardise names in tree to WCVP and remove synonyms (issues with e.g. Tabernaemontana/Nerium duplicates)
# - Remove duplicate tips and non accepted name tips
# - Add unrepresented species as polytomies. Note that not all species are succesfully matched to a genus
# A gentianales version of the smb tree, with names standardised

standardised_smb_tree = ape::read.tree(file.path('..','Genus','temp_outputs','standardised_gentianales_smb_tree.tre'))
standardised_smb_tree = ape::drop.tip(standardised_smb_tree,c('NON_FAMILY_TIP'))

duplicated_tips = which(duplicated(standardised_smb_tree$tip.label))
# Check for repeated tips
print('duplicated tips will be removed, preserving a single instance')
print(length(duplicated_tips))
species_tree = ape::drop.tip(standardised_smb_tree,duplicated_tips)

accepted_species_list <- readLines(file.path('inputs','species_list.txt'))

replace_underscore_with_space_in_name<- function(given_tree_name){
  new = gsub("_", " ",given_tree_name)
  return(new)
}
# Remove non family/accepted tips
unaccepted_names_in_tree = c()
for(n in species_tree$tip.label){
  if(! replace_underscore_with_space_in_name(n) %in% accepted_species_list){
    unaccepted_names_in_tree = c(unaccepted_names_in_tree,n)
  }
}
unaccepted_names_in_tree
accepted_tree = ape::drop.tip(species_tree, unaccepted_names_in_tree)

replace_space_with_underscore_in_name<- function(given_name){
  new = gsub(" ", "_",given_name)
  return(new)
}

set_labels_on_tree_to_acc_name<- function(tree){
  tree$tip.label = as.character(lapply(tree$tip.label,replace_underscore_with_space_in_name))
  return(tree)
}

add_species_as_polytomies <- function(tree, names) {
  
  names_not_in_tree = c()
  namesin_tree = c()
  for(name in names){
    if(!replace_space_with_underscore_in_name(name) %in% tree$tip.label){
      names_not_in_tree = c(names_not_in_tree,name)
    }
    if(replace_space_with_underscore_in_name(name) %in% tree$tip.label){
      namesin_tree = c(namesin_tree,name)
    }
  }
  
  testit::assert("Check tree names are all in the given accepted names", length(tree$tip.label) == length(namesin_tree))
  # i = 0
  for(name in names_not_in_tree){
    # print(i)
    if(!(ape::is.ultrametric(tree))){
      # See: https://search.r-project.org/CRAN/refmans/phytools/html/add.species.to.genus.html
      tree=phytools::force.ultrametric(tree)
      }
    tree = phytools::add.species.to.genus(tree,name, where='root')
    # i=i+1
  }
  # Needed to leave underscores in tree for phytools method, but now remove them for later processing steps
  tree = set_labels_on_tree_to_acc_name(tree)
  return(list('tree'=tree, 'names_not_in_tree'=names_not_in_tree, 'namesin_tree'=namesin_tree))
}

polytomy_output = add_species_as_polytomies(accepted_tree,accepted_species_list) # Takes about twenty minutes


cleaned_tree_with_polytomies = polytomy_output$tree
names_not_in_tree = polytomy_output$names_not_in_tree
namesin_tree = polytomy_output$namesin_tree

ape::write.tree(cleaned_tree_with_polytomies, file.path('outputs','final_SMB_Gentianales_species_tree.tre'))

fileConn<-file(file.path('outputs','tree summaries','species not in tree.txt'))
writeLines(c(paste(names_not_in_tree, sep='')),fileConn)
close(fileConn)

fileConn<-file(file.path('outputs','tree summaries','species in tree.txt'))
writeLines(c(paste(namesin_tree, sep='')),fileConn)
close(fileConn)

num_sp_in_tree = length(cleaned_tree_with_polytomies$tip.label)
fileConn<-file(file.path('outputs','tree summaries','summary.txt'))
writeLines(c(paste('num_sp_in_final_tree: ', num_sp_in_tree, sep='')), 
           fileConn)
close(fileConn)

p = ggtree::ggtree(cleaned_tree_with_polytomies,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs','SMB_Gentianales_Species_tree.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)
