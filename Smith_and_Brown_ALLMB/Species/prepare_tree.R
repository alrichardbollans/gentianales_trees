library(here)
library(dplyr)
source(here('file_paths.R'))
source(here('helpful_phyl_methods.R'))

## Process
# - Extract gentianales clade
# - Standardise names in tree to WCVP and remove non-species tips and tips not in study families.
# - Remove duplicate tips
# - Add unrepresented species as polytomies

# A gentianales version of the smb tree, with names standardised
# The process for this will need improving e.g. by dropping genera and removing other families.
standardised_smb_tree = ape::read.tree(file.path('..','Genus','temp_outputs','standardised_smb_tree.tre'))

duplicated_tips = which(duplicated(standardised_smb_tree$tip.label))
# Check for repeated tips
print('duplicated tips will be removed, preserving a single instance')
print(length(duplicated_tips))
species_tree = ape::drop.tip(standardised_smb_tree,duplicated_tips)

# Remove non family/accepted tips
unaccepted_names_in_tree = c()
for(n in species_tree$tip.label){
  if(! replace_underscore_with_space_in_name(n) %in% all_trait_data$accepted_name){
    unaccepted_names_in_tree = c(unaccepted_names_in_tree,n)
  }
}
unaccepted_names_in_tree
accepted_tree = ape::drop.tip(species_tree, unaccepted_names_in_tree)

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
  i = 0
  for(name in names_not_in_tree){
    print(i)
    if(!(ape::is.ultrametric(tree))){
      # See: https://search.r-project.org/CRAN/refmans/phytools/html/add.species.to.genus.html
      tree=phytools::force.ultrametric(tree)
      }
    tree = phytools::add.species.to.genus(tree,name, where='root')
    i=i+1
  }
  # Needed to leave underscores in tree for phytools method, but now remove them for later processing steps
  tree = set_labels_on_tree_to_acc_name(tree)
  return(list('tree'=tree, 'names_not_in_tree'=names_not_in_tree, 'namesin_tree'=namesin_tree))
}

polytomy_output = add_species_as_polytomies(accepted_tree,all_trait_data$accepted_name)


cleaned_tree_with_polytomies = polytomy_output$tree
names_not_in_tree = polytomy_output$names_not_in_tree
namesin_tree = polytomy_output$namesin_tree

ape::write.tree(cleaned_tree_with_polytomies, prepared_smb_species_tree_path)

fileConn<-file(file.path('outputs','tree summaries','species not in tree.txt'))
writeLines(c(paste(names_not_in_tree, sep='')),fileConn)
close(fileConn)

fileConn<-file(file.path('outputs','tree summaries','species in tree.txt'))
writeLines(c(paste(namesin_tree, sep='')),fileConn)
close(fileConn)

