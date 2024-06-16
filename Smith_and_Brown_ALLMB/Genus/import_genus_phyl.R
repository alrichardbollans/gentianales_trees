library(here)

## Process
# - Extract gentianales clade
# - Standardise names in tree to WCVP and remove synonyms (issues with e.g. Tabernaemontana/Nerium duplicates)
# - Where nodes are already a genus, make this into a tip and remove descendant species
# - For remaining genera, find the most recent common ancestors for their descendant species and rename this to the genus
# - If this MRCA is also the mrca for another genus, make this into a separate tip with zero length from the node
# - Drop any tips that aren't (accepted) genera
# - Remove duplicate tips

### Some smb methods

get_smb_genus_name_from_tree <-function(given_tree_name){
  new = stringr::str_split(given_tree_name,'_')[[1]][1]
  return(new)
}

resolve_tree_names_to_genus <- function(tree){
  tree$tip.label = as.character(lapply(tree$tip.label,get_smb_genus_name_from_tree))
  return(tree)
}

remove_duplicated_tips <- function(tree){
  # Check for repeated tips
  duplicated_tips = which(duplicated(tree$tip.label))
  if (length(duplicated_tips)>0){
    # Check for repeated tips
    print('duplicated tips will be removed, preserving a single instance')
    print(length(duplicated_tips))
    tree = ape::drop.tip(tree,duplicated_tips)
  }
  
  return(tree)
}

### Importing
## Get gentianales tree
smb_tree_ = ggtree::read.tree(file.path('..','v0.1','ALLMB.tre'))
gentianales_tree = ape::extract.clade(smb_tree_,c('Gentianales.rn.d8s.tre'))
ape::write.tree(gentianales_tree, file=file.path('..','SMB_ALLMB_Gentianales.tre'))
p = ggtree::ggtree(gentianales_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('..','SMB_ALLMB_Gentianales.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)
# 
# for (x in gentianales_tree$node.label){
#   if (grepl('Gal',x, fixed = TRUE)){
#     print(x)
#   }
# }

# Then reimport the tree with names standardised by wcvpy
# Matched with 'fuzzy' to avoid species being matched to genera
# There were quite a few genera synonyms in there! e.g. Labordia, Chazaliella etc..
standardised_smb_tree = ggtree::read.tree(file.path('temp_outputs','standardised_gentianales_smb_tree.tre'))
accepted_genus_list <- readLines(file.path('inputs','genus_list.txt'))

### Checks
testit::assert("Check name standardisation preserves tree tips", length(gentianales_tree$tip.label) == length(standardised_smb_tree$tip.label))
# This plot just confirms some of the branches creating a big polytomy from the main Rubiaceae branch
check = ape::keep.tip(gentianales_tree, c('Genipa_americana','Tarenna_leioloba','Tridentea_virescens','Raritebe_palicoureoides','Amphidasya_ambigua','Riqueuria_avenia', 'Polyura_geminata','Peripeplus_bracteosus', 'Neblinathamnus_argyreus', 'Neblinathamnus_brasiliensis'))
plot(check)

# This plot looks at the polyphyly of loganiaceae/gelsemiaceae
check = ape::keep.tip(gentianales_tree, c('Strychnos_cathayensis','Mostuea_surinamensis','Gelsemium_sempervirens','Geniostoma_huttonii'))
plot(check)

# Remove poor names
standardised_smb_tree = ape::drop.tip(standardised_smb_tree,c('NON_FAMILY_TIP'))

## Find nodes and tips of accepted genera
common_genera = intersect(accepted_genus_list,standardised_smb_tree$node.label)
common_genera
#### Add genus nodes as tips

# Function modified from http://blog.phytools.org/2012/11/adding-single-tip-to-tree.html
bind_node_tip<-function(tree,node_label){
  # Essentially makes a copy of the node as a tip (edge length is 0)
  node_index <-match(node_label, tree$node.label)+length(tree$tip.label)
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=node_label,
            edge.length=0,
            Nnode=1)
  class(tip)<-"phylo"
  obj<-ape::bind.tree(tree,tip,where=node_index)
  return(obj)
}

get_descendant_tips <- function(tree, node_label){
  node_index <-match(node_label, tree$node.label)+length(tree$tip.label)
  tips = geiger::tips(tree, node_index)
  tips = tips[! tips %in% c(node_label)] # Don't include node label in descendants
  # print(node_label)
  # print(tips)
  return(tips)
}

get_species_tips_for_genus <- function(tree, genus){
  # Gets all species tips that start with given genus
  tips = c()
  for(tip in tree$tip.label){
    if(startsWith(tip,paste(genus,"_",sep=""))){
      tips = c(tips,tip)
      
    }
  }
  
  tips = tips[! tips %in% c(genus)] # Don't include node label in descendants
  return(tips)
}

add_node_tip_and_remove_descendant <- function(tree, node_label){
  # For a given node label, makes it into a tip and then drops associated descendant species
  if(node_label %in% tree$node.label){
    tips_to_drop = get_species_tips_for_genus(tree,node_label)
    new_tree = ape::drop.tip(tree,c(node_label), trim.internal=FALSE, collapse.singles=FALSE) # if name is already a tip, drop it
    new_tree = bind_node_tip(new_tree,node_label)
    new_tree = ape::drop.tip(new_tree,tips_to_drop)
    return(new_tree)
  }else{
    print('WARNING: Following node no longer in tree. May have become a tip')
    print(node_label)
    
    return(tree)
  }
  
}

# For nodes that are known accepted genera, make them into tips and remove descendants
node_tip_tree = standardised_smb_tree
for (genus in common_genera) {
  node_tip_tree<-add_node_tip_and_remove_descendant(node_tip_tree,genus)
}

# Now check tips that are genera.
new_tip_genera = intersect(accepted_genus_list,node_tip_tree$tip.label)
new_tip_genera

testit::assert("Check extracting nodes preserves important nodes", length(common_genera) < length(new_tip_genera))
# plot(node_tip_tree)


# For remaining genera that weren't already nodes, find the most recent common ancestors for their species and
# give associated node the genus name
renamed_node_tip_tree = node_tip_tree
remaining_genera = setdiff(accepted_genus_list, common_genera) # genera left to process
testit::assert("Check all common_genera are now tips", length(setdiff(common_genera, renamed_node_tip_tree$tip.label)) ==0)
number_of_tips = ape::Ntip(renamed_node_tip_tree)
known_mrcas = c()
genera_to_add_as_polytomes = c()
for(genus in remaining_genera){

  associated_species = get_species_tips_for_genus(renamed_node_tip_tree, genus)
  
  if (!is.null(associated_species) && length(associated_species)>0){
    if (length(associated_species)>1){
      # When there is at least one associated species, find the mrca and make that the genus node
      
      mrca = phytools::findMRCA(renamed_node_tip_tree,tips = associated_species,type='node')
      
      if( mrca %in% known_mrcas){
        #  If mrca is shared with another genus, handle this later
        genera_to_add_as_polytomes = c(genera_to_add_as_polytomes, genus)
        
        
        
      } else{
        known_mrcas = c(known_mrcas, mrca)
        mrca_node_index = mrca - number_of_tips
        
        renamed_node_tip_tree$node.label[mrca_node_index] = genus
      }
      
    }else{
      # If there is only one associated species, make that the genus tip
      
      renamed_node_tip_tree$tip.label[renamed_node_tip_tree$tip.label==associated_species[1]] = genus
    }
    
  }
}

for(genus in genera_to_add_as_polytomes){
  associated_species = get_species_tips_for_genus(renamed_node_tip_tree, genus)
  mrca = phytools::findMRCA(renamed_node_tip_tree,tips = associated_species,type='node')
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=genus,
            edge.length=0,
            Nnode=1)
  class(tip)<-"phylo"
  
  renamed_node_tip_tree=ape::bind.tree(renamed_node_tip_tree,tip,where=mrca)
}


final_node_tip_tree = renamed_node_tip_tree
for (genus in remaining_genera) {
  final_node_tip_tree<-add_node_tip_and_remove_descendant(final_node_tip_tree,genus)
}

plot(final_node_tip_tree)

genera_not_in_tree = setdiff(common_genera,final_node_tip_tree$tip.label)
remaining_genera_not_in_tree = setdiff(remaining_genera,final_node_tip_tree$tip.label)

# Make all tip names into respective genus names and remove tips that aren't accepted genera.
# tree_with_names_resolved_to_genus = resolve_tree_names_to_genus(node_tip_tree)
tree_with_names_resolved_to_genus = final_node_tip_tree
new_common_genera = intersect(accepted_genus_list,final_node_tip_tree$tip.label)

non_genus_tips_to_remove = tree_with_names_resolved_to_genus$tip.label[! tree_with_names_resolved_to_genus$tip.label %in% c(new_common_genera)]
non_genus_tips_to_remove
genus_tree = ape::drop.tip(tree_with_names_resolved_to_genus, non_genus_tips_to_remove)
new_common_genera1 = intersect(accepted_genus_list,genus_tree$tip.label)
testit::assert("Check extracting nodes preserves important nodes", length(new_common_genera) == length(new_common_genera1))
plot(genus_tree)
p = ggtree::ggtree(genus_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs','SMB_Gentianales_Genus_tree.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)

deduplicated_genus_tree = remove_duplicated_tips(genus_tree)
# deduplicated_genus_tree = phytools::force.ultrametric(deduplicated_genus_tree, method='extend')
testit::assert("Check extracting nodes preserves important nodes", length(new_common_genera1) == length(intersect(accepted_genus_list,deduplicated_genus_tree$tip.label)))
p = ggtree::ggtree(deduplicated_genus_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs','SMB_Gentianales_Genus_tree_deduplicated.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)
ape::write.tree(deduplicated_genus_tree, file=file.path('outputs','final_SMB_Gentianales_genus_tree.tre'))



num_genera_in_tree = length(deduplicated_genus_tree$tip.label)

fileConn<-file(file.path('outputs','tree summaries','summary.txt'))
writeLines(c(paste('num_genera_in_final_tree: ', num_genera_in_tree, sep='')), 
           fileConn)
close(fileConn)

