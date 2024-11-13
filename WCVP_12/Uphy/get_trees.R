# devtools::install_github("jinyizju/U.PhyloMaker") #https://github.com/jinyizju/U.PhyloMaker

megatree = ape::read.tree(file.path('inputs', 'plant_megatree.tre'))# from https://github.com/megatrees/plant_20221117/tree/main

splist = read.csv(file.path('inputs', 'species_family_list.csv'))

## TODO: Currently hybrids have been removed as I tried a few fixes that didnt work. some sort of fix for hybrids
gen_list = read.csv(file.path('inputs', 'genus_family_list.csv'))

result = U.PhyloMaker::phylo.maker(sp.list = splist, tree = megatree, 
                          gen.list = gen_list, nodes.type = 1, scenario = 3)

gentianales_sp_tree = result$phylo


ape::write.tree(gentianales_sp_tree, file.path('outputs', 'Species',
                                                        'Uphylomaker_species_tree.tre'))

num_sp_in_tree = length(gentianales_sp_tree$tip.label)
fileConn<-file(file.path('outputs','Species','summary.txt'))
writeLines(c(paste('num_sp_in_final_tree: ', num_sp_in_tree, sep='')), 
           fileConn)
close(fileConn)

p = ggtree::ggtree(gentianales_sp_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs','Species',
                               'Uphylomaker_species_tree.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)

# plot(gentianales_sp_tree)

### Genus tree
# Genus tips are found from the MRCA of species in the tree. 
# This is similar to the build.nodes.2 option of U.phylomaker 

accepted_genus_list = gen_list$genus
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

aden = get_species_tips_for_genus(gentianales_sp_tree,'Adenorandia')


port = get_species_tips_for_genus(gentianales_sp_tree,'Portlandia')
port_node = phytools::findMRCA(gentianales_sp_tree,tips = port,type='node')
port_tree = ape::extract.clade(gentianales_sp_tree, port_node)
plot(port_tree)

leuco = get_species_tips_for_genus(gentianales_sp_tree,'Leucolophus')

port_node = phytools::findMRCA(gentianales_sp_tree,tips = c(aden,leuco, port),type='node')
port_tree = ape::extract.clade(gentianales_sp_tree, port_node)
port_sp_tree = ape::keep.tip(gentianales_sp_tree, c(aden,leuco, port))
plot(port_sp_tree)

## Find information on genera in the tree
# If there is only one associated species tip, make that the genus tip
number_of_tips = ape::Ntip(gentianales_sp_tree)
known_mrcas = c()
genera_to_add_as_polytomes = c()
remaining_genera = c()
genera_with_single_tips = c()
final_node_tip_tree = gentianales_sp_tree
for(genus in accepted_genus_list){
  
  associated_species = get_species_tips_for_genus(final_node_tip_tree, genus)
  
  if (!is.null(associated_species) && length(associated_species)>0){
    if (length(associated_species)>1){
      # When there is at least one associated species, find the mrca and make that the genus node
      
      mrca = phytools::findMRCA(final_node_tip_tree,tips = associated_species,type='node')
      
      if( mrca %in% known_mrcas){
        #  If mrca is shared with another genus, handle this later
        genera_to_add_as_polytomes = c(genera_to_add_as_polytomes, genus)
        
        
        
      } else{
        known_mrcas = c(known_mrcas, mrca)
        mrca_node_index = mrca - number_of_tips
        
        final_node_tip_tree$node.label[mrca_node_index] = genus
        remaining_genera = c(remaining_genera, genus)
      }
      
    }else{
      # If there is only one associated species, make that the genus tip
      
      final_node_tip_tree$tip.label[final_node_tip_tree$tip.label==associated_species[1]] = genus
      genera_with_single_tips = c(genera_with_single_tips, genus)
      }
    
  }
}



# Where a genus shares an mrca with another genus, it is added as a polytomy at that mrca node
for(genus in genera_to_add_as_polytomes){
  associated_species = get_species_tips_for_genus(final_node_tip_tree, genus)
  mrca = phytools::findMRCA(final_node_tip_tree,tips = associated_species,type='node')
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=genus,
            edge.length=0,
            Nnode=1)
  class(tip)<-"phylo"
  
  final_node_tip_tree=ape::bind.tree(final_node_tip_tree,tip,where=mrca)
}

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

# Now add the remaining genera as tips at the respective mrca nodes
for (genus in remaining_genera) {
  final_node_tip_tree<-add_node_tip_and_remove_descendant(final_node_tip_tree,genus)
}

# Remove species tips from genera added as polytomies
for (genus in genera_to_add_as_polytomes){
  tips_to_drop = get_species_tips_for_genus(final_node_tip_tree,genus)
  final_node_tip_tree = ape::drop.tip(final_node_tip_tree,tips_to_drop)
}

genus_tree = final_node_tip_tree
# plot(genus_tree)

non_genus_tips = setdiff(genus_tree$tip.label, accepted_genus_list)
testit::assert("Check no erroneous tips", length(non_genus_tips) == 0)
genera_not_in_tree = setdiff(accepted_genus_list, genus_tree$tip.label)
genera_in_tree = intersect(accepted_genus_list,genus_tree$tip.label)

duplicated_tips = which(duplicated(genus_tree$tip.label))
testit::assert("Check no duplicates", length(duplicated_tips) == 0)

p = ggtree::ggtree(genus_tree,layout="circular") +
  ggtree::geom_tiplab2(size=2, show.legend=FALSE)
ggplot2::ggsave(file=file.path('outputs', 'Genus','Uphylomaker_genus_tree.jpg'),width=20, height=16,
                dpi = 300, limitsize=FALSE)
ape::write.tree(genus_tree, file=file.path('outputs', 'Genus','Uphylomaker_genus_tree.tre'))



num_genera_in_tree = length(genus_tree$tip.label)

fileConn<-file(file.path('outputs', 'Genus','Uphylomaker_genus_tree_summary.txt'))
writeLines(c(paste('num_genera_in_final_tree: ', num_genera_in_tree, sep='')), 
           fileConn)
close(fileConn)

