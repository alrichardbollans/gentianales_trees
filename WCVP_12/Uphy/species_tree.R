megatree = ape::read.tree(file.path('inputs', 'plant_megatree.tre'))

splist = read.csv(file.path('inputs', 'species_family_list.csv'))

## some sort of fix for hybrids
gen_list = read.csv(file.path('inputs', 'genus_family_list.csv'))

result = U.PhyloMaker::phylo.maker(sp.list = splist, tree = megatree, 
                          gen.list = gen_list, nodes.type = 1, scenario = 3)

gentianales_sp_tree = result$phylo

plot(gentianales_sp_tree)

### Genus tree?
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
number_of_tips = ape::Ntip(gentianales_sp_tree)
known_mrcas = c()
genera_to_add_as_polytomes = c()
for(genus in accepted_genus_list){
  
  associated_species = get_species_tips_for_genus(gentianales_sp_tree, genus)
  
  if (!is.null(associated_species) && length(associated_species)>0){
    if (length(associated_species)>1){
      # When there is at least one associated species, find the mrca and make that the genus node
      
      mrca = phytools::findMRCA(gentianales_sp_tree,tips = associated_species,type='node')
      
      if( mrca %in% known_mrcas){
        #  If mrca is shared with another genus, handle this later
        genera_to_add_as_polytomes = c(genera_to_add_as_polytomes, genus)
        
        
        
      } else{
        known_mrcas = c(known_mrcas, mrca)
        mrca_node_index = mrca - number_of_tips
        
        gentianales_sp_tree$node.label[mrca_node_index] = genus
      }
      
    }else{
      # If there is only one associated species, make that the genus tip
      
      gentianales_sp_tree$tip.label[gentianales_sp_tree$tip.label==associated_species[1]] = genus
    }
    
  }
}

associated_species = get_species_tips_for_genus(gentianales_sp_tree, 'Tabernaemontana')
mrca = phytools::findMRCA(gentianales_sp_tree,tips = associated_species,type='node')
g = ape::extract.clade(gentianales_sp_tree, mrca)
plot(g)
