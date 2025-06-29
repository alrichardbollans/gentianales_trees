decompose_tree <- function(out_folder, tree_file){
  ground_tree = ape::read.tree(file.path(out_folder, tree_file))
  pvrd <- PVR::PVRdecomp(ground_tree)
  # get the vectors and their names into a data frame
  eigenVec.df <- data.frame(pvrd@Eigen$vectors,row.names=pvrd@phylo$tip.label)
  write.csv(eigenVec.df, file = file.path(out_folder,'all_eigenvectors.csv'))
  eigenValue.df <- data.frame(pvrd@Eigen$values)
  
  # As in Bachman, Steven Philip, et al. "Extinction risk predictions for the world's flowering plants to support their conservation." bioRxiv (2023): 2023-08.
  # use the Broken Stick method to select the eigenvectors
  # Discussed in Diniz-Filho, J. A. F., Sant'Ana, C. E. R. d., & Bini, L. M. (1998). An Eigenvector Method for Estimating Phylogenetic Inertia. Evolution, 52(5), 1247-1262. https://doi.org/10.2307/2411294
  # Original citation Jackson, Donald A. "Stopping rules in principal components analysis: a comparison of heuristical and statistical approaches." Ecology 74.8 (1993): 2204-2214.
  # The code calculates Broken Stick values cumulatively in reverse order (e.g., for j=1, it sums the last eigenvalue's expected variance). Sorting these values in decreasing order coincidentally aligns them with the correct order.
  eigenvalues <- eigenValue.df$pvrd.Eigen.values
  number_eigenvalues <- length(eigenvalues)
  broken_stick_model_df <- data.frame(j=seq(1:number_eigenvalues), p=0)
  # The first value of p is set to 1/n, which is the expected proportion of variance explained by the first eigenvector under the broken stick model.
  broken_stick_model_df$p[1] <- 1/number_eigenvalues 
  #The subsequent values of p are calculated iteratively by adding 1/(n+1-i) to the previous value, following the broken stick distribution.
  for (i in 2:number_eigenvalues)
  {
    broken_stick_model_df$p[i] = broken_stick_model_df$p[i-1] + (1/(number_eigenvalues + 1 - i))
  }
  #Finally, the p values are scaled to percentages by multiplying by 100 and dividing by n.
  broken_stick_model_df$p <- 100*broken_stick_model_df$p/number_eigenvalues
  
  # get a data frame with selected vectors
  eigenValue.df$percent.expected <- 100*eigenvalues/sum(eigenvalues) # Percent expected explained variance by the eigenvector according to eigenvalue
  eigenValue.df$broken.stick <- sort(broken_stick_model_df$p, decreasing = T) # Expected explained according to p value
  eigenValue.df$diff <- eigenValue.df$percent.expected - eigenValue.df$broken.stick
  selected_eigenValue.df <- eigenValue.df[eigenValue.df$diff > 0,] # Use those vectors where the variance explained according to eigenvalue is > according to p value
  num_selected = nrow(selected_eigenValue.df) # first X vectors selected
  
  explained_var= sum(selected_eigenValue.df$percent.expected) # X% of the variance
  
  param_df = data.frame(broken_stick_number=c(num_selected),explained_var=c(explained_var))
  
  
  write.csv(param_df,file.path(out_folder,"broken_stick_parameters.csv"))
}


decompose_tree('outputs/Species', 'Uphylomaker_species_tree.tre')
