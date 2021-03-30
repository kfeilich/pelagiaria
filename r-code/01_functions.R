# 01_functions.R

# Author: Kara Feilich
# Date modified: 26 November 2018
# This script contains custom functions used for running the analysis, listed
# in alphabetical order by function name.

#' For a paired tree and rank envelope output, this function returns whether the
#' datum at each relative time point is lower (=1) than the rank envelope, or 
#' whether it isn't (=0)
#' @param tree A phylo object
#' @param rank.env.output  An object of class "envelope_test"
#' 
#' @return A dataframe showing whether each time point is lower than envelope
checkIsLower <- function(tree, rank.env.output){
  # First, put together all of the relevant variables in a data.frame for ease
  # of analysis
  time.env.output <- data.frame("relative.time" = rank.env.output$r, 
                                "data.curve" = rank.env.output$data_curve,
                                "lower" = rank.env.output$lower,
                                "upper" = rank.env.output$upper)
  # Calculate the height of the tree (root age)
  tree.height <- max(phytools::nodeHeights(tree))
  # Calculate the actual time from each relative time
  time.env.output$actual.time <- (1-time.env.output$r) * tree.height
  # Initialize each check of whether the data are lower than the envelope
  # with 0 -- indicating that the data are not lower than the envelope
  time.env.output$is.it.lower <- 0
  # For the timepoints that actually ARE lower, make is.it.lower = 1
  time.env.output$is.it.lower[which(
    time.env.output$data.curve < time.env.output$lower)] <-1
  # Outputs the whole dataframe
  time.env.output$row.index <- as.numeric(seq(1:nrow(time.env.output)))
  return(time.env.output)
}

#' For a paired tree and rank envelope output, this function returns whether the
#' datum at each relative time point is higher (=1) than the rank envelope, or 
#' whether it isn't (=0)
#' @param tree A phylo object
#' @param rank.env.output  An object of class "envelope_test"
#' 
#' @return A dataframe showing whether each time point is lower than envelope
checkIsHigher <- function(tree, rank.env.output){
  # First, put together all of the relevant variables in a data.frame for ease
  # of analysis
  time.env.output <- data.frame("relative.time" = rank.env.output$r, 
                                "data.curve" = rank.env.output$data_curve,
                                "lower" = rank.env.output$lower,
                                "upper" = rank.env.output$upper)
  # Calculate the height of the tree (root age)
  tree.height <- max(phytools::nodeHeights(tree))
  # Calculate the actual time from each relative time
  time.env.output$actual.time <- (1-time.env.output$r) * tree.height
  # Initialize each check of whether the data are lower than the envelope
  # with 0 -- indicating that the data are not lower than the envelope
  time.env.output$is.it.higher <- 0
  # For the timepoints that actually ARE lower, make is.it.lower = 1
  time.env.output$is.it.higher[which(
    time.env.output$data.curve > time.env.output$upper)] <-1
  # Outputs the whole dataframe
  time.env.output$row.index <- as.numeric(seq(1:nrow(time.env.output)))
  return(time.env.output)
}


#' Matches a phylogenetic tree tips and a trait dataframe 
#' 
#' clean.tree() cleans the names in order to match a single Pelagia tree with
#' the aspect.ratios data set
#' @param tree A phylo object
#' @param traits The dataframe of a single continuous trait variable
#' 
#' @returns a cleaned phylo object
clean.tree <- function(tree, traits) {
  cleaned.tree <- tree
  # checks data in tree, traits
  matches <- geiger::name.check (cleaned.tree, traits) 
    # drops non-matching taxa
  cleaned.tree <- ape::drop.tip (cleaned.tree, c(matches$tree_not_data))  
  return(cleaned.tree)
}  

#' calls clean.tree() on list of trees
#' 
#' clean.trees() is a wrapper to call clean.tree on a list of phylo objects, 
#' and returns a list of outputs.
#' @param trees A list of phylo objects, in this case, the posterior trees
#' @param traits The dataframe of a single continuous trait variable
#' 
#' @return A list of cleaned phylo objects
clean.trees <- function(trees, traits) {
  cleaned.trees <- sapply(X = trees, FUN = clean.tree, traits, simplify = FALSE)
  return(cleaned.trees)
}  

#' For the output of checkIsLower, return runs of 2 or more is.it.lower == 1
#' 
#' This function shows runs of 2 or more timepoints that were lower than the 
#' rank envelope.
#' @param is.it.lower.result A dataframe output by checkIsLower()
#' 
#' @return A dataframe with run start and end points for sig. lower runs
findLowerRuns <- function(is.it.lower.result){
  is.lower.runs <- with(rle(is.it.lower.result$is.it.lower), {
    ok <- values == 1 & lengths > 2
    ends <- cumsum(lengths)
    starts <- ends - lengths + 1
    data.frame(starts, ends)[ok, ]
  })
  is.lower.runs$starts <- is.it.lower.result[is.lower.runs$starts, "actual.time"]
  is.lower.runs$ends <- is.it.lower.result[is.lower.runs$ends, "actual.time"]
  return(is.lower.runs)
}

#' For the output of checkIsHigher, return runs of 2 or more is.it.higher == 1
#' 
#' This function shows runs of 2 or more timepoints that were lower than the 
#' rank envelope.
#' @param is.it.higher.result A dataframe output by checkIsLower()
#' 
#' @return A dataframe with run start and end points for sig. lower runs
findHigherRuns <- function(is.it.higher.result){
  is.higher.runs <- with(rle(is.it.higher.result$is.it.higher), {
    ok <- values == 1 & lengths > 2
    ends <- cumsum(lengths)
    starts <- ends - lengths + 1
    data.frame(starts, ends)[ok, ]
  })
  is.higher.runs$starts <- is.it.higher.result[is.higher.runs$starts, "actual.time"]
  is.higher.runs$ends <- is.it.higher.result[is.higher.runs$ends, "actual.time"]
  return(is.higher.runs)
}

#' Returns the value of the MDI statistic from the output of dtt1().
#' 
#' getMDI() returns the p-value from the output of the rank envelope test
#' from Murrell 2018
#' 
#' @param dtt1.output An object of class "envelope_test"
#' @return a number (the MDI statistic)
getMDI <- function(dtt1.output){
  mdi <- dtt1.output$MDI
  return(mdi)
}

#' Returns a MDI statistic p-value from the output of dtt1().
#' 
#' getMDIPvalue() returns the p-value from the output of the rank envelope test
#' from Murrell 2018
#' 
#' @param dtt1.output An object of class "envelope_test"
#' @return a number (the p value)
getMDIPvalue <- function(dtt1.output){
  p.value <- dtt1.output$MDIpVal
  return(p.value)
}

#' Returns p-value interval bounds from the output of a rank envelope test
#' 
#' getPinterval() returns the p-value interval from the output of the rank 
#' envelope test from Murrell 2018
#' 
#' @param rank.env.output An object of class "envelope_test"
#' @return a dataframe with upper and lower bounds for the p value
getPinterval <-function(rank.env.output){
  p.int <- rank.env.output$p_interval
  p.int.df <- data.frame("Lower" = p.int[1], "Upper" = p.int[2])
  return(p.int.df)
}

#' Returns a p value from the output of a rank envelope test
#' 
#' getPvalue() returns the p-value from the output of the rank envelope test
#' from Murrell 2018
#' 
#' @param rank.env.output An object of class "envelope_test"
#' @return a number (the p value)
getPvalue <- function(rank.env.output){
  p.value <- rank.env.output$p
  return(p.value)
}

#'  This is a custom ggplot theme to meet journal style guidelines
theme_procB <- function() {
  theme_bw(base_size = 12, base_family = "Times New Roman")  %+replace% 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


#' Returns all the results from running the modified dtt with rank envelope tests
#' 
#' Calls on the functions above to return a list of all of the outputs from the
#' DTT, MDI, and rank envelope tests for a given continuous trait and a list of
#' trees
#' 
#' @param trait_vector a named vector of trait values for the tips
#' @param list_of_trees a multiphylo object
#' @param n_sim the number of simulations for the rank envelope test
#' @return a list of all of the results listed above, as well as any binned times
#' that were significantly lower than the rank envelope, if there were any
run_modified_dtt <- function(trait_vector, list_of_trees, n_sim = 10000, 
                             test_type = c("less", "greater", "two.sided")) { 
  
  trait_vector_clean <- na.omit(trait_vector)
  cleaned_trees <- clean.trees(list_of_trees, trait_vector_clean)
  
  dtt.results.all <- purrr::map(cleaned_trees, function(x) 
    dtt1(x, trait_vector_clean, 
         plot=FALSE, 
         nsim = n_sim, 
         calculateMDIp=TRUE))
  
  mdi.values <- purrr::map_dbl(dtt.results.all, getMDI)
  mdi.p.values <- purrr::map_dbl(dtt.results.all, getMDIPvalue)
  mdi.output <- data.frame("MDI" = mdi.values, 
                           "MDI.pval" = mdi.p.values)
  
  rank.envelope.results <- purrr::map(dtt.results.all, rank_env_dtt, 
                                      Plot = F, test = test_type)
  
  # Find significant times #######################################################

  # First, we need maximum root age of all of the trees
  unnamed.trees <- cleaned_trees  
  names(unnamed.trees) <- NULL # Need to clear off the names or result is wonky
  root.ages <- purrr::map_dbl(unnamed.trees, 
                              function(x) max(phytools::nodeHeights(x)))
  maximum.root.age <- max(root.ages)
  
  # Call checkIsLower and/or checkIsHigher on all tree-rankenv pairs. 
  # purrr::map2 iterates the function over each x, y pair of 
  # (cleaned.trees, rank.envelope.results). Then, summarize significant time 
  # points by tree for each million years ###########
  
  if (test_type %in% c("less", "two.sided")) {
    is.it.lower.results <- purrr::map2(cleaned_trees, rank.envelope.results,
                                     checkIsLower)
    lower.runs <- purrr::map(is.it.lower.results, findLowerRuns)
    names(is.it.lower.results) <- as.numeric(seq(1:length(is.it.lower.results)))
    all.lower.timepoints <- purrr::map_dfr(is.it.lower.results, rbind, .id = 'tree.index')
    all.lower.timepoints$time.bins <- cut(all.lower.timepoints$actual.time, 
                                    breaks = c(0, seq(1,maximum.root.age+1)), 
                                    labels=FALSE)
    
    # Then, aggregate the significantly lower points by tree and by bin
    
    agg.lower.timepts <- aggregate(x= list(is.lower = all.lower.timepoints$is.it.lower),
                                   by = list(Tree = all.lower.timepoints$tree.index, 
                                             Bin = all.lower.timepoints$time.bins), 
                                   FUN = sum)
    
    # Set it so multiple points from the same tree in the same bin only count as 1
    agg.lower.timepts$is.lower[which(agg.lower.timepts$is.lower > 0)] <- 1   
    # Get dataframe of only significant points
    lower.timepts <- agg.lower.timepts[which(agg.lower.timepts$is.lower == 1),]
    if (nrow(lower.timepts) > 1) {
      
      # Aggregate by bin
      lower.bins <- aggregate(x=list(NumTrees = lower.timepts$is.lower), 
                              by = list(Bin = lower.timepts$Bin), FUN = sum)
    } else{
      lower.bins <- NA
    }
  }
  
  if (test_type %in% c("greater", "two.sided")) {
    is.it.higher.results <- purrr::map2(cleaned_trees, rank.envelope.results,
                                      checkIsHigher)
    higher.runs <- purrr::map(is.it.higher.results, findHigherRuns)
    names(is.it.higher.results) <- as.numeric(seq(1:length(is.it.higher.results)))
    all.higher.timepoints <- purrr::map_dfr(is.it.higher.results, rbind, .id = 'tree.index')
    all.higher.timepoints$time.bins <- cut(all.higher.timepoints$actual.time, 
                                          breaks = c(0, seq(1,maximum.root.age+1)), 
                                          labels=FALSE)
    # Then, aggregate the significantly higher points by tree and by bin
    agg.higher.timepts <- aggregate(x= list(is.higher = all.higher.timepoints$is.it.higher),
                           by = list(Tree = all.higher.timepoints$tree.index, 
                                     Bin = all.higher.timepoints$time.bins), 
                           FUN = sum)
    # Set it so multiple points from the same tree in the same bin only count as 1
    agg.higher.timepts$is.higher[which(agg.higher.timepts$is.higher > 0)] <- 1   
    # Get dataframe of only significant points
    higher.timepts <- agg.higher.timepts[which(agg.higher.timepts$is.higher == 1),]
    if (nrow(higher.timepts) > 1) {
      # Aggregate by bin
      higher.bins <- aggregate(x=list(NumTrees = higher.timepts$is.higher), 
                            by = list(Bin = higher.timepts$Bin), FUN = sum)
    } else {
      higher.bins <- NA
    }
  }
  
  # NOTE: bin labels mean that the time point was in the interval 
  # (1-binlabel, bin label)    
  
  results <- switch(test_type, 
                    "less" =  list(dtt_results = dtt.results.all,
                                   mdi_output = mdi.output,
                                   rank_env_results = rank.envelope.results,
                                   all_lower_timepoints = all.lower.timepoints,
                                   agg_lower_timepoints = agg.lower.timepts,
                                   lower_bins = lower.bins,
                                   cleaned_trees = cleaned_trees),
                    "greater" = list(dtt_results = dtt.results.all,
                                     mdi_output = mdi.output,
                                     rank_env_results = rank.envelope.results,
                                     all_higher_timepoints = all.higher.timepoints,
                                     agg_higher_timepoints = agg.higher.timepts,
                                     higher_bins = higher.bins,
                                     cleaned_trees = cleaned_trees),
                    "two.sided" = list(dtt_results = dtt.results.all,
                                       mdi_output = mdi.output,
                                       rank_env_results = rank.envelope.results,
                                       all_lower_timepoints = all.lower.timepoints,
                                       agg_lower_timepoints = agg.lower.timepts,
                                       all_higher_timepoints = all.higher.timepoints,
                                       agg_higher_timepoints = agg.higher.timepts,
                                       higher_bins = higher.bins,
                                       lower_bins = lower.bins,
                                       cleaned_trees = cleaned_trees))
  
  results
}

#' Plots proportion of  trees with a timepoint sig. lower than the rank envelope
#' 
#' @param dtt_output An object returned by run_modified_dtt
#' @return a ggplot 

plot_dtt_islower <- function(dtt_output) {

  timepoints.plot <- ggplot(data = dtt_output$lower_bins,
                            aes(x = Bin, y = NumTrees/length(dtt_output$cleaned_trees))) +
    geom_col(position=position_nudge(x=0.5), fill = "grey50") +  # Nudge b/c of bin bounds
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0)) +  
    labs(      labs(x='Millions of years ago (binned)', 
                    y='Proportion of\nPosterior Trees',
                    title = "Timepoints Below the Rank Envelope") ) +
    scale_x_reverse(limits = c(67,0))+
    geom_vline(xintercept = 66, colour = 'red', linetype = "dashed") +
    theme_procB()
  
  timepoints.plot 
}

#' Plots proportion of  trees with a timepoint sig. higher than the rank envelope
#' 
#' @param dtt_output An object returned by run_modified_dtt
#' @return a ggplot 

plot_dtt_ishigher <- function(dtt_output) {

    timepoints.plot <- ggplot(data = dtt_output$higher_bins,
                              aes(x = Bin, y = NumTrees/length(dtt_output$cleaned_trees))) +
      geom_col(position=position_nudge(x=0.5), fill = "grey50") +  # Nudge b/c of bin bounds
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0)) +  
      labs(x='Millions of years ago (binned)', 
           y='Proportion of\nPosterior Trees',
           title = "Timepoints Above the Rank Envelope") +
      scale_x_reverse(limits = c(67,0))+
      geom_vline(xintercept = 66, colour = 'red', linetype = "dashed") +
      theme_procB()
  
  timepoints.plot 
}


# Functions ########################

#' Runs a CTT and null envelope analysis on a single tree
#' 
#' @param tree a phylo object
#' @param cat_trait a factor vector with ti
#' p labels as names
#' @param n_sim the number of simulations used to generate the null envelope
do_ctt <- function(tree, i, cat_trait, n_sim = 1000) { 
  tree <- ladderize(treedata(tree, cat_trait)$phy)
  trait <- cat_trait[tree$tip.label]
  trees_simmap <- make.simmap(tree, cat_trait , model="SYM", nsim=n_sim)
  
  # Generate (or simulates) a 'changes through time' plot from a set of stochastic map character histories
  object <- ctt(trees_simmap) # segment number to match nulo
  
  # use the function sim.multiCtt to simulate various rather than a single CTT
  Q<-trees_simmap[[1]]$Q    # Note, Qs for all simulated trees are the same. 
  nulo<-sim.multiCtt(tree,Q,nsim=n_sim)
  save(object,nulo,file = paste0("data/cttout/",i,".RData"))
  list(object = object, nulo = nulo)
  
}

# Find nulo bounds
#' For a set of null simulations and a ctt from data, this function returns the 
#' bounds of the null envelope and the number or rate of changes from the data ctt 
#' 
#' @param null A multiCtt object (from null simulations)
#' @param actual A ctt object (from the data)
#' @param type The kind of comparison (number of changes, or rate of changes)
#' @param alpha The alpha value for a two tailed test
#' 
#' @return A list of time bins, env. bounds, and data values
check_ctt_env <- function(ctt_output, type = c("number", "rate"), alpha = 0.05) { 
  actual <- ctt_output$object
  null <- ctt_output$nulo
  # Find envelope parameters
  segments <- null[[1]]$segments
  nchanges <- sapply(null, function(x) x$nchanges)
  edge.length <- sapply(null, function(x) x$edge.length)
  null_summary <- list(segments = segments, 
                       nchanges = rowMeans(nchanges),
                       edge.length=rowMeans(edge.length),
                       tree=null[[1]]$tree)
  class(null_summary)<-"ctt"
  
  # Find indices of null bound cases
  lower<-max(floor(alpha/2*length(null)),1)  # evals to 2 for alpha = 0.05
  upper<-min(ceiling((1-alpha/2)*length(null)),ncol(nchanges))  # evals to 98 for alpha = 0.05
  polygon_bin_bounds<-max(nodeHeights(null[[1]]$tree))-as.vector(t(segments)) # Time points of bins, duplicated?
  polygon_bin_bounds<-c(polygon_bin_bounds,polygon_bin_bounds[length(polygon_bin_bounds):1])
  
  if (type == "number") {
    y.lower <- apply(nchanges, 1, sort)[lower, ] 
    y.upper <- apply(nchanges,1,sort)[upper,]
  } else if (type == "rate") {
    y.lower <- apply(nchanges/edge.length, 1, sort)[lower, ] 
    y.upper <- apply(nchanges/edge.length,1,sort)[upper,]
  } else { stop("Please provide valid arg(type)")  }
  
  y.lower<-as.vector(rbind(y.lower,y.lower))
  y.upper<-as.vector(rbind(y.upper,y.upper))
  yy<-c(y.lower,y.upper[length(y.upper):1])
  
  # Find parameters from data
  h <- max(nodeHeights(actual$tree))
  time_bin_bounds <- h - as.vector(t(actual$segments)) 
  data_y_vals <- if(type=="number") rbind(actual$nchanges,actual$nchanges) else 
    rbind(actual$nchanges/actual$edge.length,actual$nchanges/actual$edge.length)
  
  list(time_bin_bounds = time_bin_bounds, 
       polygon_bin_bounds= polygon_bin_bounds, 
       polygon_env_bounds = yy, 
       env_maxima = y.upper[c(TRUE,FALSE)],
       env_minima = y.lower[c(TRUE,FALSE)],
       data_y_vals = data_y_vals)
}


# Find nulo bounds
#' For a set of null simulations and a ctt from data, this function returns the 
#' bounds of the null envelope and the number or rate of changes from the data ctt 
#' 
#' @param null A multiCtt object (from null simulations)
#' @param actual A ctt object (from the data)
#' @param type The kind of comparison (number of changes, or rate of changes)
#' @param alpha The alpha value for a two tailed test
#' 
#' @return A list of time bins, env. bounds, and data values
lower_or_higher <- function(ctt_env) {
  data_values <- as.vector(ctt_env$data_y_vals[1,])
  min_values <- ctt_env$env_minima
  max_values <- ctt_env$env_maxima
  
  output <- numeric(length(data_values))
  for (i in 1:20) { 
    if (data_values[i] > max_values[i]) {
      output[i] <- 1
    } else if (data_values[i] < min_values[i]) {
      output[i] <- -1
    } else {
      output[i] <- 0
    }
  }
  
  output
}


### REWRITE OF CTT ###################################################


# computing the mean number of character changes through time from a set of stochastic map trees
## written by Liam J. Revell 2017

ctt<-function(trees,...){
  if(!(inherits(trees,"multiSimmap")))
    stop("trees should be an object of class \"multiSimmap\".")
  tree<-as.phylo(trees[[1]])
  segments <- max(nodeHeights(tree))
  changes<-sapply(trees,getChanges)
  h<-max(nodeHeights(tree))
  b<-segments
  segs<-cbind(seq(0,h-h/b,h/b),
              seq(1/b*h,h,h/b))
  nchanges<-rep(0,b)
  for(i in 1:length(changes)){
    for(j in 1:length(changes[[i]])){
      ind<-which((changes[[i]][j]>segs[,1])+
                   (changes[[i]][j]<=segs[,2])==2)
      nchanges[ind]<-nchanges[ind]+1/length(changes)
    }
  }
  LTT<-ltt(tree,plot=FALSE)
  LTT<-cbind(LTT$ltt[2:(length(LTT$ltt)-1)],
             LTT$times[2:(length(LTT$ltt)-1)],
             LTT$times[3:length(LTT$ltt)])
  ii<-1
  edge.length<-rep(0,b)
  for(i in 1:nrow(segs)){
    done.seg<-FALSE
    while(LTT[ii,2]<=segs[i,2]&&done.seg==FALSE){
      edge.length[i]<-edge.length[i]+
        LTT[ii,1]*(min(segs[i,2],LTT[ii,3])-
                     max(segs[i,1],LTT[ii,2]))
      if(LTT[ii,3]>=segs[i,2]) done.seg<-TRUE
      if(LTT[ii,3]<=segs[i,2]) ii<-if(ii<nrow(LTT)) ii+1 else ii
    }
  }
  object<-list(segments=segs,nchanges=nchanges,edge.length=edge.length,tree=tree)
  class(object)<-"ctt"
  object
}	

plot.ctt<-function(x,...){
  h<-max(nodeHeights(x$tree))
  args<-list(...)
  if(!is.null(args$type)){ 
    type<-args$type
    args$type<-NULL
  } else type<-"rate"
  if(!is.null(args$show.tree)){
    show.tree<-args$show.tree
    args$show.tree<-NULL
  } else show.tree<-FALSE
  if(!is.null(args$add)){ 
    add<-args$add
    args$add<-NULL
  } else add<-FALSE
  if(is.null(args$ylim)) 
    args$ylim<-if(type=="number")c(0,max(x$nchanges)) else 
      c(0,max(x$nchanges/x$edge.length))
  if(is.null(args$xlim))
    args$xlim<-c(max(x$segments),min(x$segments))
  if(is.null(args$lwd)) args$lwd<-2
  if(is.null(args$xlab)) args$xlab<-"time since the present"
  if(is.null(args$ylab)) 
    args$ylab<-if(type=="number") "mean number of changes" else
      "mean number of changes / total edge length"
  args$type<-"l"
  args$x<-h-as.vector(t(x$segments))
  args$y<-if(type=="number") rbind(x$nchanges,x$nchanges) else 
    rbind(x$nchanges/x$edge.length,x$nchanges/x$edge.length)
  if(!add) do.call(plot,args)
  else do.call(lines,args)
  if(show.tree) plotTree(x$tree,add=TRUE,ftype="off",lwd=1,
                         color=make.transparent("blue",0.1),mar=par()$mar,
                         direction="leftwards",xlim=args$xlim)
}

sim.ctt<-function(tree,Q,anc=NULL,nmaps=100,...){
  x<-as.factor(sim.history(tree,Q,anc=anc,message=FALSE)$states)
  while(length(levels(x))!=ncol(Q)) 
    x<-as.factor(sim.history(tree,Q,anc=anc,message=FALSE)$states)
  flush.console()
  cat("Starting stochastic mapping with simulated data vector.... ")
  flush.console()
  trees<-make.simmap(tree,x,Q=Q,nsim=nmaps,message=FALSE)
  cat("Done.\n")
  flush.console()
  ctt(trees)
}

sim.multiCtt<-function(tree,Q,anc=NULL,nmaps=100,nsim=100,...){
  object<-replicate(nsim,sim.ctt(tree,Q,anc=anc,nmaps=nmaps,...),simplify=FALSE)
  class(object)<-"multiCtt"
  object
}

getChanges<-function(tree){
  states<-sort(unique(getStates(tree)))
  nc<-sapply(tree$maps,length)-1
  b<-which(nc>0)
  nc<-nc[b]
  xx<-vector()
  H<-nodeHeights(tree)
  for(i in 1:length(b)){
    for(j in 1:nc[i]){
      ss<-names(tree$maps[[b[i]]])[j+1]
      x<-rep(H[b[i],1]+cumsum(tree$maps[[b[i]]])[j],2)
      xx<-c(xx,setNames(x[1],
                        paste(names(tree$maps[[b[i]]])[j:(j+1)],
                              collapse="->")))
    }
  }
  xx
}

print.ctt<-function(x,...){
  cat("Object of class \"ctt\" consisting of:\n")
  cat("   (1) a matrix (segments) with the beginning & ending time of each segment.\n")
  cat("   (2) a vector (nchanges) with the mean number of changes in each segment.\n")
  cat("   (3) a vector (edge.length) containing the total edge length of each segement.\n")
  cat("   (4) an object of class \"phylo\".\n\n")
}

print.multiCtt<-function(x,...){
  cat(paste(length(x),"objects of class \"ctt\" in a list.\n\n"))
}

plot.multiCtt<-function(x,...){
  if(hasArg(alpha)) alpha<-list(...)$alpha
  else alpha<-0.05
  segments<-x[[1]]$segments
  nchanges<-sapply(x,function(x) x$nchanges)
  if(hasArg(type)) type<-list(...)$type
  else type<-"rate"
  edge.length<-sapply(x,function(x) x$edge.length)
  obj<-list(segments=segments,nchanges=rowMeans(nchanges),
            edge.length=rowMeans(edge.length),tree=x[[1]]$tree)
  class(obj)<-"ctt"
  lower<-max(floor(alpha/2*length(x)),1)
  upper<-min(ceiling((1-alpha/2)*length(x)),ncol(nchanges))
  xx<-max(nodeHeights(x[[1]]$tree))-as.vector(t(segments))
  xx<-c(xx,xx[length(xx):1])
  y.lower<-if(type=="number") apply(nchanges,1,sort)[lower,] else
    if(type=="rate") apply(nchanges/edge.length,1,sort)[lower,]
  y.upper<-if(type=="number") apply(nchanges,1,sort)[upper,] else
    if(type=="rate") apply(nchanges/edge.length,1,sort)[upper,]
  y.lower<-as.vector(rbind(y.lower,y.lower))
  y.upper<-as.vector(rbind(y.upper,y.upper))
  yy<-c(y.lower,y.upper[length(y.upper):1])
  args<-list(...)
  if(!is.null(args$alpha)) args$alpha<-NULL
  if(is.null(args$col)) args$col<-"blue"
  if(is.null(args$ylim)) args$ylim<-range(yy)
  args$x<-obj
  do.call(plot,args)
  polygon(xx,yy,col=make.transparent("grey",0.4),border=0)
}



