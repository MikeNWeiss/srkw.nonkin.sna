#' Load the 2015 SRKW Preferred Association Network
#'
#' @return A network object of the 2015 preferred associations with nodal attributes
#'
#' @details The 2015 preferred association network is based on sightings data of SRKW matrilines from May to September of 2015. Preferred associations were determined by computing half weigth indexes of association and then using a randomization test to determine significant dyads (Bejder et al. 1998).
#'
#' @references Bejder L., Fletcher D., and Brager S. 1998. "A method for testing association patterns of animals." Animal Behavior, vol. 56, no. 3, pp.719-725
Load2015Net <- function(){ #load in the 2015 whale network
  adjacency.matrix <- as.matrix(read.csv("~/srkw.nonkin.sna/data/2015.adjmat.csv", row.names=1))
  n <- network(adjacency.matrix, directed=F) #Loads the adjacency matrix as a network
  # add node attributes
  n%v%"male" <- c('m', 'nm', 'm', 'nm', 'nm', 'm', 'nm', 'nm', 'm','m','m','m','m','m','nm','m','nm','nm','nm','nm','m','m','m','m','m')
  n%v%"prf" <- c("nf", "f", "f", "nf", "nf", "nf", "f", "f", "f", "nf", "nf", "nf", "nf", "nf", "f", "nf", "f","f", "nf", "nf", "nf", "nf", "nf", "nf", "nf")
  n%v%"color" <- c("blue", "red", "purple", "white", "white", "blue", "red", "red", "purple", "blue", "blue", "blue", "blue", "blue", "red", "blue", "red", "red", "white", "white", "blue", "blue", "blue", "blue", "blue")
  n%v%"shape" <- c(rep(1000, 7), rep(4, 5), rep(3, 13))
  n%v%"pod" <- c(rep("J", 7),rep("K", 5), rep("L", 13))
  n%v%"n" <- c(3,5,6,6,3,3,1,5,7,4,2,1,6,2,1,2,9,5,3,2,1,1,1,1,1)
  n #return the network with attributes
}

#' Create the ergm used in the original analysis
#'
#' @param net A network object containing nodal attributes pod (pod identity), male (presence of an adult male), and prf (presence of post-reproductive females), such as returned by Load2015Net.
#' @return An ergm object containing parameters for edges, nodematch("pod"), nodefactor(c("male", "prf")), gwesp(0.5, fixed=T), and isolates (which will be fixed at -Inf for the 2015 network)
#' @details This function is essentially a shortcut to generate the ergm used in the analysis of the 2015 network that I performed. The model terms used are edges, uniform pod homophily, node factor by social roles, isolates (fixed at -Inf), and GWESP with alpha=0.5. The MCMC sample size and interval are both set to 5,000.
whale.ergm <- function(net){
  model <- ergm(net ~ edges + nodematch("pod") + nodefactor(c("male", "prf")) + gwesp(0.5, fixed=T) + isolates, control=control.ergm(MCMC.samplesize=5000, MCMC.interval=5000)) #parameter estimation for whale network
  model
}

#' Converts the whale network from a network object to a graph object
#'
#' @param net A network object, containing the nodal attributes used by Load2015Net()
#' @return A graph object version of the network
#' @details This function converts networks to a graph objects to allow for community detection using the igraph package.
WhaleNetGraph <- function(net){#turn whale network into a graph object for community detection
  g <- graph.adjacency(as.matrix(net), mode="undirected") #Make the graph
  V(g)$male <- ifelse(net%v%"male"=="m", 1, 0) #Import vertex attributes
  V(g)$prf <- ifelse(net%v%"prf"=="f", 1, 0)
  V(g)$n <- net%v%"n"
  V(g)$pod <- net%v%"pod"
  g #Return graph object
}

#' Goodness of fit diagnostics for ERGM based on community structure
#'
#' @param model An ergm object
#' @param sample.size The number of simulations to perform. Defaults to 500
#' @details This function does a goodness of fit test for an ergm, based on the general idea used in the ergm gof function, but incorporating large scale network characteristics, like number of connected components, communtiy structure, and modularity. It is currently written for use
#' @return A plot of densities of number of components, communities, and modularity based on the ergm, compared to initial network.
community.gof <- function(model, sample.size=500){ #In addition to the gof function built into the ergm package, this function was written to assess how well a model reproduces the community structure, modularity, and number of connected components in a network. This function, unlike most functions in this package, is fairly general.
  output <- matrix(nrow=sample.size, ncol=3)
  colnames(output) <- c("n.components", "n.communities", "modularity")
  for(i in 1:sample.size){
    sim.net <- simulate(model)
    sim.g <- graph.adjacency(as.matrix(sim.net), mode="undirected")
    comm <- optimal.community(sim.g)$membership
    output[i, "n.communities"] <- max(comm) #Number of communities
    output[i, "modularity"] <- modularity(sim.g, comm) #Modularity of community division
    output[i, "n.components"] <- components(sim.g)$no # Number of connected Components
  }
  plot(density(output[,"n.communities"], bw=0.3), main="Communities", xlab="n.communities", ylab="density")
  abline(v=max(optimal.community(NetGraph(model$network))$membership))
  plot(density(output[,"n.components"], bw=0.3), main="Components", xlab="n.components", ylab="density")
  abline(v=communities(NetGraph(model$network))$no)
  plot(density(output[,"modularity"]), main="Modularity", xlab="modularity", ylab="density")
  abline(v=modularity(optimal.community(NetGraph(model$network))))
}

#' Test for unusually high number of cross-pod dyad types based on an ergm via simulation
#'
#' @param model The ergm to simulate from
#' @param reps The number of networks to simulate. Defaults to 10,000.
#' @return A matrix containing the number of total cross-pod dyads, and dyads of each type, for each simulation
simulate.nonkin <- function(model, reps=10000){#This function simulates from an ergm of the whale network and stores information about the types of cross-pod dyads that occur. The "model" argument indicates which ergm to simulate from. The "reps" argument defines number of simulations. The default is 10,000 simulations.
  output <- matrix(nrow=reps, ncol=4) #Make your output matrix
  colnames(output) <- c("cp.tot", "cp.tandem", "cp.m.empty",  "cp.m.m") #Name your columns
  for(i in 1:reps){
    sim.net <- simulate(model) #Simulate from model
    output[i,"cp.tot"] <- summary(sim.net~nodemix("pod"))[2] + summary(sim.net~nodemix("pod"))[4] + summary(sim.net~nodemix("pod"))[5] #Total number of cross-pod bonds in simulated network
    summary <- summary(sim.net~nodemix(c("pod", "prf", "male"), base = c(1:10, 15, 20, 21, 26, 27, 28, 36, 44, 45, 53:55))) #Make a big summary of just cross-pod bonds
    output[i, "cp.tandem"] <-  summary[7]+summary[10]+summary[15]+summary[19]+summary[21]+summary[25] #Tandem
    output[i, "cp.m.empty"] <- summary[12]+summary[23]+summary[29]+summary[33] #M-Empty
    output[i, "cp.m.m"] <- summary[11]+summary[22]+summary[26] #M-M
  }
  output #Return output matrix after all reps complete
}

#' Produce a violin plot of the results of simulate.nonkin
#'
#' @param simresults A matrix produced by simulate.nonkin
#' @param point.col The color of the points used to indicate observed values. Defaults to red.
#' @param vio.col The color of the violin plots indicating simulation results. Defaults to grey.
#' @return A violin plot of the results, with observed values plotted in red
simresults.vioplot <- function(simresults, point.col="red", vio.col="grey"){ #Plot results of simulation. simresults is a matrix produced by simulate.nonkin()
  simresults.sub <- simresults[simresults[,"cp.tot"]!=0,] #Remove cases with no cross-pod ties
  p.tandem <- simresults.sub[,"cp.tandem"]/simresults.sub[,"cp.tot"] #Porion tandem
  p.m.empty <- simresults.sub[,"cp.m.empty"]/simresults.sub[,"cp.tot"] #Portion m.empty
  p.m.m <- simresults.sub[,"cp.m.m"]/simresults.sub[,"cp.tot"] #Portion m.m
  vioplot(p.tandem, p.m.empty, p.m.m, wex=0.8, col=vio.col) #Violin plots
  points(c(1:3), c(4/9, 4/9, 1/9), cex=1.5, col=point.col, pch=19) #Add points corresponding to observed values
}
