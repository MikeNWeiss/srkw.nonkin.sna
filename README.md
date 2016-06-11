# srkw.nonkin.sna
An R package for analyzing SRKW social networks using ERGMs

These are specific functions used in my thesis to analyze the non-kinship bonds found in the southern resident population of 
killer whales in 2015.

The main functions are the simulate.nonkin and whale.community.gof functions.
The simulate.nonkin function repeatedly simulates from an ERGM of the whale network and tabulates the different 
kinds of cross-pod bonds in each simulation.
The whale.community.gof function is a goodness of fit function, much like the one built into the ergm package. However, 
instead of the default statistics used in the gof.ergm function,  this function simulates from a model and records number of
connected components, and also performs optimal-modularity community detection using the igraph package, recording the number
of communities and the modularity. These statistics are then compared to the original network to determine how well the model
reproduces large-scale community features of the network. For this function, there is the also the NetGraph function that
turns the whale network object into an igraph object for community detection functions.

In addition, the package contains a function to plot the results of a simulation (simresults.vioplot) and a function with
no arguments that loads the 2015 SRKW preferred associates network (Load2015Net), which is stored as a .csv in the package. 
In addition, the whale.ergm function produces an ergm object with the network statistics used in my analysis.

Many of these functions, particularly the community goodness of fit function, can likely be pruned/extended to be more
generally usable, but as they are these are specific to the network I worked with. 
