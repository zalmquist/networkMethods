#' Lab 1
#'
#' \code{lab1} returns simple description of the first lab in Soc 88412.
#'
#' This is a simple description of the lab, use help(lab1) to review the
#' examples for this lab
#'
#' @param Empty
#' @examples
#' ## Practice R code
#' matrix(c(1,2,3,4),nc=2,nr=2)
#' ## More to come!
#' \dontrun{
#'
#' }
lab1 <- function() {"This is Lab one."}


#' Lab 2
#'
#' \code{lab2} returns simple description of the first lab in Soc 88412.
#'
#' This is a simple description of the lab, use help(lab1) to review the
#' examples for this lab
#'
#' @param Empty
#' @examples
#' library(sna)                      # Load the sna library
#' data(contig_1993)
#' data(mids_1993)
#'
#' #-Network visualization with gplot----------------------------------------------
#' #
#' # Begin by plotting contiguity among nations in 1993 (from the Correlates of
#' # War project)
#' gplot(contig_1993)                              # The default visualization
#' gplot(contig_1993, usearrows=FALSE)             # Turn off arrows manually
#' gplot(contig_1993, gmode="graph")               # Can also tell gplot the data
#' # is undirected
#'
#' # We can add labels to the vertices - network.vertex.names reports them
#' gplot(contig_1993, gmode="graph", label=network.vertex.names(contig_1993))

#' # This plot is too large/dense for the default settings to work.  Let's refine
#' # them.
#' gplot(contig_1993, gmode="graph", boxed.labels=FALSE, label.cex=0.5,
#'       label.col=4, label=network.vertex.names(contig_1993)) # Shrink labels,
#' # remove boxes, recolor
#'
#' # Here's an example of directed data - militarized interstate disputes (MIDs)
#' # for 1993
#' gplot(mids_1993, boxed.labels=FALSE, label.cex=0.5, label.col=4,
#'       label=network.vertex.names(mids_1993))       # Basic display, with labels
#'
#' # All those isolates can get in the way - we can suppress them using
#' # displayisolates
#' gplot(mids_1993, displayisolates=FALSE, boxed.labels=FALSE, label.cex=0.5,
#'       label.col=4,label=network.vertex.names(mids_1993))
#'
#' # The default layout algorithm is that of Frutchterman-Reingold (1991), can use
#' # others
#' gplot(mids_1993, displayisolates=FALSE, boxed.labels=FALSE, label.cex=0.5,
#'       label.col=4,label=network.vertex.names(mids_1993),
#'       mode="circle")   # The infamous circle
#' gplot(mids_1993, displayisolates=FALSE, boxed.labels=FALSE, label.cex=0.5,
#'       label.col=4,label=network.vertex.names(mids_1993),
#'       mode="mds")      # MDS of position similarity
#'
#' # When a layout is generated, the results can be saved for later reuse:
#' coords <- gplot(contig_1993)                  # Capture the magic of the moment
#' coords                                        # Show the vertex coordinates
#'
#' #Saved (or a priori) layouts can be used via the coord argument:
#' gplot(mids_1993, boxed.labels=FALSE, label.cex=0.5, label.col=4, coord=coords,
#'       label=network.vertex.names(mids_1993))        # Relive the magic
#'
#' # When the default settings are insufficient, interactive mode allows for
#' # tweaking
#' coords <- gplot(contig_1993, interactive=TRUE)   # Modify and save
#' gplot(contig_1993, coord=coords, gmode="graph")  # Should reproduce the modified
#' # layout
#'\dontrun{
#' #
#' #Three-dimensional visualization with gplot3d (requires the rgl package)--------
#' # Note that this requires X11
#'
#' # Note: if you haven't done so, you can install the rgl package by typing
#' # install.packages("rgl") at the command prompt (assuming you are online).
#' # If that package is not installed, you'll get a boring error message instead
#' # of exciting visualization.
#'
#' gplot3d(contig_1993, label=network.vertex.names(contig_1993))   # Experience
#' # the future!
#' # Other layouts are possible here, too:
#' gplot3d(contig_1993, label=network.vertex.names(contig_1993),mode="kamadakawai")
#'
#' #For more information....
#' ?gplot
#' ?gplot.layout
#'
#' #For more information....
#' ?gplot3d
#' ?gplot3d.layout
#' }
lab2 <- function() {"This is Lab two."}


#' Lab 3
#'
#' \code{lab3} returns simple description of the third lab in Soc 88412.
#'
#' This is a simple description of the lab, use help(lab3) to review the
#' examples for this lab
#'
#' @param Empty
#' @examples
#' library(sna)                      # Load the sna library
#' data(mids_1993)
#' data(contig_1993)
#' #
#' #Basic centrality indices: degree, betweenness, and closeness-------------------
#' #
#' # We begin with the simplest case: degree
#' degree(mids_1993)                                        # Default: total degree
#' ideg <- degree(mids_1993, cmode="indegree")              # Indegree for MIDs
#' odeg <- degree(mids_1993, cmode="outdegree")             # Outdegree for MIDs
#' all(degree(mids_1993) == ideg+odeg)                      # In + out = total?
#'
#' # Once centrality scores are computed, we can handle them using standard R
#' # methods:
#' plot(ideg, odeg, type="n", xlab="Incoming MIDs", ylab="Outgoing MIDs")
#' abline(0, 1, lty=3)
#' text(jitter(ideg), jitter(odeg), network.vertex.names(contig_1993), cex=0.75,
#'      col=2)   #Plot index by odeg
#'
#' \dontrun{
#' #Plot simple histograms of the degree distribution:
#' pdf("simpleHistIndOut.pdf")
#' par(mfrow=c(2,2))                                       # Set up a 2x2 display
#' hist(ideg, xlab="Indegree", main="Indegree Distribution", prob=TRUE)
#' hist(odeg, xlab="Outdegree", main="Outdegree Distribution", prob=TRUE)
#' hist(ideg+odeg, xlab="Total Degree", main="Total Degree Distribution",
#'      prob=TRUE)
#' dev.off()
#'}
#'
#' # Centrality scores can also be used with other sna routines, e.g., gplot
#' gplot(mids_1993, vertex.cex=(ideg+odeg)^0.5/2, vertex.sides=50,
#'       boxed.labels=FALSE,label.cex=0.4,
#'       vertex.col=rgb(odeg/max(odeg),0,ideg/max(ideg)),
#'       label=network.vertex.names(mids_1993))
#'
#' # Betweenness and closeness are also popular measures
#' bet <- betweenness(contig_1993, gmode="graph")       # Geographic betweenness
#' bet
#' gplot(contig_1993, vertex.cex=sqrt(bet)/25, gmode="graph")   # Use w/gplot
#' clo <- closeness(contig_1993)                        # Geographic closeness
#' clo                                                  # A large world after all?
#'
#' # Can use sna routines to explore alternatives to the common measures....
#' closeness2 <- function(x){            # Create an alternate closeness function!
#'   geo <- 1/geodist(x)$gdist         # Get the matrix of 1/geodesic distance
#'   diag(geo) <- 0                    # Define self-ties as 0
#'   apply(geo, 1, sum)                # Return sum(1/geodist) for each vertex
#' }
#' clo2 <- closeness2(contig_1993)       # Use our new function on contiguity data
#' hist(clo2, xlab="Alt. Closeness", prob=TRUE)    # Much better behaved!
#' cor(clo2, bet)                                  # Correlate with betweenness
#' plot(clo2, bet)                            # Plot the bivariate relationship
#'
#' \dontrun{
#' #For more information....
#' ?betweenness
#' ?bonpow
#' ?closeness
#' ?degree
#' ?evcent
#' ?graphcent
#' ?infocent
#' ?prestige
#' ?stresscent
#' }
#' #
#' #Simple hypothesis tests for NLIs----------------------------------------------
#' #
#' library(network)                               #Load network if needed
#' data(emon)                                     #Load Drabek et al. data
#'
#' #Extract ties from the Cheyenne EMON communicating at least "every few hours"
#' g<-as.sociomatrix(emon[[1]],"Frequency")       #Need to get the frequency info
#' g<-symmetrize((g>0)&(g<4))                     #Note the reverse coding!
#'
#' #Get some potential covariates
#' drs<-emon[[1]]%v%"Decision.Rank.Score"         #Get decision rank (see man page)
#' crs<-emon[[1]]%v%"Command.Rank.Score"          #Get command rank
#'
#' #Calculate some basic centrality measures
#' deg<-degree(g,gmode="graph")
#' bet<-betweenness(g,gmode="graph")
#' clo<-closeness(g,gmode="graph")
#'
#' #Raw correlations
#' cor(cbind(deg,bet,clo),cbind(drs,crs))
#'
#' #Classical tests (using asymptotic t distribution)
#' cor.test(deg,drs)
#' cor.test(bet,drs)
#' cor.test(clo,drs)
#'
#' #Permutation tests
#' perm.cor.test<-function(x,y,niter=5000){  #Define a simple test function
#'   c.obs<-cor(x,y,use="complete.obs")
#'   c.rep<-vector()
#'   for(i in 1:niter)
#'     c.rep[i]<-cor(x,sample(y),use="complete.obs")
#'   cat("Vector Permutation Test:\n\tObserved correlation: ",c.obs,"\tReplicate quantiles (niter=",niter,")\n",sep="")
#'   cat("\t\tPr(rho>=obs):",mean(c.rep>=c.obs),"\n")
#'   cat("\t\tPr(rho<=obs):",mean(c.rep<=c.obs),"\n")
#'   cat("\t\tPr(|rho|>=|obs|):",mean(abs(c.rep)>=abs(c.obs)),"\n")
#'   invisible(list(obs=c.obs,rep=c.rep))
#' }
#' perm.cor.test(deg,drs)                     #Non-parametric tests of correlation
#' perm.cor.test(bet,drs)
#' perm.cor.test(clo,drs)
#'
#' \dontrun{
#' #For more information....
#' ?emon
#' ?cor.test
#' ?t.test
#' ?sample
#' }
#' #
#' #Using NLIs as regression covariates--------------------------------------------
#' #
#' pstaff<-emon[[1]]%v%"Paid.Staff"                     # Get more EMON covariates
#' vstaff<-emon[[1]]%v%"Volunteer.Staff"
#' govt<-((emon[[1]]%v%"Sponsorship")!="Private")
#'
#' #Very simple model: decision rank is linear in size, degree, and govt status
#' mod<-lm(drs~deg+pstaff+vstaff+govt)
#' summary(mod)
#' anova(mod)                                            #Some useful lm tools
#' AIC(mod)
#'
#' #Does total size change the picture?
#' mod2<-lm(drs~deg+I(pstaff+vstaff)+govt)                #Pre-add sizes
#' summary(mod2)
#'
#' #Try with alternative measures....
#' mod3<-lm(drs~bet+pstaff+vstaff+govt)                   #Betweenness
#' summary(mod3)
#' mod4<-lm(drs~clo+pstaff+vstaff+govt)                   #Closeness
#' summary(mod4)
#' AIC(mod,mod3,mod4)                                     #Closeness wins!
#'
#'
#' \dontrun{
#'#For more information....
#' ?lm
#' ?anova
#' ?AIC
#' }
lab3 <- function() {"This is Lab three."}

