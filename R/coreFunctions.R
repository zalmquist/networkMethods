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
#' a <- 3							# assignment
#' a								# evaluation
#'
#' sqrt(a)						# perform an operation
#' b <- sqrt(a)						# perform operation and save
#' b
#'
#' a == a                                         # A is A?
#' a != b                                         # A is not B
#'
#' ls()	# list objects in global environment
#' #help(sqrt)						# help w/ functions
#' #?sqrt							# same thing
#' #help.start()						# lots of help
#' #help.search("sqrt")					# what am I looking for?
#' #apropos("sqr")					# it's on the tip of my tongue...
#'
#' rm(a)							# remove an object
#' # Creating vectors using the concatenation operator
#' a <- c(1,3,5)						# create a vector by concatenation
#' a
#' a[2]								# select the second element
#' b <- c("one","three","five")			# also works with strings
#' b
#' b[2]
#' a <- c(a,a)						# can apply recursively
#' a
#' a <- c(a,b)						# mixing types - who will win?
#' a								# there can be only one!
#'
#' # Sequences and replication
#' a <- seq(from=1,to=5,by=1)				# from 1 to 5 the slow way
#' b <- 1:5						# a shortcut!
#' a==b								# all TRUE
#' rep(1,5)						# a lot of 1s
#' rep(1:5,2)						# repeat an entire sequence
#' rep(1:5,each=2)					# same, but element-wise
#' rep(1:5,times=5:1)					# can vary the count of each element
#'
#' # Any and all (with vectors)
#' a <- 1:5						# create a vector
#' a>2								# some TRUE, some FALSE
#' any(a>2)						# are any elements TRUE?
#' all(a>2)						# are all elements TRUE?
#'
#' # From vectors to matrices
#' a <- matrix(1:25, nr=5, nc=5)			# create a matrix the "formal" way
#' a
#' a[1,2]							# select a matrix element (two dimensions)
#' a[1,]							# shortcut for ith row
#' all(a[1,]==a[1,1:5])				# show the equivalence
#' a[,2]							# can also perform for columns
#' a[2:3,3:5]						# select submatrices
#' a[-1,]							# nice trick: negative numbers omit cells!
#' a[-2,-2]						# get rid of number two
#'
#' b <- cbind(1:5,1:5)					# another way to create matrices
#' b
#' d <- rbind(1:5,1:5)					# can perform with rows, too
#' d
#' try(cbind(b,d))						# no go: must have compatible dimensions!
#' dim(b)							# what were those dimensions, anyway?
#' dim(d)
#' NROW(b)
#' NCOL(b)
#' cbind(b,b)						# here's a better example
#'
#' t(b)								# can transpose b
#' cbind(t(b),d)						# now it works
#' # Most arithmetic operators are applied element-wise
#' a <- 1:5
#' a + 1							# addition
#' a * 2							# multiplication
#' a / 3							# division
#' a - 4							# subtraction
#' a ^ 5							# you get the idea...
#'
#' a + a							# also works on pairs of vectors
#' a * a
#' a %*% a						# note, not element-wise!
#' a + 1:6						# problem: need same length
#'
#' a <- rbind(1:5,2:6)					# same principles apply to matrices
#' b <- rbind(3:7,4:8)
#' a + b
#' a / b
#'
#' a %*% t(b)						# matrix multiplication
#'
#' # Logical operators (generally) work like arithmetic ones
#' a > 0
#' a == b
#' a != b
#' !(a == b)
#' (a>2) | (b>4)
#' (a>2) & (b>4)
#' (a>2) || (b>4)					# beware the "double-pipe"!
#' (a>2) && (b>4)					# (and the "double-ampersand"!)
#'
#' # Ditto for many other basic transformations
#' log(a)
#' exp(b)
#' sqrt(a+b)						# note that we can nest statements!
#' log((sqrt(a+b)+a)*b)				# as recursive as we wanna be
#' # R has many other data types.  One important type is the list.
#' a <- list(1:5)
#' a								# not an ordinary vector...
#' a <- list(1:5,letters[1:3])			# can we mix types and lengths?
#' a								# yes!
#' b <- matrix(1:3,3,3)
#' a <- list(1:5,letters[1:3],b)			# anything can be stuffed in here
#' a
#' a[[1]]							# retrieve the first item
#' a[[2]][3]						# the letter "c"
#' (a[[3]])[1,3]						# it's really just recursion again
#' a <- list(boo=1:4,hoo=5)				# list elements are often named
#' names(a)						# get the element names
#' a[["boo"]]						# ask for it by name
#' a$hoo							# use "$" to get what you want
#' # a+3								# whoops - not a vector!
#' a[[1]]+3						# that works
#' a[[2]] <- a[[2]]*4					# can also perform assignment
#' a$woo <- "glorp"					# works with "$"
#' a[["foo"]] <- "shazam"				# prolonging the magic
#' a
#'
#' # Another useful generalization: the data frame
#' d <- data.frame(income=1:5,sane=c(T,T,T,T,F),name=LETTERS[1:5])  # Store multiple types
#' d
#' d[1,2]							# acts a lot like a matrix!
#' d[,1]*5
#' d[-1,]
#' names(d)						# also acts like a list
#' d[[2]]
#' d$sane[3]<-FALSE
#' d
#' d[2,3]							# hmm - our data got factorized!
#' d$name <- LETTERS[1:5]				# eliminate evil factors by overwriting
#' d[2,3]
#' d <- data.frame(income=1:5,sane=c(T,T,T,T,F),name=LETTERS[1:5],stringsAsFactors=FALSE)
#' d								# another way to fix it
#'
#' d <- as.data.frame(cbind(1:5,2:6))		# can create from matrices
#' d
#' is.data.frame(d)					# how can we tell it's not a matrix?
#' is.matrix(d)						# the truth comes out
#'
#' # When two dimensions are not enough: arrays
#' a <- array(1:18, dim=c(2,3,3))			# now in 3D
#' a
#' a[1,2,3]						# selection works like a matrix
#' a[1,2,]
#' a[1,,]
#' a[-1,2:3,1:2]
#' a*5								# ditto for element-wise operations
#' a <- array(dim=c(2,3,2,5,6))			# can have any number of dimensions
#' dim(a)
#' # Many packages have built-in data for testing and educational purposes
#' #data()							# list them all
#' #data(package="base")				# all base package
#' #?USArrests						# get help on a data set
#' data(USArrests)					# load the data set
#' head(USArrests)						# view the object
#' # R's workhorse is the "plot" command
#' plot(USArrests$Murder,USArrests$UrbanPop)
#' plot(USArrests$Murder,USArrests$UrbanPop,log="xy")	# log-log scale
#' plot(USArrests$Murder,USArrests$Assault,xlab="Murder",ylab="Assault",main="My Plot")
#'
#' # Can also add text
#' plot(USArrests$Murder,USArrests$Assault,xlab="Murder",ylab="Assault",main="My Plot",type="n")
#' text(USArrests$Murder,USArrests$Assault,rownames(USArrests),cex=.5)
#'
#' # Histograms and boxplots are often helpful
#' hist(USArrests$Murder)
#'
#' boxplot(USArrests)
#' ## More to come!
#' \dontrun{
#' # We won't use them right now, but here are some useful commands:
#' ?read.table							# a workhorse routine
#' ?read.csv							# a specialized CSV version
#' ?scan								# a more low-level variant
#' apropos("read")						# list various "read" commands
#' ?load								# loads objects in native R format
#' ?save								# saves objects in native R format
#' ?write.table							# counterpart to read.table
#' apropos("write")						# various "write" functions
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


#' Lab 4
#'
#' \code{lab4} returns simple description of the fourth lab in Soc 88412.
#'
#' This is a simple description of the lab, use help(lab4) to review the
#' examples for this lab
#'
#' @param Empty
#' @examples
#'
#' #
#' #-A reminder on getting started....---------------------------------------------
#' #
#' library(sna)                      # Load the sna library
#'
#' # Note: today's exercise requires the numDeriv package.  It's normally
#' # considered optional, so the system doesn't install it by default.  You can
#' # be sure that you have the package ready to roll by typing
#' #
#' #     install.packages("numDeriv")
#' #
#' # before beginning the exercise.
#' #
#' #Network autocorrelation models - the generative side---------------------------
#' #
#' # Let's begin by constructing some simulated data:
#' w1<-rgraph(100)                                        #Draw an AR matrix
#' w2<-rgraph(100)                                        #Draw an MA matrix
#' x<-matrix(rnorm(100*5),100,5)                          #Draw some covariates
#' r1<-0.2                                                #Set the model parameters
#' r2<-0.1
#' sigma<-0.1
#' beta<-rnorm(5)
#' nu<-rnorm(100,0,sigma)                                 #Draw the disturbances
#'
#' # Assemble y from its components, for several scenarios:
#' ex<-nu                                         #Draw "raw" errors
#' yx<-x%*%beta+ex                                #Compute y w/no special effects
#'
#' ex1<-nu                                        #Draw "raw" errors
#' yx1<-qr.solve(diag(100)-r1*w1,x%*%beta+ex1)    #Compute y w/AR
#'
#' ex2<-qr.solve(diag(100)-r2*w2,nu)              #Draw errors w/MA
#' yx2<-x%*%beta+ex2                              #Compute y w/no special effects
#'
#' ex12<-qr.solve(diag(100)-r2*w2,nu)             #Draw errors w/MA
#' yx12<-qr.solve(diag(100)-r1*w1,x%*%beta+ex12)  #Compute y w/AR
#'
#' #
#' #Network autocorrelation models - the inferential side--------------------------
#' #
#' # Note: some of these models take a while to fit.  This is normal, if
#' # regrettable.
#'
#' # Fit models to the "normal" case - should find that simple is best
#' fit.x.x<-lnam(yx,x=x)
#' fit.x.x1<-lnam(yx,x=x,W1=w1)
#' fit.x.x2<-lnam(yx,x=x,W2=w2)
#' fit.x.x12<-lnam(yx,x=x,W1=w1,W2=w2)
#' summary(fit.x.x)
#' cbind(fit.x.x$beta,beta)                                    #Should be close
#'
#' # Fit models to the AR case - should find that AR is best
#' fit.x1.x<-lnam(yx1,x=x)
#' fit.x1.x1<-lnam(yx1,x=x,W1=w1)
#' fit.x1.x2<-lnam(yx1,x=x,W2=w2)
#' fit.x1.x12<-lnam(yx1,x=x,W1=w1,W2=w2)                       #Can be slow
#' summary(fit.x1.x1)
#' cbind(fit.x1.x$beta,beta)                                   #Should be scary bad
#' cbind(fit.x1.x1$beta,beta)                                  #Should be close
#'
#' # Fit models to the MA case - should find that MA is best
#' fit.x2.x<-lnam(yx2,x=x)
#' fit.x2.x1<-lnam(yx2,x=x,W1=w1)
#' fit.x2.x2<-lnam(yx2,x=x,W2=w2)
#' fit.x2.x12<-lnam(yx2,x=x,W1=w1,W2=w2)
#' summary(fit.x2.x2)
#' cbind(fit.x2.x$beta,beta)                                   #Should be so-so
#' cbind(fit.x2.x2$beta,beta)                                  #Should be close
#'
#' # Fit models to the ARMA case - should find that ARMA is best
#' fit.x12.x<-lnam(yx12,x=x)
#' fit.x12.x1<-lnam(yx12,x=x,W1=w1)
#' fit.x12.x2<-lnam(yx12,x=x,W2=w2)
#' fit.x12.x12<-lnam(yx12,x=x,W1=w1,W2=w2)
#' summary(fit.x12.x12)
#' cbind(fit.x12.x$beta,beta)                                    #Should be awful
#' cbind(fit.x12.x12$beta,beta)                                  #Should be close
#'
#' # Finally, note that lnam has a plot method
#' plot(fit.x12.x12)
#'
lab4 <- function(){"This is Lab four."}


#require(RCurl)
#url<-"https://docs.google.com/spreadsheets/d/1_nuyyVE_bBX-p4Rrqa2WDi6LOEcdeZn5HETGg4HoJDY/pub?output=csv"
#myCsv <- getURL(url)
#myData<-read.csv(textConnection(myCsv))






