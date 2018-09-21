"plot.KSperm" <- 
function(x, ...)
{
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    Dvec <- x$Dvec
    obsD <- x$obsD
    xmax = max(max(Dvec), obsD)
    par(mfrow = c(1,1))
    plot(ecdf(Dvec), verticals=TRUE, do.points=FALSE, ann=FALSE,
        col="blue2", lwd=2, xlim = c(0, xmax))
    abline(v=0, lty="solid", col="black")
    abline(v=obsD, lty="dashed", lwd=2, col="red")
    title(main = "LC Confirm Inference: Ignorable X-covariates?", 
        ylab = "Cumulative Probability",
        sub = paste("Observed D =", round(obsD, 4), ", Pmax =", x$Pmax))
    if (x$Type == 1)
        title(xlab = "Kolmogorov-Smirnov D-statistics for NULL LTDs")
    else
        title(xlab = "Kolmogorov-Smirnov D-statistics for NULL LRCs")
}

"print.KSperm" <-
function(x, ...)
{
    cat("\nKSperm Object: NULL Distribution of the Kolmogorov-Smirnov")
    cat("\nD-statistic when given X-covariates are ignorable.\n")
    cat("\nData Frame:", x$dframe, "\n")
    cat("Outcome Variable:", x$yvar, "\n")
    cat("Treatment Variable:", x$trtm, "\n")
    if (x$Type == 1)
        cat("Effect-Size estimates: Local Treatment Differences (LTDs)\n")
    else
        cat("Effect-Size estimates: Local Rank Correlations (LRCs)\n")
    cat("Number of Clusters per Replication:", x$nclus, "\n")
    cat("Number of Replications:", x$reps, "\n")
    # cat("Number of Random NULL D-statistics:", length(x$Dvec), "\n" )
    cat("\n    Observed Kolmogorov-Smirnov D-statistic =", x$obsD)
    cat("\n    Sorted NULL D-statistic values =\n\n")
    print(x$Dvec)
    cat("\n    The simulated p.value is thus less than Pmax =", x$Pmax, "\n\n")
}

"KSperm" <-
function(x, reps=100)
{
    # Compute Random Permutation Distribution of Kolmogorov-Smirnov D-statistics.
    if (missing(x) || !inherits(x, "confirm"))  
        stop("First argument to KSperm() must be a confirm() object.")
    type = x$Type
    nclus = x$nclus
    LCdist <- x$LCdist  # needed Cluster, Outcome & Treatment Info; observed lstat NOT used here
    names(LCdist) <- c("c", "y", "t", "lstat")
    units = length(LCdist[,4])
    dfconf <- x$dfconf  # saved Permutation Distribution of lstat = LTDs or LRCs
    obsD <- x$KSobsD$statistic  
    olist <- list(hiclus = x$hiclus, dframe = x$dframe, trtm = x$trtm, yvar = x$yvar,
        Type=type, reps=reps, nclus=nclus, units=units)	
    Dvec <- as.vector(rnorm(reps))
    for(i in 1:reps) { 
        cperm <- LCdist$c[order(as.vector(rnorm(units)))]	# Permuted Cluster IDs
        pdf <- as.data.frame(cbind(cperm, LCdist[,3:4]))    # y & t columns NOT Permuted
        names(pdf) <- c("c","y","t")
        if (type == 2) {
            dfLRC <- do.call( rbind, lapply( split(pdf, pdf$c),
    	        function(x) {
    	            x <- na.omit(x)
                    if (length(x[,1]) < 3)
                        LRC = NA
                    else
                        LRC = round( cor(x$y, x$t, method = "spearman"), 8)
                    data.frame(c=x$c[1], LRC=LRC, w=length(x$t))
                } )
            )
            pdf = merge(pdf, dfLRC[,1:2], by.x="c", by.y="c")
        }
        else {
            dfLTD <- do.call( rbind, lapply( split(pdf, pdf$c),
                function(x) {
    		        n1 = sum(x$t)
    		        n0 = length(x$t) - n1
    		        if(n1 == 0 || n0 == 0)
                        LTD = NA
                    else
                        LTD = round( sum(as.numeric(x$y * x$t))/n1 -
                                sum(as.numeric(x$y * (1-x$t)))/n0, 8)
                    data.frame(c=x$c[1], LTD=LTD, w=length(x$t))
                } )
            )
            pdf = merge(pdf, dfLTD[,1:2], by.x="c", by.y="c")
        }
        names(pdf) <- c("c", "y", "t", "lstat")
        pdf <- as.data.frame(pdf$lstat)
        names(pdf) <- "lstat"		
	    kso <- ks.test(pdf$lstat, dfconf$lstat)	  
	    Dvec[i] <- kso$statistic			
    }
    Dvec <- Dvec[order(Dvec)]  # ...rep samples of NULL D-distribution for ignorable X-covariates
    rnk = (rank(c(obsD, Dvec)))[1]
    maxp = round(((reps^2 + reps - 1 - (reps - 1)*rnk) / reps^2 ), 4)
    olist <- c(olist, list(obsD=obsD, Dvec=Dvec, Pmax=maxp))	
    class(olist) <- "KSperm"
    olist
}
