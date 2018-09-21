"plot.ivadj" <-
function (x, maxsiz = 0.15, ...) 
{
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(lty = 1, col = "black")
    plot(x$LCtable$PS, x$LCtable$LAO, ann = FALSE, type = "n")
    symbols(x$LCtable$PS, x$LCtable$LAO, circles = sqrt(x$LCtable$w), inches = maxsiz, add = TRUE)
    abline(x$ivfit, lty = 2, lwd=2, col = "red")	
    lines(x$smfit, lty = 2, lwd=2, col = "blue2")
	if (x$Type == 1)
        title(main = paste("IV Adjustment for", length(x$LCtable$LAO), "Clusters"), 
            xlab = "Within-Cluster Treatment Fraction (Propensity)", 
            ylab = "Observed LAO", sub = "Symbol Area proportional to Cluster Size")
    else
        title(main = paste("IV Adjustment for", length(x$LCtable$LAO), "Clusters"), 
            xlab = "Within-Cluster Relative Exposure (Propensity)", 
            ylab = "Observed LAO", sub = "Symbol Area proportional to Cluster Size")			
}

"print.ivadj" <-
function (x, ...) 
{
    cat("\nivadj: Instrumental Variable (IV) Adjustment via Clustering\n")
    cat("\nCluster Tabulation:\n===================\n\n") 
    print(x$LCtable)
    cat("\nSummary of smooth.spline() fit:\n===============================\n\n") 
    print(x$smfit)
    cat("\nSummary of lm() fit:\n====================\n") 	
    print(x$ivsum)
    cat("\nListing of lm() predictions:\n============================\n\n") 	
    print(x$ivpred)
}

"ivadj" <-
function (x) 
{
    if (missing(x) || (!inherits(x, "ltdagg") && !inherits(x, "lrcagg")))  
        stop("First argument to ivadj() must be a ltdagg() or lrcagg() object.")
    if (inherits(x, "lrcagg")) {
        type = 2
	    LCtable <- x$LRCtabl
        expmin <- min(LCtable$PS)
        expmax <- max(LCtable$PS)
        LCtable$PS <- (LCtable$PS - expmin)/(expmax - expmin) # Relative Exposure PS
    } else {
        type = 1
        LCtable = x$LTDtabl
	}
    if (length(LCtable$LAO) < 5) {
        cat("\nIV Adjustment should not be attempted with fewer than 5 Clusters.\n\n")
        return(NULL)
    }
    ivfit <- lm(LCtable$LAO ~ LCtable$PS, weights = LCtable$w)
    ivsum <- summary(ivfit)
    ivpred <- predict(ivfit, data.frame(LCtable$PS), se.fit = TRUE, type="response")
    smfit <- smooth.spline(LCtable$PS, y=LCtable$LAO, w=LCtable$w)
    olist <- list(LCtable = LCtable, Type=type, ivfit=ivfit, ivsum=ivsum,
        ivpred=ivpred, smfit=smfit)
    class(olist) <- "ivadj"
    olist
}
