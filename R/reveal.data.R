"reveal.data" <-
function(x, clus.var="Clus", effe.var="eSiz")
{
    # Prepare for LC Reveal analyses: Form a data.frame by appending the LTD or LRC treatment
    # effect-size measure from ltdagg() or lrcagg() as well as a Cluster membership-number
    # variable to a copy of data.frame specified in LCsetup(). 
    if (missing(x) || (!inherits(x, "ltdagg") && !inherits(x, "lrcagg")))  
        stop("First argument to reveal.data() must be a ltdagg() or lrcagg() output object.")
    if (inherits(x, "lrcagg")) {
        type = 2
	    LCdist <- x$LRCdist
    }
    else {
        type = 1
	    LCdist <- x$LTDdist
    }
    LCdist <- LCdist[order(LCdist$ID),]
    inpdf <- get(x$dframe)
    onams <- c( clus.var, effe.var, names(inpdf)) 
    outdf <- as.data.frame(cbind(LCdist[,1], LCdist[,5], inpdf))
    names(outdf) <- onams
    # class(outdf) <- "data.frame"
    outdf
}
