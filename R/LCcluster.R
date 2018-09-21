"plot.LCcluster" <-
function (x, ...) 
{
    if (x$method == "diana")
        plot(x$hclus, main = "Diana Divisive Dendrogram", 
            sub = paste("Divisive Coefficient = ", round(x$hclus$dc, 
            digits = 2)))
    else
        plot(x$hclus, main = paste( x$method, "Agglomerative Dendrogram"))
}

"print.LCcluster" <-
function (x, ...) 
{
    cat("\n\nLCcluster object: Hierarchical Clustering for LC\n")
    cat("\nData Frame input:", x$dframe)
    if( x$method == "diana" ) cat("\nClustering algorithm used: diana")
    cat("\nCovariate X variables:")
    print(x$xvars, quote = FALSE)
    if( x$method != "diana" ) print(x$hclus)
    else cat("\nDivisive Coefficient = ", round(x$hclus$dc, digits = 2), "\n\n")    
}

"LCcluster" <-
function (dframe, xvars, method = "ward.D") 
{
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("First argument to LCcluster() must be an existing data.frame name.")
    if (missing(xvars)) 
        stop("Second argument to LCcluster() must be a list of X variables.")
    # Center and Rotate X-coordinates of Experimental Units...
    xpc <- prcomp(dframe[, xvars], scale. = TRUE, rank. = length(xvars))
    # Calculate Mahalanobis Coordinates...
    for ( i in 1:length(xvars)) {
        z = xpc$x[,i] / xpc$sdev[i]
        if ( i == 1 )
            xmat = z
        else
            xmat = cbind( xmat, z )
    }
    dim(xmat) <- c(length(z),length(xvars))
    if (method == "diana") {
        hclus <- diana(dist(xmat), metric = "euclidean", 
            stand = TRUE, keep.diss = FALSE, keep.data = FALSE)
    }
    else {
        hclus <- hclust(dist(xmat), method = method)
    }
    dframe <- deparse(substitute(dframe))
    olist <- list(dframe = dframe, xvars = xvars, method = method, 
        hclus = hclus)
    class(olist) <- "LCcluster"
    olist
}
