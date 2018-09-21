"LCsetup" <-
function (hclobj, dframe, trex, yvar) 
{
    if (missing(hclobj) || (!inherits(hclobj$hclus, "diana") && 
        !inherits(hclobj$hclus, "hclust"))) 
        stop("First argument to LCsetup must be a diana or hclust object.")	
    hclobj <- deparse(substitute(hclobj))
    if (missing(dframe) || !inherits(dframe, "data.frame")) 
        stop("Second argument to LCsetup must be an existing Data Frame.")
    if (missing(trex)) 
        stop("Third argument to LCsetup must name the Treatment or Exposure variable.")
    trex <- deparse(substitute(trex))
    if (!is.element(as.character(trex), dimnames(dframe)[[2]])) 
        stop("Treatment or Exposure must be an existing Data Frame variable.")
    tvec <- dframe[,trex]
    Kmax <- floor(length(tvec)/12)    # Guideline: Maximum Number of Clusters
    NumLevels <- length(table(tvec))
    if (NumLevels < 2)
        stop("Treatment or Exposure Level is identical for all Experimental Units.")
    if (missing(yvar)) 
        stop("Fourth argument to LCsetup must name the y-Outcome variable.")
    yvar <- deparse(substitute(yvar))
    if (!is.element(as.character(yvar), dimnames(dframe)[[2]])) 
        stop("Specified y-Outcome must be an existing Data Frame variable.")
    z <- na.omit(data.frame(cbind(dframe[,yvar], tvec)))
    names(z) <- c("y", "t")
    e <- new.env()
    e$NumLevels <- NumLevels
    e$Kmax <- Kmax
    if (NumLevels > 2) { 
        cat("\nThe Treatment variable is an Exposure with", NumLevels, "levels.")
        cat("\nLocal Treatment Difference (LTD) analyses are not applicable here.")
        cat("\nOnly Local Rank Correlations (LRCs) can be formed Within Clusters.\n\n")
        LCmean = round( cor(z$y, z$t, method = "spearman"), 8)
        e$LRCmin <- min(0, LCmean)
        e$LRCmax <- max(0, LCmean)
        }
    else { 
        cat("\nThe Treatment variable has two levels.")
        cat("\nLocal Rank Correlation (LRC) analyses are not applicable here.")
        cat("\nOnly Local Treatment Differences (LTDs) can be formed Within Clusters.\n\n")
        if (max(tvec) != 1 || min(tvec) != 0 || class(tvec) == "character") 
            stop("The two levels of the Treatment variable must be the integers: 1 or 0.")  
        n1 = sum(z$t)
        n0 = length(z$t) - n1
        if( n1 == 0 || n0 == 0 ) 
            stop("At least one Treatment Cohort contains only units with NA values.")
        LCmean = round( sum(z$y * z$t)/n1 - sum(z$y * (1-z$t))/n0, 8)
        e$LTDmin <- min(0, LCmean)
        e$LTDmax <- max(0, LCmean)		
        }
    e$aggdf <- data.frame(Label = "TEMP", Blocks = 1, LRCmean = 0, LRCstde = 0)
    e$boxdf <- data.frame(LCmean, 1)
    names(e$boxdf) <- c("LCstat", "K")
    dframe <- deparse(substitute(dframe))
    e$pars <- cbind(hclobj, dframe, trex, yvar)
    e
}
