#############################################################################
## package 'secr'
## read.mask.R
## 2022-11-15 separate file
#############################################################################

read.mask <- function (file = NULL, data = NULL, spacing = NULL, columns = NULL, ...)
{
    if (is.null(data) & !is.null(file)) {
        fl <- nchar(file)
        SS <- tolower(substring(file,fl-3,fl)) == '.csv'
        if (SS) {
            data <- read.csv (file)
            if ('HABITAT' %in% names(data))
                data <- data[data$HABITAT == 1,]
        }
        else {
            data <- read.table (file, ...)
        }
    }
    else if (is.null(data))
        stop("require one of 'file' or 'data'")
    if (length(dim(data))!=2)
        stop ("require dataframe or matrix for 'data' input to read.mask")
    
    coln <- colnames(data)
    ixy <- match(c('x','y'), coln)
    if (any(is.na(ixy))) ixy <- 1:2
    mask <- as.data.frame(data[,ixy])
    names(mask) <- c('x', 'y')
    if (any(!apply(mask, 2, is.numeric)))
        stop ("non-numeric x or y coordinates")
    if (any(is.na(mask)))
        stop ("missing value(s) in x or y coordinates")
    
    class(mask) <- c('mask', 'data.frame')
    
    ## add covariates
    if (ncol(data) > 2) {
        df <- as.data.frame(data[,-ixy, drop = FALSE])
        if (!is.null(columns)) {
            ##            if (!all(columns %in% names(mask)))
            ## bug fixed 2014-05-31
            if (!all(columns %in% names(df)))
                stop ("columns missing from input")
            df <- df[,columns, drop=FALSE]
        }
        if (ncol(df)>0)
            covariates(mask) <- df
    }
    if (is.null(spacing))
    {
        sp      <- as.matrix(dist(as.matrix(mask)))
        spacing <- apply(sp,1,function(x) min(x[x>0]))
        ## 2020-08-25 bug fixed
        ## spacing <- mean (spacing, na.rm=T)
        spacing <- median (spacing, na.rm=T)
    }
    
    area    <- spacing^2 / 10000
    
    xl <- range(mask$x) + spacing/2 * c(-1,1)
    yl <- range(mask$y) + spacing/2 * c(-1,1)
    
    attr(mask,'type')    <- 'user'
    attr(mask,'meanSD')  <- getMeanSD(mask)
    attr(mask,'area')    <- area
    attr(mask,'spacing') <- spacing
    attr(mask,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    attr(mask,'polygon') <- NULL
    mask
}
