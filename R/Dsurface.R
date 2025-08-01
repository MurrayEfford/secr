############################################################################################
## package 'secr'
## Dsurface.R
## 2011-10-21, modified through 2011-11-10
## 2012-10-25 bug fixed in getDensityArray with multi-session mask
## 2014-03-18 fixedbeta allowed
## 2014-10-13 predictD generalized for noneuc
## 2016-05-13 plot.Dsurface no longer requires 'covariate' specified for multisession input
## 2018-04-29 drop dim of predictD output
## 2020-11-05 rectangularMask rewritten to allow disjunct mask blocks
## 2022-11-28 spotHeight bug fixed xy <- data.frame(locator(1))
## 2025-07-16 sigmaxy variation of noneuc
############################################################################################

secr_predictD <- function (object, regionmask, group, session,
                se.D = FALSE, cl.D = FALSE, alpha = 0.05, 
                parameter = 'D', aslist = FALSE) {
    ## For one session and group at a time
    ## not exported; used also by region.N()
    
    if (length(parameter)==0) return (NULL)
    parameter <- match.arg(parameter, .localstuff$spatialparametersD, several.ok = TRUE)
    if (aslist) {
        ## 2025-07-23 allow multiple parameters and return list
        out <- vector(mode = 'list', length(parameter))
        names(out) <- parameter
        for (p in parameter) {
            out[[p]] <- secr_predictD(object, regionmask, group, session,
                                  se.D, cl.D, alpha, parameter = p,
                                 aslist = FALSE)
        }
        out
    }
    else {
        ## return all-1's if not relevant
        if ((parameter %in% .localstuff$spatialparameters) &&
            !(parameter %in% secr_getuserdistnames(object$details$userdist)))
            return(rep(1, nrow(regionmask)))
        sessionlevels <- session(object$capthist)
        grouplevels <- secr_group.levels(object$capthist,object$groups)
        if (is.null(session))
            session <- sessionlevels[1]
        if (is.null(group))
            groupn <- 1
        else
            groupn <- match(group, grouplevels)  # convert to numeric 2024-12-24
        if (is.null(regionmask))
            regionmask <- object$mask
        if (ms(regionmask))
            regionmask <- regionmask[[session]]
        
        ## must use mean and SD for scaling from original mask(s)
        ## this is a list for ms masks
        if (ms(object))
            meanSD <- attr(object$mask[[session]], 'meanSD')
        else
            meanSD <- attr(object$mask, 'meanSD')
        
        ## 2011-11-08 allow for 'mashing'
        if (ms(object))  {
            n.mash <- attr (object$capthist[[session]], 'n.mash')
            n.clust <- length(n.mash)
        }
        else {
            n.mash <- attr (object$capthist, 'n.mash')
            n.clust <- length(n.mash)
        }
        ## 2012-07-24 allow for unmash model fit
        unmash <- object$details$unmash
        if (is.null(unmash))
            unmash <- FALSE
        if (is.null(n.mash) | unmash)
            n.clust <- 1
        
        if (is.null(object$model$D) && parameter == 'D') {    ## implies CL && !relativeD
            ## no density model (conditional likelihood fit)
            temp <- derived(object, se.D = (se.D | cl.D)) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate'] / n.clust
            if (se.D) {
                attr(D, 'SE') <- temp['D', 'SE'] / n.clust
            }
            if (cl.D) {
                z <- abs(qnorm(1-alpha/2))
                attr(D, 'lcl') <- temp['D', 'lcl'] / n.clust
                attr(D, 'ucl') <- temp['D', 'lcl'] / n.clust
            }
            return (D)
        }
        
        ## user-defined density model
        else if (secr_userD(object) && (parameter == 'D')) {
            designD <- object$details$userDfn
            if (!is.function(designD))
                stop ("details$userDfn must be a function")
            if (se.D | cl.D)
                warning ("SE and confidence intervals may be unavailable with userDfn")
            if (n.clust>1)
                warning ("no adjustment for mashing when using userDfn")
            ## does not use link$D when calling user function userDfn
            D <- secr_getD (designD, object$fit$par, regionmask, object$parindx,
                       object$link$D, object$fixed, grouplevels, sessionlevels,
                       'D')
            ###########################################
            ## 2022-05-24 temp fix for unresolved issue
            ## return(D[,groupn,session])
            return(D[,groupn+1,session])
            ###########################################
        }
        ## linear density model on link scale
        else {
            
            newdata <- D.designdata (regionmask, object$model[[parameter]],
                                     grouplevels, sessionlevels, sessioncov = object$sessioncov,
                                     meanSD = meanSD)
            dimD <- attr(newdata, "dimD")
            ## if newdata has more than one group or session...
            if (prod(dimD[2:3]) > 1) {
                ## select a single group by number 
                groupOK <- (groupn == (1:length(grouplevels))) | (dimD[2]==1)
                groupOK <- rep(rep(groupOK, each = dimD[1]), dimD[3])
                ## select a single session
                sessionOK <- if (is.character(session))
                    session == sessionlevels
                else
                    session == 1:dimD[3]
                sessionOK <- rep(sessionOK, each = prod(dimD[1:2]))
                newdata <- newdata[sessionOK & groupOK,]
            }
            if (length(grouplevels)>1 && !('g' %in% names(newdata))) {
                # fiddle to allow groups in derivedDsurface 2024-12-24
                newdata$g <- rep(group, nrow(newdata))
                newdata$g <- factor(newdata$g, levels = grouplevels)
            }
            class(newdata) <- c('mask', 'data.frame')
            attr (newdata, 'area') <- attr(regionmask, 'area')
            
            
            #############################################
            ## allow for fixed beta parameters
            beta <- secr_complete.beta(object)
            beta.vcv <- secr_complete.beta.vcv(object)
            #############################################
            
            indx <- object$parindx[[parameter]]
            betaD <- beta[indx]
            if (object$model[[parameter]] == ~1) {
                D <- untransform(betaD, object$link[[parameter]])
                D <- max(D,0) / n.clust
                return ( rep(D, nrow(regionmask) ))
            }
            else {
                vars <- all.vars(object$model[[parameter]])
                if (any(!(vars %in% names(newdata))))
                    stop ("one or more model covariates not found")
                newdata <- as.data.frame(newdata)
                
                if (is.null(object$details[['f']]) || parameter != 'D') {
                    mat <- secr_general.model.matrix(
                        object$model[[parameter]], 
                        data = newdata, 
                        gamsmth = object$smoothsetup[[parameter]], 
                        contrasts = object$details$contrasts)
                    lpred <- as.numeric(mat %*% betaD)   ## 2018-04-29 drop matrix class, dim attribute
                    temp <- untransform(lpred, object$link[[parameter]])
                    temp <- pmax(temp, 0) / n.clust
                    
                    if (se.D | cl.D) {
                        vcv <- beta.vcv [indx,indx]
                        selpred <- sapply(1:nrow(mat), function(i)
                            mat[i,, drop=F] %*% vcv %*% t(mat[i,, drop=F]))^0.5
                        if (se.D) {
                            attr(temp, 'SE') <- se.untransform (lpred, selpred, object$link[[parameter]]) / n.clust
                            attr(temp, 'SE')[temp<=0] <- NA
                        }
                        if (cl.D) {
                            z <- abs(qnorm(1-alpha/2))
                            attr(temp, 'lcl') <- untransform (lpred - z * selpred, object$link[[parameter]]) / n.clust
                            attr(temp, 'ucl') <- untransform (lpred + z * selpred, object$link[[parameter]]) / n.clust
                            attr(temp, 'lcl')[temp<=0] <- NA
                            attr(temp, 'ucl')[temp<=0] <- NA
                        }
                    }
                }
                else {
                    f <- object$details[['f']]
                    lpred <- f(newdata[, vars[1]], beta[indx])
                    temp <- untransform(lpred, object$link[[parameter]])
                    temp <- pmax(temp, 0) / n.clust
                    if (se.D | cl.D) {
                        warning("se.D and cl.D not available for user-specified density function")
                    }
                }
                return (temp)
            }
        }
    }
}
############################################################################################

rectangularMaskold <- function (mask) {
    if (ms(mask)) {
        lapply(mask, rectangularMaskold)
    }
    else {
        temp <- expand.grid (x=sort(unique(mask$x)), y=sort(unique(mask$y)))
        temp <- read.mask(data = temp, spacing = spacing(mask))
        OK <- match(interaction(temp), interaction(mask))
        if (!is.null(covariates(mask))) {
            covariates(temp) <- covariates(mask)[OK,, drop = FALSE]
            rownames(covariates(temp)) <- 1:nrow(temp)
        }
        class(temp) <- class(mask)
        attr(temp, 'polygon') <-     attr(mask, 'polygon')
        attr(temp, 'poly.habitat') <-     attr(mask, 'poly.habitat')
        attr(temp, 'type') <- 'user'
        attr(temp, 'meanSD') <- attr(mask, 'meanSD')
        attr(temp, 'area') <-  attr(mask, 'area')
        attr(temp, 'boundingbox') <- attr(mask, 'boundingbox')
        attr(temp, 'OK') <- !is.na(OK)
        temp
    }
}
############################################################################################

rectangularMask <- function (mask) {
    if (ms(mask)) {
        lapply(mask, rectangularMask)
    }
    else {
        # rewritten 2020-11-05
        sp <- spacing(mask)
        minx <- min(mask$x)
        maxx <- max(mask$x)
        miny <- min(mask$y)
        maxy <- max(mask$y)
        nx1 <- round((maxx-minx)/sp)
        ny1 <- round((maxy-miny)/sp)
        
        # build new mask
        temp <- expand.grid (x = 0:nx1, y = 0:ny1)
        temp <- sweep(temp*sp, MARGIN = 2, STATS = c(minx,miny), FUN = "+")
        temp <- read.mask(data = temp, spacing = sp)
        
        # which new points occur in old mask?
        maski <- data.frame(x = round((mask$x-minx)/sp), y = round((mask$y-miny)/sp))
        # maski <- maski$y * (nx1+1) + maski$x
        maski <- maski$y * (nx1+1) + maski$x + 1    ## 2021-05-16
        OK <- logical(nrow(temp))
        OK[maski] <- TRUE
        
        # add covariates?
        if (!is.null(covariates(mask))) {
            covariates(temp) <- covariates(mask)[rep(1,nrow(temp)),,drop = FALSE]
            covariates(temp)[,] <- NA
            covariates(temp)[maski,] <- covariates(mask)[,]
            rownames(covariates(temp)) <- 1:nrow(temp)
        }
        
        # attributes
        class(temp) <- class(mask)
        attr(temp, 'type') <- 'user'
        attr(temp, 'meanSD') <- attr(mask, 'meanSD')   ## same as old mask? cf getMeanSD(mask)
        attr(temp, 'polygon') <- attr(mask, 'polygon')           ## same as old mask  
        attr(temp, 'poly.habitat') <- attr(mask, 'poly.habitat') ## same as old mask
        attr(temp, 'area') <-  attr(mask, 'area')                ## same as old mask
        attr(temp, 'boundingbox') <- attr(mask, 'boundingbox')   ## same as old mask
        attr(temp, 'OK') <- OK
        temp
    }
}
############################################################################################

predictDsurface <- function (object, mask = NULL, se.D = FALSE, cl.D = FALSE, alpha = 0.05,
                             parameter = 'D') {
    parameter <- match.arg(parameter, .localstuff$spatialparametersD)
    sessionlevels <- session(object$capthist)
    grouplevels <- secr_group.levels(object$capthist, object$groups)
    if (is.null(mask))
        mask <- object$mask
    densitylist <- vector('list')
    if (ms(mask)) {
        if (is.null(names(mask)))
            names(mask) <- sessionlevels
        if (any(names(mask) != sessionlevels))
            stop("names(object$mask) conflicts with session(object$capthist)")
    }
    for (session in sessionlevels) {
        if (ms(mask))
            sessmask <- mask[[session]]
        else
            sessmask <- mask
        D <- data.frame (seq = 1:nrow(sessmask))
        for (group in grouplevels) {
            predicted <- secr_predictD(object, sessmask, group, session, se.D, cl.D, 
                                  alpha, parameter)
            D[,paste(parameter,group,sep='.')] <- as.numeric(predicted)
            if (se.D) {
                D[,paste('SE',group,sep='.')] <- as.numeric(attr(predicted, 'SE'))
            }
            if (cl.D) {
                D[,paste('lcl',group,sep='.')] <- as.numeric(attr(predicted, 'lcl'))
                D[,paste('ucl',group,sep='.')] <- as.numeric(attr(predicted, 'ucl'))
            }
        }
        if (is.null(covariates(sessmask)))
            tempcov <- D
        else
            tempcov <-  cbind (covariates(sessmask), D)
        tempcov$seq <- NULL ## drop dummy

        ## impose null observations
        if (!is.null(attr(sessmask, 'OK')))
             tempcov[!attr(sessmask,'OK'),] <- NA
        if (ms(mask))
            covariates(mask[[session]]) <- tempcov
        else
            covariates(mask) <- tempcov
    }
    if (ms(mask)) {
        ## drop 'data.frame' because fools ms(), but worrying...
        class (mask) <- c('Dsurface', 'mask', 'list')
        for (i in 1:length(mask))
            class(mask[[i]]) <- c('Dsurface','mask', 'data.frame')
    }
    else
        class (mask) <- c('Dsurface', 'mask', 'data.frame')
    attr (mask, 'parameter') <- parameter
    mask
}
############################################################################################

plot.Dsurface <- function (x, covariate, group = NULL, plottype = 'shaded',
     scale = 1, ...) {
    # 2016-05-13
    if (missing(covariate)) {
        covariate <- attr(x, 'parameter')
        if (is.null(covariate)) covariate <- 'D'  ## for backwards compatibility
    }
    if (ms(x)) {
        breaklist <- lapply(x, plot, covariate, group, plottype, scale, ...)
        invisible(breaklist)
    }
    else {
        if (is.null(group))
            group <- 0
        # moved 2016-05-13
        # if (missing(covariate)) {
        #     covariate <- attr(x, 'parameter')
        #     if (is.null(covariate)) covariate <- 'D'  ## for backwards compatibility
        # }
        if (length(covariate)>1)
            stop ("whoa... just one at a time")
        if (covariate %in% c(.localstuff$spatialparametersD,'SE','lcl','ucl')) {
            covariate <- paste(covariate, group, sep='.')
        }
        if (!(covariate %in% names(covariates(x))))
            stop ("covariate ", covariate, " not found")
        covariates(x)[,covariate] <- covariates(x)[,covariate] * scale
        if (plottype %in% c('contour','persp')) {
            xval <- sort(unique(x$x))
            yval <- sort(unique(x$y))
            if (nrow(x) != length(xval)*length(yval)) {
                x <- rectangularMask(x)
                if(nrow(x) != length(xval)*length(yval))
                    stop ("failed to convert irregular mask to rectangle")
            }
            zmat <- matrix(covariates(x)[,covariate], nrow = length(xval))
            if (plottype == 'contour')
                contour(x=xval, y=yval, z=zmat, ...)
            else
                persp(x=xval, y=yval, z=zmat, ...)
        }
        else {
            class(x) <- c('mask','data.frame')
            ## use scale = 1 because already scaled by now
           
            covlevels <- plot(x, covariate = covariate, dots = (plottype == 'dots'),
                              scale = 1, ...)
            if (!is.null(covlevels)) invisible(covlevels)
        }
    }
}
############################################################################################

# re-written for noneuc 2025-07-16
Dsurface.as.data.frame <- function (x, scale = 1) {
    covnames <- names(covariates(x))
    parameter <- attr(x, 'parameter')
    if (is.null(parameter)) parameter <- 'D'    ## for backwards compatibility
    OK <- grepl(parameter, covnames) |
        (substring(covnames,1,3) == 'SE.') |
        (substring(covnames,1,4) %in% c('lcl.', 'ucl.'))
    covnames <- covnames[OK]
    densities <- covariates(x)[,covnames] * scale
    df <- cbind(x, densities)
    names(df) <- c('x','y',covnames)
    df
}
############################################################################################

print.Dsurface <- function (x, scale = 1, ...) {
    if (ms(x)) {
        out <- vector('list')
        for (session in names(x)) {
            cat ('Session ', session, '\n')
            print(x[[session]], scale, ...)
            out[[session]] <- x[[session]]
        }
        names(out) <- names(x)
        out
    }
    else {
        df <- Dsurface.as.data.frame(x, scale)
        print(df, ...)
    }
    invisible(df)
}
############################################################################################

summary.Dsurface <- function (object, scale = 1, ...) {
    if (ms(object)) {
        temp <- lapply(object, summary, scale, ...)
        class(temp) <- c('summary.Dsurface', 'list')
      temp
    }
    else {
        covnames <- names(covariates(object))
        prefix <- attr(object, 'parameter')
        if (is.null(prefix)) prefix <- 'D'    ## for backwards compatibility
        prefix <- paste(prefix, '.', sep='')        
        covnames <- covnames[substring(covnames,1,8) == prefix]
        densities <- covariates(object)[,covnames,drop=FALSE]
        densities <- densities * scale
        if (dim(densities)[2] > 1)
            densities$Total <- apply(densities,1,sum)

         class(object) <- c('mask','data.frame')
         covariates(object) <- NULL
         tempmasksummary <- summary(object,...)

        list (
            mask = tempmasksummary,
            density = apply (densities, 2, summary))
    }
}
############################################################################################

spotHeight <- function (object, prefix = NULL, dec = 2, point = FALSE, text = TRUE,
                        sep = ', ', session = 1, scale = 1, ...) {
    if ((!inherits(object, 'mask')) | (is.null(covariates(object))))
        stop ("requires plotted Dsurface or mask with covariates")
    ## Esc or click outside to break
    decxy <- 0   ## decimal places for x-y coordinates
    if (is.null(prefix)) {
      if (inherits(object,'Dsurface')) {
        prefix <- attr(object, 'parameter')
        if (is.null(prefix)) prefix <- 'D'    ## for backwards compatibility
        prefix <- paste(prefix, '.', sep='')
      }
    else
      prefix <- ''
    }
    # can only deal with one session
    if (ms(object))
        object <- object[[session]]

    covnames <- names(covariates(object))
    if (all(prefix == ""))
        OK <- rep(TRUE, times=length(covnames))
    else
        OK <- pmatch(prefix, covnames)
     if (all(is.na(OK)))
        stop("prefix does not match any covariates")
    out <- vector('list')
    i <- 0
    repeat {
        xy <- data.frame(locator(1))
        if (is.null(xy) || nrow(xy) == 0)
            break
        maskrow <- nearesttrap(xy, object)
        if (distancetotrap(xy,object[maskrow,]) > (2 * spacing(object)))
            break

        xy <- object[maskrow,,drop = FALSE]  ## centre
        if (point)
            points(xy, ...)
        D <- covariates(object)[maskrow,OK, drop = FALSE]
        for (j in 1:ncol(D)) {
            if (is.numeric(D[,j]))
                D[,j] <- format(round(D[,j]*scale, dec), nsmall = dec, trim = TRUE)
        }
        out[[i<-i+1]] <- cbind(xy,D)
        Dstr <- paste(D, collapse = sep)
        if (text)
            text(xy[1], xy[2], Dstr, ...)
        cat ('xy ', paste(round(xy,decxy), collapse = ', '), " : ", Dstr, '\n')
        flush.console()
    }
    invisible (do.call(rbind, out))
}
############################################################################################

secr_getDensityArray <- function (x, paddedlength = NULL) {
    if (ms(x)) {
        ## find maximum mask points (may vary between sessions)
        maxnrow <- max(sapply(x, nrow))
        nsession <- length(x)
        densities <- lapply(x, secr_getDensityArray, maxnrow)
        do.call(abind, densities)
    }
    else {
        covnames <- names(covariates(x))
        OK <- substring(covnames,1,2) == 'D.'
        covnames <- covnames[OK]
        densities <- covariates(x)[,covnames, drop = FALSE]
        nDcol <- ncol(densities)
        if (is.null(paddedlength))
            paddedlength <- nrow(x)
        dmat <- array(dim=c(paddedlength,nDcol))
        dmat[1:nrow(x),] <- unlist(densities)
#        if (!is.null(paddedlength))
#            densities <- c(densities, rep(NA, paddedlength-length(densities)))
#        array(densities, dim=c(nrow(x), length(covnames), 1))
        array(dmat, dim=c(dim(dmat), 1))
    }
}

############################################################################################
