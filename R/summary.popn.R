#############################################################################
## package 'secr'
## summary.popn.R
## 2022-08-23 moved from methods.R 
#############################################################################

summary.popn <- function(object, collapse = FALSE, ...) {
    ## 2019-06-10
    if (ms(object)) {
        if (collapse) {
            pop1 <- do.call(rbind, object)
            temp <- summary(pop1)
            temp$collapsed <- TRUE
            # if (!is.null(covariates(pop1))) {
            #     temp$covar <- summary(covariates(pop1), ...)
            # } else 
            temp$covar <- NULL
            pop2 <- popIDsplit(object)   # dim = c(length(ID), nsess, 2); see plot.popn.R
            move <- function (sessxy) sqrt(diff(sessxy[,1])^2 + diff(sessxy[,2])^2)
            m <- t(apply(pop2,1,move))
            temp$nanimals <- nrow(pop2)
            temp$nsessions <- ncol(pop2)
            pop2 <- pop2[,,1]
            pop2[!is.na(pop2)] <- 1 
            pop2[is.na(pop2)] <- 0
            pop2 <- t(apply(pop2, 1, function(x) {x[x==0 & cumsum(x)>0] <- -1; x}))
            # 0 not yet recruited; 1 alive; -1 dead
            recruited <- function(x) c(x[1], diff(x) == 1)    # 0 -> 1
            died <- function(x) c(0, diff(x) == -2)           # 1 -> -1
            temp$recruits <- apply(apply(pop2,1,recruited),1,sum)
            temp$deaths <- apply(apply(pop2,1,died),1,sum)
            temp$movements <- m
            temp$status <- pop2
            temp
        }
        else {
            temp <- lapply(object, summary.popn)
            class(temp) <- c('summary.popn', 'list')
        }
        temp
    }
    else {
        if (is.null(object$x) | is.null(object$y))
            stop ("not a valid popn")
        nd <- length(object$x)
        if (length(object$x) != length(object$y))
            stop  ("not a valid popn")
        
        if (!is.null(covariates(object))) {
            sumcovar <- summary(covariates(object), ...)
        } else sumcovar <- NULL
        
        popnclass <- 'popn'
        
        ## rearranged 2014-09-06
        temp <- list (
            popnclass = popnclass,
            nanimals = nrow(object),
            xrange = range(object$x),
            yrange = range(object$y),
            boundingbox = attr(object, 'boundingbox',exact = TRUE),
            covar = sumcovar,
            collapsed = FALSE
        )
        class(temp) <- 'summary.popn'
        temp
    }
    
}
############################################################################################

print.summary.popn <- function (x, ...) {
    if (ms(x)) {
        lapply (x, print.summary.popn)
    }
    else {
        cat ('Object class        ', x$popnclass, '\n')
        cat ('Number of animals   ', x$nanimals, '\n')
        if (x$collapsed) {
            cat ('Sessions            ', x$nsessions, '\n')
            cat ('Animals  by session ', paste(apply(x$status==1,2,sum,na.rm=TRUE), collapse=', '), '\n')
            cat ('Recruits by session ', paste(x$recruits, collapse = ', '), '\n')
            cat ('Deaths   by session ', paste(x$deaths, collapse = ', '), '\n')
            cat ('Average move        ', round(mean(as.numeric(x$movements), na.rm = TRUE),2), '\n')
            cat ('\n')
        }      
        cat ('x-range m        ', x$xrange, '\n')
        cat ('y-range m        ', x$yrange, '\n')
        cat ('Bounding box     ','\n')
        print (x$boundingbox, ...)
        if (!is.null(x$covar)) {
            cat ('Summary of covariates', '\n')
            print(x$covar, ...)
        }
        cat('\n')
    }
}
############################################################################################

