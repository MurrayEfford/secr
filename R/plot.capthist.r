##############################################################################
## package 'secr'
## plot.capthist.R
## 2013-11-20
## 2015-10-11 type = 'sightings'
## 2016-10-08 secr 3.0
## 2016-12-08 type = "centres"
##############################################################################

plot.capthist <- function(x, rad = 5,
   hidetraps = FALSE, tracks = FALSE,
   title = TRUE, subtitle = TRUE,
   add = FALSE,
   varycol = TRUE, icolours = NULL, randcol = FALSE,
   lab1cap = FALSE, laboffset = 4,
   ncap = FALSE,
   splitocc = NULL, col2 = 'green',
   type = c("petal", "n.per.detector", "n.per.cluster", "sightings", "centres", "telemetry"),
   cappar = list(cex=1.3, pch=16, col='blue'),
   trkpar = list(col='blue', lwd=1),
   labpar = list(cex=0.7, col='black'),
   ...)

    # see also version in d:\single sample with stems=F, mst=F 2009 02 22

{
    ## recursive if list of capthist
    if (ms(x)) {
        if ((prod(par()$mfrow) < length(x)) & !add)
            warning("screen layout does not allow for all sessions and some plots may be lost;",
                    " set par mfrow")
        sapply (x, plot.capthist,
                rad = rad, hidetraps = hidetraps, tracks = tracks,
                title = title, subtitle = subtitle, add = add, varycol = varycol, icolours =
                    icolours, randcol = randcol, lab1cap = lab1cap, laboffset =
                    laboffset, ncap = ncap, splitocc = splitocc, col2 = col2,
                type = type, cappar = cappar, trkpar = trkpar, labpar = labpar, ...)
    }
    else {
        plotproxcapt <- function (xy, occ, icol, emphasis = FALSE) {
            oxy <- order(occ)  # sort by occasion
            xy <- xy[oxy,]
            occ <- occ[oxy]
            xy[,1] <- xy[,1] + cos(occ * 2 * pi / nocc) * rad
            xy[,2] <- xy[,2] - sin(occ * 2 * pi / nocc) * rad
            
            if (!is.null(splitocc)) {
                colr <- ifelse(occ<splitocc,cappar$col, col2)
                # trkpar$x <- xy
                # trkpar$col <- colr
                # do.call(lines, trkpar)
                # cappar$x <- xy
                # cappar$col <- colr
                # do.call(points, cappar)
                par(trkpar)
                if (tracks) lines (xy, col = colr)
                par(cappar)
                points (xy, col = colr)
            }
            else {
                # trkpar$x <- xy
                # trkpar$col <- colr
                # do.call(lines, trkpar)
                # cappar$x <- xy
                # cappar$col <- colr
                # do.call(points, cappar)
                par(trkpar)
                if (varycol) par(col=icol)   
                if (tracks) lines (xy)
                par(cappar)
                if (varycol) par(col=icol)   
                points (xy)
                if (emphasis)
                    points (xy, col='black', bg=par()$col, pch=21)
            }
        }
        plotpolygoncapt <- function (xy, icol, emphasis = FALSE) {
            par(trkpar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            if (tracks) lines (xy)
            par(cappar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            points (xy)
            if (emphasis)
                points (xy, col='black', bg=par()$col, pch=21)
        }
        labcapt <- function (n) {
            if ( detectr[1] %in% c('proximity', 'count', 'polygonX',
                'transectX', 'signal', 'signalnoise', 'polygon',
                'transect', 'unmarked', 'presence') ) {
                warning ("labels not implemented for this detector type")
            }
            else {
                xn <- apply(abs(x[n,,,drop=FALSE]),2,sum)
                o1 <- sum(cumsum(xn)==0)+1   # first occasion
                t1 <- which(x[n,o1,]>0)                 # first trap site
                # cat(n, ' ',xn, ' o1 ', o1, ' t1 ', t1, '\n')
                dx  <- (cos((1:nocc) * 2 * pi / nocc) * rad)[o1]
                dy  <- (sin((1:nocc) * 2 * pi / nocc) * rad)[o1]
                par(labpar)
                if (varycol) par(col=n)   # override
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                text (traps$x[t1]+dx+laboffset[1], traps$y[t1]-dy+laboffsety, row.names(x)[n])
            }
        }
        labhead <- function (n, df) {
                par(labpar)
                if (varycol) par(col=n)   # override
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                text (head(df[[n]],1)$x++laboffset[1], head(df[[n]],1)$y+laboffsety,
                      row.names(x)[n])
        }
        ncapt <- function (x) {
            if (detectr[1] %in% .localstuff$polydetectors) {
               stop ("ncap does not work with polygon and similar detectors")
            }
            temp <- t(apply (abs(x),c(2,3),sum)) # capts/trap/day)
                
            dx  <- rep(cos((1:nocc) * 2 * pi / nocc) * rad, rep(nrow(traps),nocc))
            dy  <- rep(sin((1:nocc) * 2 * pi / nocc) * rad, rep(nrow(traps),nocc))
            par(labpar)
            par(adj=0.5)
            OK <- temp>0
            text ((traps$x[row(temp)]+dx)[OK], (traps$y[row(temp)]-dy)[OK], as.character(temp[OK]))
        }

        plotsignal <- function (df, minsignal, maxsignal,n) {
            # df is a dataframe for one animal
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day

            ## function (signal, occasion, trap, minsignal, maxsignal, n)
            .localstuff$i <- .localstuff$i+1
            dx <- rep( (cos((.localstuff$i) * 2 * pi / n) * rad), nrow(df))
            dy <- rep( (sin((.localstuff$i) * 2 * pi / n) * rad), nrow(df))
            sq <- order(df$signal)     # plot darkest points last
            df <- df[sq,]
            df$trap <- as.character(df$trap)
            if (maxsignal>minsignal)
                greycol <- grey(0.7 * (1 - (df$signal-minsignal)/(maxsignal-minsignal)))
            else
                greycol <- grey(0.5)
            if (tracks) lines (traps$x[df$trap]+dx, traps$y[df$trap]-dy, col=greycol)
            par(cappar)
            points (traps[df$trap,'x']+dx, traps[df$trap,'y']-dy, col = greycol)
        }
        plotsightings <- function (x) {
            Tu <- Tu(x)
            Tu0 <- Tu
            marking <- Tu
            if (is.null(Tu)) stop ("sightings type requires sighting data Tu")
            markocc <- markocc(traps(x))
            Tu0[Tu!=0] <- NA
            Tu0[,markocc>0] <- NA
            Tu[Tu==0] <- NA
            marking[,markocc<1] <- NA
            dx  <- rep((cos((1:nocc) * 2 * pi / nocc) * rad), each = nrow(Tu))
            dy  <- rep((sin((1:nocc) * 2 * pi / nocc) * rad), each = nrow(Tu))
            dx0 <- dx; dx0[is.na(Tu0)] <- NA
            
            if (all(detector(traps(x)) %in% c('polygon'))) {
                centres <- split(traps(x), polyID(traps(x)))
                ## assume each polygon closed, so first vertex redundant
                centres <- lapply(centres, function(xy) apply(xy[-1,,drop=FALSE], 2, mean))
                trapxy <- data.frame(do.call(rbind,centres))
                names(trapxy) <- c('x','y')
            }
            else trapxy <- traps(x)
            
            par(labpar)
            text (rep(trapxy$x, nocc) + dx, rep(trapxy$y, nocc) - dy, Tu, cex = 0.8)
            par(cappar)
            points (rep(trapxy$x, nocc) + dx0, rep(trapxy$y, nocc) - dy, pch = 1, cex = 0.7)
            ## marking occasions shown as dot
            dx[is.na(marking)] <- NA
            points (rep(trapxy$x, nocc) + dx, rep(trapxy$y, nocc) - dy, pch=16, cex=0.4)
        }
        
        plotcentres <- function (x) {
            xtraps <- traps(x)
            meanxya <- function(xy) apply(xy, 2, mean)
            meanxy <- function(trp) apply(xtraps[trp,],2,mean)
            if (all(detector(traps(x)) %in% 'telemetry')) {
                xyl <- telemetryxy(x)
                xyl <- lapply(xyl, meanxya)
                xy <- do.call(rbind, xyl)
                
            }
            else if (!any(detector(traps(x)) %in% .localstuff$polydetectors)) {
                ## 2021-05-19 sortorder not material
                trplist <- split(trap(x, sortorder = 'ksn'), animalID(x, sortorder = 'ksn'))
                xyl <- lapply(trplist, meanxy)
                xy <- do.call(rbind, xyl)
            }
            else {
                xyl <- telemetryxy(x)
                if (is.null(xyl))
                    xyl <- split(xy(x), animalID(x, sortorder = "ksn"))
                xyl <- lapply(xyl, meanxya)
                xy <- do.call(rbind, xyl)
            }
            if (rad>0) xy <- xy + (runif(2*nrow(xy))-0.5) * rad
            points (xy, ...)
        }
        
        ###########
        ## MAINLINE
        x <- check3D(x)
        opal <- palette() ; on.exit(palette(opal))
        type <- match.arg(type)
        traps <- traps(x)
        detectr <- expanddet(x)
        ## if (is.null(rownames(x))) { 2020-09-09
        if (is.null(rownames(x)) && nrow(x)>0) {
                warning ("capthist has no rownames; using 1:nrow")
            rownames(x) <- 1:nrow(x)
        }
        if (all(detectr == 'telemetry')) {
            type <- 'telemetry'
            # warning("assuming type = 'telemetry' as all data from telemetry")
        }
        if (type =='telemetry')
            nocc <- sum(detectr == 'telemetry')
        else
            nocc <- sum(detectr != 'telemetry')
        
        if (type %in% c('petal','centres'))
            cappar <- replacedefaults (list(cex=1.3, pch=16, col='blue'), cappar)
        if (type == 'sightings')
            cappar <- replacedefaults (list(cex=1, pch=16, col='blue'), cappar)
        if (type %in% c('n.per.cluster','n.per.detector'))
            cappar <- replacedefaults (list(cex = 3, pch = 21), cappar)

        trkpar <- replacedefaults (list(col='blue', lwd=1), trkpar)
        labpar <- replacedefaults (list(cex=0.7, col='black'), labpar)
        initialpar <- par(cappar)
        if (!add) {
            if (type=="telemetry") {
                xyl <- telemetryxy(x)
                xy <- do.call(rbind, xyl)
                xl <- range(xy[,1])
                yl <- range(xy[,2])
                tr <- expand.grid(x=xl, y=yl)
                class(tr) <- c('traps', 'data.frame')
                plot(tr, hidetr=TRUE, ...)
            }
            else {
                plot(traps, hidetr=hidetraps, ...)
            }
        }

        if (nrow(x) == 0) {
            warning("no detections in capthist object")
            type <- 'null'
        }

        if (is.null(icolours)) icolours <- topo.colors((nrow(x)+1)*1.5)
        if (varycol) {
            if (randcol) icolours <- sample(icolours)
            test <- try (palette(icolours))  ## too many?
            if (inherits(test, 'try-error'))
                stop ("requested too many colours; ",
                      "try with varycol = FALSE")
            icol <- 0
        }
        else {
            # splitcol?
        }
        
        if (type == 'petal') {

            if ((nocc == 1) & ! (detectr[1] %in% c('signal','signalnoise'))) rad <- 0

            if ( detectr[1] %in% .localstuff$polydetectors ) {
                ## occasions not distinguished
                lxy <- split (xy(x), animalID(x, names = FALSE, sortorder = "ksn"))
                mapply (plotpolygoncapt, lxy, 1:length(lxy))
            }
            else if ( detectr[1] %in% c('signal','signalnoise') )
            {
                .localstuff$i <- 0
                temp <- data.frame(
                    ID = animalID(x, sortorder = "ksn"), 
                    occ = occasion(x, sortorder = "ksn"),
                    trap = trap(x, sortorder = "ksn"), 
                    signal = signal(x))
                lsignal <- split(temp, animalID(x, names = FALSE, sortorder = "ksn"))
                lapply(lsignal, plotsignal, minsignal = min(temp$signal),
                    maxsignal = max(temp$signal), n=nrow(x))
            }
            else  { 
                xydf <- as.data.frame(traps(x)[trap(x, sortorder = 'snk'),])
                occ <- occasion(x, sortorder = 'snk')
                OK <- detectr[occ] != 'telemetry'
                ID <- factor(animalID(x), levels = rownames(x))
                lxy <- split(xydf[OK,], ID[OK])
                occ <- split(occ[OK], ID[OK])
                if (tracks & any(unlist(sapply(occ, duplicated))))
                    warning("track for repeat detections on same occasion",
                            " joins points in arbitrary sequence")
                mapply(plotproxcapt, lxy, occ, 1:length(lxy), telemetered(x))
            }

            if (lab1cap) {
                if ( detectr[1] %in% .localstuff$polydetectors ) {
                    lxy <- split (xy(x), animalID(x, names = FALSE, sortorder = "ksn"))
                    sapply(1:nrow(x), labhead, df=lxy)
                }
                else {
                    sapply(1:nrow(x), labcapt)
                }
            }

            if (ncap) { ncapt(x)}
        }
        else if (type %in% c('n.per.cluster','n.per.detector')) {
            if (type == 'n.per.detector') {
                ## never yields zeros
                temp <- table(trap(x, sortorder = 'snk'), animalID(x, sortorder = 'snk'))>0
                nj <- apply(temp,1,sum)
                centres <- traps(x)[names(nj),]
            }
            else if (type == 'n.per.cluster') {
                nj <- cluster.counts(x)
                centres <- cluster.centres(traps)
                # hide zeros, if present
                centres <- centres[nj>0,]
                nj <- nj[nj>0]
            }
            else stop ("unrecognised type")

            if (is.null(icolours)) {
                icolours <- topo.colors(max(nj)*1.5)
            }
            palette (icolours)
            npal <- length(icolours)
            if (max(nj) < npal) {
                cols <- npal-nj
            }
            else {
                cols <- round(npal * (1-nj/max(nj)))
            }
            if (cappar$pch == 21)
                fg <- 'black'
            else
                fg <- cols
            par(cappar)
            points(centres, col = fg, bg = cols, pch = cappar$pch, cex = cappar$cex)

            if (ncap) {
                par(labpar)
                par(adj=0.5)
                vadj <- diff(par()$usr[3:4])/500  ## better centring!
                text(centres$x, centres$y + vadj, nj)
            }

            ## should export data for legend:
            ## count classes 0, 1-, 2-,...,-max
            ## and corresponding colours
            tempcol <- npal- (1:max(nj))
            output <- data.frame(
                legend = 1:max(nj),
                col = tempcol,
                colour = icolours[tempcol],
                stringsAsFactors = FALSE
            )
        }
        else if (type == 'sightings') {
            plotsightings(x)
        }
        else if (type == 'centres') {
            plotcentres(x)
        }
        else if (type == 'telemetry') {
            lxy <- telemetryxy(x)
            seqnum <- match(names(lxy), rownames(x))
            caught <- apply(abs(x[seqnum,detectr!='telemetry',,drop=FALSE]),1,sum)>0
            mapply (plotpolygoncapt, lxy, seqnum, caught)
        }
        else
            if (type != 'null') stop ("type not recognised")
        

        ####################################################
        ## Titles

        if (type == 'telemetry') {
            nocc <- sum(detectr == 'telemetry')
            nd <- sum(abs(x)[,detectr == 'telemetry',])
            nanimal <- length(telemetryxy(x))
        }
        else {
            nocc <- sum(detectr != 'telemetry')
            nd <- sum(abs(x)[,detectr != 'telemetry',])
            nanimal <- sum(apply(abs(x)[,detectr != 'telemetry',,drop=FALSE],1,sum)>0)
        }
        pl <- if (nocc>1) 's' else ''  ## plural marker
        if (type == 'sightings') {
            markocc <- markocc(traps(x))
            Tu <- Tu(x)
            nd <- sum(Tu)
            nocc <- sum(markocc<1)
        }
        if (type == 'telemetry') {
            nd <- sum(sapply(lxy,nrow))
            nocc <- sum(detector(traps(x))=="telemetry")
            nanimal <- length(lxy)
        }
        if (is.logical(title)) {
            txt <- ifelse (is.null(session(x)), paste(deparse(substitute(x)),
                       collapse=''), session(x))
            title <- ifelse(title, txt, '')
        }
        if (title != '') {
            par(col='black')
            mtext(side=3,line=1.2, text = title, cex=0.7)
        }
        if (is.logical(subtitle)) {

            subtitle <- if (subtitle) {
                if (type == 'sightings') {
                    if (any(markocc<0))
                        paste0(nocc, ' sighting occasion', pl, ', ', nd, ' sightings of marked and unmarked animals')
                    else
                        paste0(nocc, ' sighting occasion', pl, ', ', nd, ' sightings of unmarked animals')
                }
                else if (type == 'telemetry')
                    paste0(nocc, ' occasion', pl, ', ', nd, ' fixes,',
                           nanimal, ' animals')
                else
                    paste0(nocc, ' occasion', pl, ', ', nd, ' detections, ',
                           nanimal, ' animals')
            }
                else ''
        }
        if (subtitle != '') {
            par(col='black')
            mtext(text = subtitle, side=3,line=0.2, cex=0.7)
        }
        ####################################################

        par(initialpar)   # restore
        if (type %in% c('n.per.detector','n.per.cluster'))
            invisible(output)
        else
            invisible(nd)
    }
}
############################################################################################

plotMCP <- function(x, add = FALSE, col = 'black', fill = NA, lab1cap = FALSE,
                    laboffset = 4, ncap = FALSE, ...) {
    plotone <- function (df, col, fill) {
        if (nrow(df)>0) {
            par(fg=col)
            polygon(df, col=fill)
        }
    }
    mcp <- function (df) {
        if (nrow(df)>1) {
        df <- df[chull(df[,1], df[,2]),]
        rbind(df, df[1,])
        }
        else df
    }
    if (ms(x)) {
        lapply(x, plotMCP, add = add, col = col, ...)
    }
    else {
        xyl <- telemetryxy(x)
        if (is.null(xyl)) {
            if (all(detector(traps(x)) %in% c('polygon','polygonX')))
                xyl <- split(xy(x), animalID(x, sortorder = 'ksn'))
            else {
                df <- as.data.frame(traps(x)[trap(x, names = FALSE, sortorder = 'snk'),])
                xyl <- split(df, animalID(x, sortorder = 'snk'))
            }
        }
        if (!add)
            plot(traps(x), ...)
        fg <- par()$fg
        on.exit(par(fg=fg))
        if (missing(col)) col <- 'black'
        xymcp <- lapply(xyl, mcp)
        if (length(xymcp)>0) {
        mapply(plotone, xymcp, col, fill)

        if (lab1cap | ncap) {
            labhead <- function(xy, nam, num, col) {
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                if (ncap)
                    text (xy[1,1] + laboffset[1], xy[1,2] + laboffsety, num, col=col)
                else
                    text (xy[1,1] + laboffset[1], xy[1,2] + laboffsety, nam, col=col)
            }
            mapply(labhead, xymcp, names(xyl), sapply(xyl,nrow), col)
        }
        }
        invisible(xyl)
    }
}
############################################################################################
# plotMCP(AB2004,col=1:12, gridl=F)

## 2016-10-17, 2016-11-14
occasionKey <- function (capthist, noccasions, rad = 3, x, y, px = 0.9, py = 0.9, 
                         title = 'Occasion', ...) {
    if (missing(x)) {
        ux <- par()$usr[1:2]
        x <- ux[1] + px * diff(ux)
    }
    if (missing(y)) {
        uy <- par()$usr[3:4]
        y <- uy[1] + py * diff(uy)
    }
    if (!missing(capthist)) {
        if (ms(capthist))
            stop ("occasionKey requires single-session capthist")
        noccasions <- ncol(capthist)
    }
    else {
        if (missing (noccasions))
            stop ("occasionKey requires one of capthist or noccasions")
    }
    dx  <- cos((1:noccasions) * 2 * pi / noccasions) * rad
    dy  <- sin((1:noccasions) * 2 * pi / noccasions) * rad
    points (x,y, cex = 0.5)
    text (x, y+rad*2.3, title, ...)
    text (x+dx, y-dy, 1:noccasions, ...)
}
# plot(captdata, border = 30)
# occasionKey(nocc=5, rad = 10, cex = 0.8)
############################################################################################

