############################################################################################
## package 'secr'
## make.marking.R
## 2026-11-22 DRAFT: capthist for telemetrytype 'marking' in one step (see make.marking.Rd)
############################################################################################

## The arguments are in three groups: the cameras (traps); the collared sample, one row per
## collared animal, all keyed by ID (detections, telemetry, marked, covariates); and the
## sightings of animals not individually identified, counts by detector and occasion (Tu, Tm).

make.marking <- function (traps, detections = NULL, telemetry = NULL, marked = NULL,
                          covariates = NULL, Tu = NULL, Tm = NULL, noccasions = NULL,
                          verify = TRUE, ...) {

    ## ---- detections as the 'captures' of make.capthist (Session, ID, Occasion, Detector) ----
    if (!is.null(detections)) {
        if (!is.data.frame(detections) || ncol(detections) < 4)
            stop ("detections should be a data frame with columns Session, ID, Occasion, ",
                  "Detector (as 'captures' of make.capthist)")
        detections <- detections[, 1:4]
        names(detections) <- c("Session", "ID", "Occasion", "Detector")
    }

    ## ---- sessions ----
    MS <- ms(traps)
    trapslist <- if (MS) unclass(traps) else list(traps)
    nsess <- length(trapslist)
    if (MS) {
        sessnames <- names(trapslist)
        if (is.null(sessnames) || any(sessnames == ""))
            stop ("multi-session traps must be named (the session names)")
    }
    else {
        s1 <- if (!is.null(detections) && nrow(detections) > 0)
                  unique(as.character(detections$Session))
              else if (is.data.frame(telemetry) && "Session" %in% names(telemetry))
                  unique(as.character(telemetry$Session))
              else "1"
        if (length(s1) != 1)
            stop ("traps is for a single session but the data have Session values: ",
                  paste(s1, collapse = ", "))
        sessnames <- s1
    }
    for (i in seq_len(nsess))
        if (!inherits(trapslist[[i]], "traps"))
            stop ("traps must be a 'traps' object or a list of them")
    ## detector type is taken from traps: 'count' (the model now in secr) or 'proximity'
    ## (binary detection; the pooled-binary model, so far only in RTMB); one type for all
    ## sessions, and no telemetry detector (it is added here)
    dettype <- unique(unlist(lapply(trapslist, function (tr) as.character(detector(tr)))))
    if (length(dettype) != 1 || !(dettype %in% c("count", "proximity")))
        stop ("traps should have detector type 'count' or 'proximity' in every session ",
              "(no telemetry detector, no other type); found: ",
              paste(dettype, collapse = ", "))
    if (dettype == "proximity")
        message ("proximity (binary) detectors: this capthist is for the pooled-binary model; ",
                 "secr.fit() does not yet fit telemetrytype 'marking' to proximity detectors")

    ## ---- split an argument by session: NULL, one object (single session), or a list with
    ## one element per session ----
    bysession <- function (x, what) {
        if (is.null(x)) return(rep(list(NULL), nsess))
        if (is.data.frame(x) || is.matrix(x)) {                 # a single object
            if (nsess == 1) return(list(x))
            stop (what, " should be a list with one element per session ",
                  "(or a data frame with a Session column)")
        }
        if (is.list(x)) {
            if (length(x) != nsess)
                stop (what, " should be a list with one element per session")
            if (!is.null(names(x)) && all(names(x) %in% sessnames)) x <- x[sessnames]
            return(unname(x))
        }
        stop (what, " should be a data frame, matrix or list")
    }
    ## a data frame with a Session column is split by it
    splitframe <- function (x, what) {
        if (is.data.frame(x) && "Session" %in% names(x)) {
            sess <- as.character(x$Session)
            if (!all(sess %in% sessnames))
                stop ("Session in ", what, " not matched in sessions: ",
                      paste(setdiff(unique(sess), sessnames), collapse = ", "))
            lapply(sessnames, function(s) {
                d <- x[sess == s, names(x) != "Session", drop = FALSE]
                if (nrow(d)) d else NULL
            })
        }
        else bysession(x, what)
    }

    det  <- if (is.null(detections)) rep(list(NULL), nsess)
            else splitframe(detections, "detections")
    tel  <- if (inherits(telemetry, "capthist")) {              # from read.telemetry()
        xy <- telemetryxy(telemetry)
        if (!ms(telemetry)) xy <- list(xy)
        if (length(xy) != nsess) stop ("telemetry capthist has a different number of sessions")
        if (!is.null(names(xy)) && all(names(xy) %in% sessnames)) xy <- xy[sessnames]
        unname(xy)
    } else splitframe(telemetry, "telemetry")
    mk   <- bysession(marked, "marked")
    cov  <- splitframe(covariates, "covariates")

    ## fixes as a named list of data frames (x, y) in sorted ID order
    fixlist <- function (z) {
        if (is.null(z)) return(structure(list(), names = character(0)))
        if (is.data.frame(z)) {
            if (!all(c("ID", "x", "y") %in% names(z))) stop ("telemetry needs columns ID, x, y")
            if (anyNA(z[, c("x", "y")])) stop ("missing values not allowed in telemetry")
            z <- split(z[, c("x", "y")], as.character(z$ID))
        }
        z <- z[order(names(z))]
        lapply(z, function (d) {
            d <- as.data.frame(d)[, 1:2]
            names(d) <- c("x", "y")
            rownames(d) <- NULL
            d
        })
    }
    anymarked <- any(!sapply(mk, is.null))
    covtemplate <- NULL
    for (d in cov) if (!is.null(d)) {
        d <- as.data.frame(d)
        d <- d[, !(names(d) %in% c("ID", "Session")), drop = FALSE]
        if (ncol(d) > 0) { covtemplate <- d[0, , drop = FALSE]; break }
    }

    ## ---- one session ----
    one <- function (tr, d, z, m, cv, sess) {
        xy  <- fixlist(z)
        ids <- names(xy)
        n   <- length(ids)
        K   <- nrow(tr)
        us  <- usage(tr)
        S   <- if (!is.null(us)) ncol(us) else if (!is.null(noccasions)) noccasions
               else stop ("need 'noccasions' because traps has no usage")
        if (is.null(us)) us <- matrix(1, K, S)
        if (nrow(us) != K) stop ("usage of traps has the wrong number of detectors")
        if (!is.null(d) && nrow(d) > 0) {
            if (!all(c("ID", "Occasion", "Detector") %in% names(d)))
                stop ("detections needs columns ID, Occasion, Detector")
            id <- as.character(d$ID)
            if (!all(id %in% ids))
                stop ("session ", sess, ": detected animal(s) without telemetry: ",
                      paste(head(setdiff(unique(id), ids), 5), collapse = ", "))
            k <- if (is.numeric(d$Detector)) as.integer(d$Detector)
                 else match(as.character(d$Detector), rownames(tr))
            s <- as.integer(d$Occasion)
            if (anyNA(k) || any(k < 1 | k > K))
                stop ("session ", sess, ": unknown Detector in detections")
            if (anyNA(s) || any(s < 1 | s > S))
                stop ("session ", sess, ": Occasion out of range in detections")
            if (any(us[cbind(k, s)] == 0))
                warning ("session ", sess, ": detections at detectors with zero usage")
        }
        else {
            id <- character(0); k <- s <- integer(0)
        }

        ## detection array: n x (S + 1) x (K + 1)
        arr <- array(0L, dim = c(n, S + 1, K + 1), dimnames = list(ids, 1:(S + 1), NULL))
        if (length(id)) {
            tab <- tabulate((match(id, ids) - 1L) + n * (s - 1L) + n * (S + 1L) * (k - 1L) + 1L,
                            nbins = length(arr))
            if (dettype == "proximity" && any(tab > 1)) {
                ## binary detection: one record per animal, occasion and detector
                warning ("session ", sess, ": ", sum(tab > 1), " animal-occasion-detector ",
                         "combination(s) with repeated detections set to 1 for proximity ",
                         "detectors")
                tab <- pmin(tab, 1L)
            }
            arr[] <- arr[] + tab
        }
        if (n) arr[, S + 1, K + 1] <- as.integer(vapply(xy, nrow, 1L))
        class(arr) <- "capthist"

        ## traps: cameras plus a copy of the first as the notional telemetry detector
        ## (as addTelemetry)
        oldtraps <- tr
        usage(oldtraps) <- NULL
        nulltraps <- subset(oldtraps, 1)
        traps(arr) <- suppressWarnings(rbind(oldtraps, nulltraps))
        detector(traps(arr)) <- c(rep(dettype, S), "telemetry")
        newusge <- matrix(0, nrow = K + 1, ncol = S + 1)
        newusge[1:K, 1:S] <- us
        newusge[K + 1, S + 1] <- 1
        markocc(traps(arr)) <- c(rep(0, S), 1)
        usage(traps(arr)) <- newusge
        telemetrytype(traps(arr)) <- "marking"
        session(arr) <- sess

        ## collar-active exposure
        if (anymarked) {
            if (is.null(m) || n == 0) m <- matrix(0, 0, S)
            else {
                m <- as.matrix(m)
                if (ncol(m) != S)
                    stop ("session ", sess, ": 'marked' should have one column per ",
                          "sighting occasion")
                if (is.null(rownames(m))) {
                    if (nrow(m) != n) stop ("session ", sess, ": 'marked' needs rownames = ID")
                    rownames(m) <- ids
                }
                if (!all(ids %in% rownames(m)))
                    stop ("session ", sess, ": 'marked' lacks collared animal(s): ",
                          paste(head(setdiff(ids, rownames(m)), 5), collapse = ", "))
                m <- m[ids, , drop = FALSE]
            }
            marked(arr) <- m
        }

        ## individual covariates
        if (!is.null(covtemplate)) {
            if (is.null(cv) || n == 0) cv <- covtemplate
            else {
                cv <- as.data.frame(cv)
                cid <- if ("ID" %in% names(cv)) as.character(cv$ID) else rownames(cv)
                if (!all(ids %in% cid))
                    stop ("session ", sess, ": covariates lack collared animal(s): ",
                          paste(head(setdiff(ids, cid), 5), collapse = ", "))
                cv <- cv[match(ids, cid), !(names(cv) %in% c("ID", "Session")), drop = FALSE]
                rownames(cv) <- ids
            }
            covariates(arr) <- cv
        }

        ## telemetry fixes
        attr(arr, "telemetryxy") <- xy
        arr
    }

    ## ---- sighting counts Tu, Tm: one matrix (K x S) per session from a matrix, a data frame
    ## or file with a leading session column (all sessions stacked, K rows per session as in
    ## addSightings()), or a list with one element per session ----
    Ks <- vapply(trapslist, nrow, 1L)
    Ss <- vapply(trapslist, function (tr) {
        us <- usage(tr)
        if (!is.null(us)) ncol(us) else if (!is.null(noccasions)) noccasions else NA_real_
    }, 1)
    sight <- function (x, what) {
        if (is.null(x)) return(rep(list(NULL), nsess))
        asmat <- function (m, i) {                       # one session's matrix, dimensions checked
            m <- as.matrix(m)
            storage.mode(m) <- "numeric"
            dimnames(m) <- NULL
            if (nrow(m) != Ks[i] || (!is.na(Ss[i]) && ncol(m) != Ss[i]))
                stop (what, " for session ", sessnames[i], " has ", nrow(m), " rows and ",
                      ncol(m), " columns; expected ", Ks[i], " detectors and ",
                      if (is.na(Ss[i])) "a constant number of" else Ss[i], " occasions")
            m
        }
        readfile <- function (f) read.table(f, ...)
        ## one file, leading session column
        if (is.character(x) && length(x) == 1) x <- readfile(x)
        if (is.data.frame(x)) {
            lab  <- as.character(x[, 1])
            ulab <- unique(lab)                          # in order of first appearance
            if (length(ulab) != nsess)
                stop (what, ": ", length(ulab), " session label(s) in the first column (",
                      paste(ulab, collapse = ", "), ") but ", nsess, " session(s)")
            if (all(ulab %in% sessnames))
                ord <- match(sessnames, ulab)            # labels are the session names
            else {
                if (nsess > 1)
                    message (what, ": session labels (", paste(ulab, collapse = ", "),
                             ") differ from the session names (",
                             paste(sessnames, collapse = ", "), "); matched by position")
                ord <- seq_len(nsess)
            }
            return (lapply(seq_len(nsess), function (i)
                asmat(x[lab == ulab[ord[i]], -1, drop = FALSE], i)))
        }
        if (is.matrix(x)) {
            if (nsess != 1) stop (what, " should be a list with one element per session")
            return (list(asmat(x, 1)))
        }
        if (is.list(x)) {
            if (length(x) != nsess)
                stop (what, " should be a list with one element per session")
            if (!is.null(names(x)) && all(names(x) %in% sessnames)) x <- x[sessnames]
            return (lapply(seq_len(nsess), function (i) {
                xi <- x[[i]]
                if (is.null(xi)) NULL
                else if (is.character(xi)) asmat(readfile(xi)[, -1, drop = FALSE], i) # file
                else asmat(xi, i)
            }))
        }
        stop (what, " should be a matrix, a data frame or file with a leading session ",
              "column, or a list")
    }

    out <- vector("list", nsess)
    for (i in seq_len(nsess))
        out[[i]] <- one(trapslist[[i]], det[[i]], tel[[i]], mk[[i]], cov[[i]], sessnames[i])

    ## ---- combine and add sightings ----
    CH <- if (MS) {
        names(out) <- sessnames    # MS.capthist() takes session names from names of arguments
        x <- do.call(MS.capthist, out)
        if (!is.null(covtemplate)) x <- shareFactorLevels(x)
        x
    } else out[[1]]
    if (!is.null(Tu) || !is.null(Tm)) {
        tu <- sight(Tu, "Tu")
        tm <- sight(Tm, "Tm")
        CH <- if (MS) addSightings(CH, unmarked = if (is.null(Tu)) NULL else tu,
                                   nonID = if (is.null(Tm)) NULL else tm, verify = FALSE)
              else    addSightings(CH, unmarked = tu[[1]], nonID = tm[[1]], verify = FALSE)
    }
    ## the count model applied to what look like presence/absence data: every sighting count
    ## of unmarked animals is 0 or 1
    if (dettype == "count") {
        tu <- if (MS) lapply(CH, secr::Tu) else list(secr::Tu(CH))
        tu <- unlist(lapply(tu, function (x) if (is.null(x)) NULL else as.vector(x)))
        if (length(tu) > 0 && sum(tu) > 0 && max(tu) <= 1)
            warning ("all counts of unmarked sightings (Tu) are 0 or 1 and detectors are ",
                     "'count': if these are presence/absence (binary) data, the count model ",
                     "is not appropriate (biased when detection is high); build with ",
                     "'proximity' detectors for the pooled-binary model", call. = FALSE)
    }
    if (verify) verify(CH)
    CH
}
