################################################################################
## package 'secr'
## bysession.R
## 2026-11-26 experimental: evaluate the sessions of a multi-session likelihood in parallel
##            (details$bysession = TRUE)
################################################################################

## One cluster is started by secr.fit (secr_bysession_start) before the likelihood is first
## evaluated and stopped when secr.fit exits (secr_bysession_stop). It is held in
## .localstuff$bs. The session likelihood functions secr_sessionLL (general) and
## secr_sessionLLfast take the session data and all other inputs as arguments. Each worker is
## assigned a fixed set of sessions (round-robin) and receives once the objects that do not
## change during the fit (the session data, designs and other fixed inputs). At each
## evaluation only the objects that depend on beta are sent (secr_bysession_LL), and the
## worker returns the log likelihood of its sessions and its rows of .localstuff$Eng.
##
## details$ncores is the total number of cores. With nsession sessions and nclust =
## min(nsession, ncores) workers, each worker uses max(1, ncores %/% nsession) threads for
## the C++ code (as list.secr.fit does for models).

## arguments of the session functions that do not change during a fit
.bsstaticnames <- c("parindx", "CL", "detectfn", "learnedresponse", "design", "design0", "fixed")

## details$bysession is FALSE (serial), TRUE (one worker per session, up to ncores) or a
## number of workers (sessions are shared among them, balanced by size)
secr_bysession_on <- function (x) {
    isTRUE(x) || (is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 1)
}

## Relative cost of a session likelihood (animals x detectors x mask points; a rough guide
## used only to balance sessions among workers)
secr_sessioncost <- function (d) {
    as.numeric(d$m) * as.numeric(d$K) * (1 + as.numeric(d$nc))
}

## Assign sessions to nclust workers, largest first to the least loaded worker
secr_assignsessions <- function (cost, nclust) {
    load <- numeric(nclust)
    idx <- vector("list", nclust)
    for (i in order(cost, decreasing = TRUE)) {
        w <- which.min(load)
        idx[[w]] <- c(idx[[w]], i)
        load[w] <- load[w] + cost[i]
    }
    lapply(idx, sort)
}

secr_bysession_start <- function (data, ncores, bysession = TRUE) {
    secr_bysession_stop()
    nsession <- length(data)
    nclust <- if (is.numeric(bysession)) min(round(bysession), nsession) else nsession
    if (nclust > ncores) {
        if (is.numeric(bysession)) warning("details$bysession exceeds ncores; using ncores workers")
        nclust <- ncores
    }
    if (nclust < 2) {
        message("details$bysession ignored: needs more than one session and more than one core")
        return(invisible(NULL))
    }
    if (ncores > parallel::detectCores())
        warning("ncores exceeds the number of available cores")
    threads <- max(1, ncores %/% nclust)
    clustertype <- getOption("secr.clustertype",
        default = if (.Platform$OS.type == "unix" && Sys.getenv("RSTUDIO") != "1")
            "FORK" else "PSOCK")
    clust <- parallel::makeCluster(nclust, type = clustertype)
    parallel::clusterEvalQ(clust, library(secr))
    parallel::clusterCall(clust, secr_bs_threads, threads)
    idx <- secr_assignsessions(vapply(data, secr_sessioncost, numeric(1)), nclust)
    message("secr.fit is using a parallel cluster with ", nclust, " workers for ", nsession,
            " sessions, ", threads, " threads per worker")
    .localstuff$bs <- list(clust = clust, idx = idx, nsession = nsession, threads = threads,
                           static = NULL, data = NULL)
    invisible(NULL)
}

secr_bysession_stop <- function () {
    bs <- .localstuff$bs
    if (!is.null(bs)) {
        try(parallel::stopCluster(bs$clust), silent = TRUE)
        .localstuff$bs <- NULL
    }
    invisible(NULL)
}

## functions run on the workers
secr_bs_threads <- function (n) {
    secr::setNumThreads(n)
    NULL
}

secr_bs_set <- function (x) {
    .localstuff$bsstatic <- x
    NULL
}

secr_bs_worker <- function (fn, dyn, ngroup, iter) {
    st <- .localstuff$bsstatic
    f <- get(fn, envir = asNamespace("secr"))
    .localstuff$Eng <- matrix(0, nrow = st$nsession, ncol = ngroup)
    .localstuff$iter <- iter
    ll <- vapply(st$data, function(d) sum(do.call(f, c(list(d), st$static, dyn))), numeric(1))
    list(ll = ll,
         sessnum = vapply(st$data, function(d) as.integer(d$sessnum), integer(1)),
         Eng = .localstuff$Eng,
         post = .localstuff$markedpost)
}

## Called from the likelihood functions in place of the serial loop over sessions.
## fn is the name of the session function; sessargs the list of its arguments other than 'data'.
secr_bysession_LL <- function (fn, data, sessargs, ngroup) {
    bs <- .localstuff$bs
    isstatic <- names(sessargs) %in% .bsstaticnames
    static <- sessargs[isstatic]
    dyn <- sessargs[!isstatic]
    if (is.null(bs$static) || !identical(bs$static, static) || !identical(bs$data, data)) {
        ## (re)send the unchanging objects; identical() is quick when the objects are the same
        ## (pointer comparison) as in repeated evaluations within a fit
        parallel::clusterApply(bs$clust,
            lapply(bs$idx, function(i) list(static = static, data = data[i], nsession = length(data))),
            secr_bs_set)
        .localstuff$bs$static <- static
        .localstuff$bs$data <- data
    }
    dyn$details$ncores <- bs$threads
    dyn$details$bysession <- FALSE
    res <- parallel::clusterCall(bs$clust, secr_bs_worker, fn, dyn, ngroup, .localstuff$iter)
    nsession <- length(data)
    ll <- numeric(nsession)
    Eng <- matrix(0, nrow = nsession, ncol = ngroup)
    for (r in res) {
        ll[r$sessnum] <- r$ll
        Eng[r$sessnum, ] <- r$Eng[r$sessnum, ]
        if (nsession %in% r$sessnum) .localstuff$markedpost <- r$post   # as the serial loop (last session)
    }
    .localstuff$Eng <- Eng
    sum(ll)
}
