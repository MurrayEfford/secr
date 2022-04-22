superN <- function (N1, phi, lambda) {
    J <- length(lambda)
    Nj <- cumprod(c(N1,lambda))[1:J]
    B <- Nj - c(0,phi * Nj)[1:J]
    list(B=B, beta =B/sum(B), Nj=Nj, superN=sum(B))
}

beta.superN <- function (phi, lambda) {
    J <- length(lambda)
    Nj <- cumprod(c(1,lambda))[1:J]
    B <- Nj - c(0,phi * Nj)[1:J]
    B/sum(B)
}
beta.superN(rep(0.7,5), rep(1.1,5))


rmultinomB <- function (N1, phi, lambda) {
    J <- length(lambda)
    Nj <- cumprod(c(N1,lambda))[1:J]
    B <- Nj - c(0,phi * Nj)[1:J]
    superN <- sum(B)
    beta <- B/sum(B)
    # list(B=B, beta =beta, Nj=Nj, superN=superN)
    rmultinom(1, superN, beta)
}
#rmultinomB(100, rep(0.7,5), rep(1.1,5))

proj <- function(N1, phi, lambda) {
    J <- length(lambda)
    f <- lambda-phi
    N <- numeric(J)
    N[1] <- rpois(1,N1)
    for (j in 2:J) {
        surv <- rbinom(1,N[j-1],phi[j-1])
        B <- rpois(1, f[j-1] * N[j-1])
        N[j] <- surv + B
    }
    N
}
# out <- matrix(NA, 100000, 5)
# for (i in 1:100000) out[i,] <- proj(100, rep(0.7,5), rep(1.1,5))
# apply(out,2,mean)
# 1.1^4
# apply(out,2,var)
