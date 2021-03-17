#Unmarked Data Examples

#Chapter 7, Applied Hierarchical Models
library(unmarked)
data(ovendata)

alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
head(alfl)

alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",package="unmarked"), row.names=1)
head(alfl.covs)
alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3, sep="")
alfl$captureHistory <- factor(alfl$captureHistory,levels=c("001", "010", "011", "100", "101", "110", "111"))

alfl.v1 <- alfl[alfl$survey==1,]
alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
head(alfl.H1, 5)

crPiFun <- function(p) {
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
  cbind("001" = (1-p1) * (1-p2) * p3,
        "010" = (1-p1) * p2     * (1-p3),
        "011" = (1-p1) * p2     * p3,
        "100" = p1     * (1-p2) * (1-p3),
        "101" = p1     * (1-p2) * p3,
        "110" = p1     * p2     * (1-p3),
        "111" = p1     * p2     * p3)
}

p <- matrix(0.2,5,3)
crPiFun(p)

rowSums(crPiFun(p))

intervalMat <- matrix(c('1','2','3'), 50, 3, byrow=TRUE)
class(alfl.H1) <- "matrix"

umf.cr1 <- unmarkedFrameMPois(y=alfl.H1,
            siteCovs=alfl.covs[,c("woody", "struct", "time.1", "date.1")],
            obsCovs=as.data.frame(list(interval=intervalMat)), obsToY=o2y, piFun="crPiFun")


##### Example gmultmix in unmarked manual ####

# Simulate data using the multinomial-Poisson model with a
# repeated constant-interval removal design.

n <- 100  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods
lam <- 3
phi <- 0.5
p <- 0.3
#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals available for detection

for(i in 1:n) {
  N[i,] <- rbinom(T, M[i], phi)
  y[i,,1] <- rbinom(T, N[i,], p)    # Observe some
  Nleft1 <- N[i,] - y[i,,1]         # Remove them
  y[i,,2] <- rbinom(T, Nleft1, p)   # ...
  Nleft2 <- Nleft1 - y[i,,2]
  y[i,,3] <- rbinom(T, Nleft2, p)
  y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
}
  umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")
  (m1 <- gmultmix(~1, ~1, ~1, data=umf1, K=30))
  
  backTransform(m1, type="lambda")        # Individuals per plot
  backTransform(m1, type="phi")           # Probability of being avilable
  (p <- backTransform(m1, type="det"))    # Probability of detection
  p <- coef(p)
  
  # Multinomial cell probabilities under removal design
  c(p, (1-p) * p, (1-p)^2 * p)
  
  #or more generally
  head(getP(m1))
  
  # Empirical Bayes estimates of super-population size
  re <- ranef(m1)
  plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)