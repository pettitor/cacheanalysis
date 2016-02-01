C <- rep(4,5)
l <- c(10,5,2,1,1,1,1,1,1,1)
l <- l/sum(l)
rho <- rep(800,5)
d <- rep(3*60,10)
b <- rep(250,10)

n2 <- length(rho)

a <- hwca(l,b,d,rho,C)

X <- hwc(l,b,d,rho,C)

x <- rep(0,dim(X)[1])
for (m in 1:dim(X)[1]) x[m] <- any(X[m,]==1)

phit <- x*l*(1-pmax(0,a-1))

