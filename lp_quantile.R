
#test
library(lpSolve)

n = dim(X)[1]; d = dim(X)[2] - 1
X_star = matrix(NA, nrow=n, ncol=n*(d+1))
for(i in 0:d){
  X_star[, (1+i*n):((i+1)*n)] <- diag(X[, i+1])
}

A1 = cbind(diag(n),0*diag(n), -1*X_star,matrix(0,n,n^2*(d+1)))
b1= y
A2 = cbind(diag(n),0*diag(n), X_star,matrix(0,n,n^2*(d+1)))
b2 = -1*y
b3 = b4 = rep(0, n^2*(d+1))

A3s = matrix(0, n^2, n)
for(i in 1:n){
  A3s[,i] = rep(1,n)
  A3s[(((i-1)*n)+1):(i*n),] = 
    A3s[(((i-1)*n)+1):(i*n),] + -1*diag(n)
}

A3beta = matrix(0,n^2*(d+1),n*(d+1))
for(i in 1:(d+1)){
  A3beta[((i-1)*n^2+1):(i*n^2),((i-1)*n+1):(i*n)] = A3s
}
A3 = cbind(matrix(0,n^2*(d+1),n), matrix(0,n^2*(d+1),n), A3beta, diag(n^2*(d+1)))

A4 = cbind(matrix(0,n^2*(d+1),n), matrix(0,n^2*(d+1),n), -1*A3beta, diag(n^2*(d+1)))

Amat = rbind(A1, A2, A3, A4);dim(Amat)
bvec = c(b1,b2,b3,b4)


tau = 0.5
lam=0

f_obj = c(rep(tau, n), rep(1-tau, n), rep(0,n*(d+1)), rep(lam,n^2*(d+1)))
f_con = Amat
f_dir = ">="
f_rhs = bvec
lpsln = lp("min", f_obj, f_con, f_dir, f_rhs)
betaest = lpsln$solution[(2*n+1):((2*n+n*(d+1)))]

plot(beta_true, pch=16, col=2, ylim=c(-5,5))
#points(beta_local, pch="l", col=3)
points(betaest, pch="p", col=4)
abline(v=seq(ntr, d*ntr, ntr), col=3, lty=3)

