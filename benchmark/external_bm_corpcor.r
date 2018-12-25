require(corpcor)
require(matrixcalc)

C_20 = as.matrix(read.table("tmpmat/C_20.csv"))
C_40 = as.matrix(read.table("tmpmat/C_40.csv"))
C_200 = as.matrix(read.table("tmpmat/C_200.csv"))
C_400 = as.matrix(read.table("tmpmat/C_400.csv"))

X_40x20   = as.matrix(read.table('tmpmat/X_40x20.csv'))
X_20x40   = as.matrix(read.table('tmpmat/X_20x40.csv'))
X_400x200 = as.matrix(read.table('tmpmat/X_400x200.csv'))
X_200x400 = as.matrix(read.table('tmpmat/X_200x400.csv'))

C = list(C_20, C_40, C_200, C_400)
X = list(X_40x20, X_20x40, X_400x200, X_200x400)

times = numeric(4)
res = numeric(4)

for (i in c(1, 2, 3, 4)){
  Xi = X[[i]]
  Ci = C[[i]]
  start = Sys.time()
  C_est = cov.shrink(Xi, lambda.var=0.0)
  times[i] = as.numeric(Sys.time()-start)
  res[i] = frobenius.norm(C_est - Ci)
}

