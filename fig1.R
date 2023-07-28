gen_data <- function(n=10, span=0.1, sd=0){
  angles = span * (1:n)
  data = cbind(cos(angles), sin(angles),angles) + cbind(matrix(rnorm(n=n*2, sd=sd),n,2),0)
  colnames(data) = c("X1","X2","angles")
  return(data.frame(data))
}

rotation_matrix <- function(theta=0){
  return(matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),2,2))
}

AR <- function(data, def=FALSE, res=FALSE){
  if(def){
    n = length(data)
    x = data[1:(n-1)]
    y = data[2:n]
    a = sum(x*y)/sum(x*x)
    if(res){
      return(sum((y-x*a)^2))
    }else{
      return(a)
    }
  }else{
    n = nrow(data)
    data = subset(data, select=c("X1", "X2"))
    X = as.matrix(data[1:(n-1),])
    Y = as.matrix(data[2:n,])
    A = solve(t(X) %*% X) %*% (t(X) %*% Y)
    if(res){
      return(sum((Y-X%*%A)^2))
    }else{
      return(A)
    }
  }
}

pred <- function(x0, A, m=10){
  tmp = NULL; x0 = matrix(x0,1,2)
  for(itr in 1:m){
    x0 = x0 %*% A
    tmp = rbind(tmp, x0)
  }
  return(data.frame(tmp))
}

## ============
## dataset
## ============
len = 9
span=0.3; n=len/span; sd=0
data = gen_data(n=n, span=span, sd=sd)  

## ============
## estimation
## ============
At = rotation_matrix(theta=span) # true A
A = AR(data=data) # with fully-observed time-series
a = AR(data=data$X1, def=TRUE) # with partially observed
f <- function(pm) AR(data=data.frame(X1=data$X1, X2=pm), res=TRUE)
.opt = optim(par=rnorm(n), fn=f, method="BFGS")
if(!.opt$convergence) cat("converged") else ("max_itr has been reached")
Ad = AR(data=data.frame(X1=data$X1, X2=.opt$par))

## ===========
## prediction
## ===========
m = 15
angles = span * (1:n)
angles_test = span * (n+(1:m))
x0 = as.matrix(subset(data, select=c("X1", "X2")))[n,]
pred.At = pred(x0=x0, A=At, m=m)
pred.A = pred(x0=x0, A=A, m=m)
x0 = c(data$X1[n], .opt$par[n])
pred.Ad = pred(x0=x0, A=Ad, m=m)

## ===========
## plot
## ===========
xl = range(angles, angles_test)
yl = c(-1.5,1.5)
legends = c("AR (conventional)", "ARS (proposal)")
colors = c("red", "blue")

pdf(file="fig1.pdf", width=7, height=5)
par(mfrow=c(1,1))
plot(angles, data$X1, xlim=xl, ylim=yl, 
     xlab="t", ylab=" ", type="b")
par(new=T)
plot(angles_test, data$X1[n] * (a^(1:m)), xlim=xl, ylim=yl, 
     xlab=" ", ylab=" ", xaxt="n", yaxt="n", type="b", col="red", pch="+")
par(new=T)
plot(angles_test, pred.Ad[,1], xlim=xl, ylim=yl, 
     xlab=" ", ylab=" ", xaxt="n", yaxt="n", type="b", col="blue", pch="+")
abline(v=angles[n])
text(x=4*span, y=1.4, labels="[In-sample]")
text(x=angles[n]+6*span, y=1.4, labels="[Out-sample]")
legend("bottomleft", legend=legends, col=colors, lty=1, pch="+", bg="white")
dev.off()