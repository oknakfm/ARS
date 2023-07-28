t <- seq(0, 2*pi, (2*pi)/16)
tfine <- seq(0, 2*pi, (2*pi)/640)
X = cbind(cos(t), sin(t))
Xfine = cbind(cos(tfine), sin(tfine))

pdf(file="fig2.pdf", width=10, height=5.5)
par(mfrow=c(1,2))
plot(Xfine, type="l", xlab="x1(t)", ylab="x2(t)", xlim=c(-1,1), ylim=c(-1,1), col="grey",
     main="entire dynamics (x1(t),x2(t))")
par(new=T)
plot(X, type="p", xlab=" ", ylab=" ", xlim=c(-1,1), ylim=c(-1,1), xaxt="n", yaxt="n")
for(k in 1:(nrow(X)-1)){
  arrows(x0 = X[k,1], y0 = X[k,2], x1 = X[k+1,1], y1 = X[k+1,2], length = 0.15, angle = 20,
         code = 1, col = par("fg"), lty = par("lty"),
         lwd = par("lwd"))
}
plot(tfine, Xfine[,1], type="l", xlab="t", ylab="z(t)=x1(t)", xlim=c(0,2*pi), ylim=c(-1,1), col="grey", 
     main="partial observation z(t)=x1(t)")
par(new=T)
plot(t, X[,1], type="p", xlab=" ", ylab=" ", xlim=c(0,2*pi), ylim=c(-1,1), xaxt="n", yaxt="n")
for(k in 1:(nrow(X)-1)){
  arrows(x0 = t[k], y0 = X[k,1], x1 = t[k+1], y1 = X[k+1,1], length = 0.15, angle = 20,
         code = 2, col = par("fg"), lty = par("lty"),
         lwd = par("lwd"))
}
dev.off()