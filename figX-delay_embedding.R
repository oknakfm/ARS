n = 16
t <- seq(0, 2*pi, 2*pi/n)
t_delayed <- append(t[2:n], t[1])

X = cbind(cos(t), sin(t))
Y = cbind(cos(t), cos(t_delayed))


pdf(file="delay_embedding.pdf", width=9, height=5)
par(mfrow=c(1,2))
plot(X, type="p", xlab="x1(t)", ylab="x2(t)", main="circular motion")

for(k in 1:(n-1)){
  arrows(x0=X[k,1], y0=X[k,2], x1=X[k+1,1], y1=X[k+1,2], length = 0.15, angle = 20,
         code = 1, col = par("fg"), lty = par("lty"))
}

plot(Y, type="p", xlab="x1(t-1)", ylab="x1(t)", main="delay embedding")
for(k in 1:(n-1)){
  arrows(x0=Y[k,2], y0=Y[k,1], x1=Y[k+1,2], y1=Y[k+1,1], length = 0.15, angle = 20,
         code = 1, col = par("fg"), lty = par("lty"))
}

dev.off()


