require(plotly)
require("parallel"); require("foreach"); require("doParallel"); require("doRNG")

## ==========
## settings
## ==========
set.seed(123)
# cmode = c("computation", "plot")
cmode = "plot"
N.Cores = detectCores() ## CPU cores for parallel computation

## ==========
## folders
## ==========
DIR = getwd()
OUT = paste0(DIR,"/Output")
dir.create(path=OUT, showWarnings=FALSE)

## ==========
## functions
## ==========
expand <- function(x){
  tmp <- x %*% t(x)
  return(append(x, tmp[upper.tri(tmp)]))
}

gen_det_seq <- function(n=10, d=3, r=2, E=NULL, x0=NULL, interaction=FALSE, 
                        n_test=10, step=10, burn_in=100){
  if(is.null(x0)) x0 = rnorm(d)
  if(is.null(E)){
    if(interaction){
      E = matrix(rnorm(d * d*(d+1)/2), d, d*(d+1)/2)
    }else{
      E = matrix(rnorm(d*d), d, d)
    }
  }
  
  X = NULL; tmp_x = x0
  for(i in 1:(burn_in + n * step)){
    if(interaction){
      tmp_x = E %*% expand(tmp_x)
    }else{
      tmp_x = E %*% tmp_x
    }
    if((i > burn_in) && ((i-burn_in) %% step == 0)) X = rbind(X, as.vector(tmp_x))
  }  
  
  Y = NULL
  for(i in 1:(n_test * step)){
    if(interaction){
      tmp_x = E %*% expand(tmp_x)
    }else{
      tmp_x = E %*% tmp_x
    }
    if(i %% step == 0) Y = rbind(Y, as.vector(tmp_x))
  }
  
  obs = if(r==1) matrix(X[,1:r],n,1) else X[,1:r]
  mis = if(d==(r+1)) matrix(c(X[,(r+1):d]),n,1) else X[,(r+1):d]
  
  return(list(
    X=X,
    obs=obs,
    mis=mis, 
    Y=Y, 
    E=E
  ))
}

ARS <- function(obs=NULL, mis=NULL, s=1, interaction=FALSE, eta=0.1, step=1, 
                maxit=10**4){
  
  n = nrow(obs)
  r = ncol(obs)
  
  if(s>0){
    param = if(is.null(mis)) rnorm(n*s, sd=sd(obs)) else as.vector(mis)

    res <- function(param){
      mis = matrix(param, n, s)
      X = cbind(obs, mis)
      if(interaction){
        X = t(apply(X, 1, expand))
      } 
      Xm = X[1:(n-step),]
      Xp = X[(1+step):n,1:(r+s)]
      E = solve(t(Xm) %*% Xm) %*% t(Xm) %*% Xp
      Yp = Xm %*% E  
      
      mse = mean((Yp - Xp)^2)
      pen = mean((mis[2:n,]-mis[1:(n-1),])^2)
      
      return(mse + eta * pen)
    }
    if(is.null(maxit)){
      X = cbind(obs, matrix(param, n, s))
    }else{
      .opt <- optim(par=param, fn=res, method="BFGS", control=list(maxit=maxit))
      X = cbind(obs, matrix(.opt$par, n, s))
    }
  }else{
    X = obs
  }
  
  if(interaction && r+s>1){
    Z = t(apply(X, 1, expand))
  }else{
    Z = X
  }
  
  Xm = Z[1:(n-step),]
  Xp = Z[(1+step):n,1:(r+s)]
  E = t(solve(t(Xm) %*% Xm) %*% t(Xm) %*% Xp)
  
  ret = if(s>0 && !is.null(maxit)){
    list(X=X, E=E, interaction=interaction, convergence=.opt$convergence)
  }else{
    list(X=X, E=E, interaction=interaction)
  }
  return(ret)
}

LorenzE <- function(alpha=10, beta=28, gamma=8/3, Delta=0.1){
  tmp1 <- cbind(diag(3), matrix(0,3,3))
  tmp2 <- rbind(c(-alpha,alpha,0,0,0,0),
        c(beta,-1, 0,0,-1,0),
        c(0,0,-gamma,1,0,0))
  return(tmp1 + Delta * tmp2)
}

prediction <- function(.ARS = NULL, n=10, step=1){
  X_all = .ARS$X
  
  pred_X = NULL; tmp_x = X_all[nrow(X_all)+1-step,]
  for(i in 1:n){
    if(.ARS$interaction){
      tmp_x = .ARS$E %*% expand(tmp_x)
    }else{
      tmp_x = .ARS$E %*% tmp_x
    }
    tmp_x = as.vector(tmp_x)
    pred_X = rbind(pred_X, tmp_x)
    X_all = rbind(X_all, tmp_x)
    tmp_x = X_all[nrow(X_all)+1-step,]
  }  
  
  return(pred_X)
}


## ==========================
## experiments
## ==========================
n_exp = 10 ## number of experiments
for(setting in 1:2){
  ## setting 1 = circular motion
  ## setting 2 = Lorenz dynamics (linear approx.)
  
  for(sd_id in 1:2){
    set.seed(123) ## specify a random seed
    if(setting==1){
      sd = c(0, 0.01)
      n = 100; n_test=n_pred=30; step=1
      d=2; r=1; x0=c(0,1); deg=0.05; 
      E = matrix(c(cos(deg), sin(deg), -sin(deg), cos(deg)),2,2)
      interaction=FALSE
      est_int = FALSE
    }else if(setting==2){
      sd = c(0, 0.01)
      n = 100; n_test=n_pred=100; step=1
      Delta = 0.01; d=3; r=2; x0=c(0.25,0.25,0.25)
      E = LorenzE(Delta=Delta)
      interaction=TRUE
      est_int = TRUE  
    }
    
    det_seq = gen_det_seq(n=n, d=d, r=r, E=E, x0=x0, interaction=interaction, 
                          n_test=n_test, step=step, burn_in=100)
    
    if("computation" %in% cmode){
      registerDoParallel(N.Cores)
      registerDoRNG(123)
      foreach(exp_id = 1:n_exp) %dopar% {
        obs = det_seq$obs + matrix(rnorm(n*r,sd=sd[sd_id]),n,r)
        
        mis = det_seq$mis
        mis = mis + rnorm(n=length(mis), sd=1)
        
        .ARS0F = ARS(obs=obs, s=0, interaction=FALSE, eta=eta)
        .ARS1F1 = ARS(obs=obs, mis=mis, s=1, interaction=est_int, eta=0, step=1)
        
        pred_0F = prediction(.ARS=.ARS0F, n=n_pred)
        pred_1F1 = prediction(.ARS=.ARS1F1, n=n_pred, step=1)
        
        E0F = apply(sapply(1:r, function(k) abs(pred_0F[,k] - det_seq$Y[,k])), 1, mean)
        E1F1 = apply(sapply(1:r, function(k) abs(pred_1F1[,k] - det_seq$Y[,k])), 1, mean)
        
        filename = paste0(OUT,"/setting",setting,"_sd",sd_id,"_exp",exp_id,".RData")
        
        save(file=filename, 
             .ARS0F, .ARS1F1, 
             E0F, E1F1, 
             pred_0F, pred_1F1, 
             obs, mis)
        
      }
      stopImplicitCluster()
      
    }
    
    if("plot" %in% cmode){
      L0F = L1F1 = L1F2 = L1F3 = NULL
      
      for(exp_id in 1:n_exp){
        filename = paste0(OUT,"/setting",setting,"_sd",sd_id,"_exp",exp_id,".RData")
        load(file=filename)
        L0F = rbind(L0F, E0F)
        L1F1 = rbind(L1F1, E1F1)
      }
      
      tm <- function(x) mean(x)
      ts <- function(x) sd(x)
      
      r1F1 = apply(L1F1/L0F, 2, tm); s1F1 = apply(L1F1/L0F, 2, ts)
      
      ind = (25/5) * 1:5
      
      ct = NULL
      ct = append(ct, paste0("n=",n,"/sd=",sd[sd_id],"\n"))
      ct = append(ct,paste0("$", signif(r1F1[ind],digits=3), "\\pm", signif(s1F1[ind],digits=3), "$ &"))
      filename = paste0(OUT,"/setting",setting,"_sd",sd_id,".txt")
      write(x=ct, file=filename)
      
      names = c("Observations (for train.)", "True (for test)", "AR", "ARS (Initial)", "ARS (Trained)")
      colors = c("grey", "black", "red", "yellow", "blue")
      pchs = c(4, 4, 1, 1, 1)
      
      filename = paste0(OUT,"/setting",setting,"_sd",sd_id,"_exp1.RData")
      load(file=filename)
      
      if(r==1){
        .ARSI = ARS(obs=obs, mis=mis, s=1, interaction=est_int, eta=0, step=1, maxit=NULL)
        pred_I = prediction(.ARS=.ARSI, n=n_pred, step=1)
        
        pdf(file=paste0(OUT,"/setting",setting,"_sd",sd_id,".pdf"), width=7, height=6)
        par(mar = c(2, 2, 2, 1))
        par(oma = c(2, 2, 2, 1))
        yl = range(obs, pred_0F[,1],pred_1F1[,1], pred_I[,1])
        xmax = (n+n_pred)
        xl = c(-18,xmax)
        plot(1:n, obs, xlim=xl, ylim=yl, xlab="k", ylab="z1", type="b", 
             col=colors[1], pch=pchs[1], main=paste0("AR/ARS computed in Instance 1 (n=",n,")"))
        par(new=T)  
        plot(n + (1:n_pred), det_seq$Y[1:n_pred,1], xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[2], pch=pchs[2],
             xaxt="n",yaxt="n")
        par(new=T)
        plot(n + (1:n_pred), pred_0F, xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[3], pch=pchs[3],
             xaxt="n",yaxt="n")
        
        par(new=T)
        plot(n + (1:n_pred), pred_I[,1], xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[4], pch=pchs[4],
             xaxt="n",yaxt="n")
        par(new=T)
        plot(n + (1:n_pred), pred_1F1[,1], xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[5], pch=pchs[5],
             xaxt="n",yaxt="n")
        
        abline(v=n)
        
        legend("bottomleft", legend=names,col=colors,pch=pchs,bg="white")
        
        dev.off()
      }else if(r==2){
        .ARSI = ARS(obs=obs, mis=mis, s=1, interaction=est_int, eta=0, step=1, maxit=NULL)
        pred_I = prediction(.ARS=.ARSI, n=n_pred, step=1)
        
        pdf(file=paste0(OUT,"/setting",setting,"_sd",sd_id,".pdf"), width=6, height=6)
        
        xl = range(obs[,1], det_seq$Y[,1], pred_0F[,1],pred_1F1[,1], pred_I[,1])
        yl = range(obs[,2], det_seq$Y[,2], pred_0F[,2],pred_1F1[,2], pred_I[,2])
        
        plot(obs[,1], obs[,2], 
             xlim=xl, ylim=yl, xlab="z1", ylab="z2", type="b", 
             col=colors[1], pch=pchs[1], main=paste0("AR/ARS computed in Instance 1 (n=",n,")"))
        
        par(new=T)
        plot(det_seq$Y[1:n_pred,1], det_seq$Y[1:n_pred,2], xlim=xl, ylim=yl, xlab=" ", ylab=" ", 
             type="b", col=colors[2], pch=pchs[2],
             xaxt="n",yaxt="n")
        
        par(new=T)
        plot(pred_0F[,1], pred_0F[,2], xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[3], pch=pchs[3],
             xaxt="n",yaxt="n")
        par(new=T)
        plot(pred_I[,1], pred_I[,2], xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[4], pch=pchs[4],
             xaxt="n",yaxt="n")
        par(new=T)
        plot(pred_1F1[,1], pred_1F1[,2], xlim=xl, ylim=yl, xlab=" ", ylab=" ", type="b", col=colors[5], pch=pchs[5],
             xaxt="n",yaxt="n")
        
        legend("topleft", legend=names,col=colors,pch=pchs,bg="white")
        dev.off()
      }
    }
  }
}