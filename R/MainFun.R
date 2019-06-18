tr<-function(X) sum(diag(X));
negloglikeF<-function(parameter,K,Y,X)
{
  n=length(Y); p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  tmp1=(n-p)*log(sigma^2)
  tmp2=sum(log(abs(eigen(V0)$value)))
  tmp3=sum(log(abs(eigen(XTV0invX)$value)))
  tmp4= t(Y) %*% P %*% Y/sigma^2;
  like=tmp1+tmp2+tmp3+tmp4;
  like=-1/2*like;
  return(-like);
}
negfirstderisigmaF<-function(parameter,X,Y,K)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  firstderi=NULL;
  for(i in 1:length(K))
  {
    tmp=-1/2*(tr(P %*% K[[i]] *2 *dr[i])-1/sigma^2 * t(Y) %*% P %*% (K[[i]]*2*dr[i]) %*% P %*% Y)
    firstderi=c(firstderi,tmp)
  }
  firstderi=c(firstderi,(-(n-p)/2/sigma^2+t(Y) %*% P %*% Y/2/sigma^4)*(2*sigma));
  -firstderi;
}

loglikeF<-function(parameter,K,Y,X)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  tmp1=(n-p)*log(sigma^2)
  tmp2=sum(log(eigen(V0)$value));
  tmp3=log(det(XTV0invX))
  tmp4= t(Y) %*% P %*% Y/sigma^2;
  like=tmp1+tmp2+tmp3+tmp4;
  like=-1/2*like;
  return(like);
}
firstderisigmaF<-function(parameter,X,Y,K)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)];#dr=round(dr,8);
  sigma=parameter[length(K)+1];
  V0=diag(n);dr=round(dr,14);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  firstderi=NULL;
  for(i in 1:length(K))
  {
    tmp=-1/2*(tr(P %*% K[[i]] *2 *dr[i])-1/sigma^2 * t(Y) %*% P %*% (K[[i]]*2*dr[i]) %*% P %*% Y)
    firstderi=c(firstderi,tmp)
  }
  firstderi=c(firstderi,(-(n-p)/2/sigma^2+t(Y) %*% P %*% Y/2/sigma^4)*(2*sigma));
  firstderi;
}

ExactNegploglikeF<-function(parameter,parameter0,K,Y,X,weight,lambda=0)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)];#dr=round(dr,8);
  dr0=parameter0[1:length(K)];#dr0=round(dr0,8);
  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  tmp1=(n-p)*log(sigma^2)
  tmp2=sum(log(abs(eigen(V0)$value)))
  tmp3=sum(log(abs(eigen(XTV0invX)$value)))
  tmp4= t(Y) %*% P %*% Y/sigma^2;
  like=tmp1+tmp2+tmp3+tmp4;
  #plike=-1/2*like-lambda*sum(weight*dr^2/2/abs(dr0));
  plike=-1/2*like-lambda*sum(weight*abs(dr));
  return(-plike);
}

ExactNegpfirstderisigmaF<-function(parameter,parameter0,X,Y,K,weight,lambda=0)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  dr0=parameter0[1:length(K)]; #dr0=round(dr0,8);

  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  firstderi=NULL;
  for(i in 1:length(K))
  {
    tmp=-1/2*(tr(P %*% K[[i]] *2 *dr[i])-1/sigma^2 * t(Y) %*% P %*% (K[[i]]*2*dr[i]) %*% P %*% Y)
    firstderi=c(firstderi,tmp)
  }
  #firstderi=firstderi-lambda*weight*dr/abs(dr0)
  firstderi=firstderi-lambda*weight*sign(dr);
  firstderi=c(firstderi,(-(n-p)/2/sigma^2+t(Y) %*% P %*% Y/2/sigma^4)*(2*sigma));
  -firstderi;
}




InforsigmaF<-function(parameter,X,Y,K)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  sigma=parameter[length(K)+1];

  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  obsInfor=matrix(NA,length(K)+1,length(K)+1);
  expInfor=matrix(NA,length(K)+1,length(K)+1);
  exactAI=matrix(NA,length(K)+1,length(K)+1);
  approxAI=matrix(NA,length(K)+1,length(K)+1);
  # For Sigma related terms #
  obsInfor[length(K)+1,length(K)+1]=(-(n-p)/sigma^2+t(Y) %*% P %*% Y/sigma^4*3)
  expInfor[length(K)+1,length(K)+1]=(2*(n-p)/sigma^2)
  exactAI[length(K)+1,length(K)+1]=(obsInfor[length(K)+1,length(K)+1]+expInfor[length(K)+1,length(K)+1])/2;
  approxAI[length(K)+1,length(K)+1]=(obsInfor[length(K)+1,length(K)+1]+expInfor[length(K)+1,length(K)+1])/2;
  for(j in 1:length(K)){
    obsInfor[length(K)+1,j]=obsInfor[j,length(K)+1]=(t(Y) %*% P %*% (K[[j]]*2*dr[j]) %*% P %*% Y/2/sigma^4)*(2*sigma);
    expInfor[length(K)+1,j]=expInfor[j,length(K)+1]=(1/2/sigma^2/tr(P %*% (K[[j]]*2*dr[j])))*(2*sigma);
    exactAI[length(K)+1,j]=exactAI[j,length(K)+1]=(obsInfor[length(K)+1,j]+expInfor[length(K)+1,j])/2;
    approxAI[length(K)+1,j]=approxAI[j,length(K)+1]=obsInfor[length(K)+1,j]
  }
  for(i in 1:length(K))
    for(j in i:length(K))
    {
      tmp1=K[[i]] * dr[i]*2;
      tmp2=K[[j]] * dr[j]*2;
      obsInfor[i,j]=obsInfor[j,i]=-1/2*tr(P %*% tmp1 %*% P %*% tmp2)+1/sigma^2 * t(Y) %*% P %*% tmp1 %*% P %*% tmp2 %*% P %*% Y   ;
      expInfor[i,j]=expInfor[j,i]=-1/2*tr(P %*% tmp1 %*% P %*% tmp2);
      if(i==j)
      {
        obsInfor[i,i]=obsInfor[i,i]+tr(P %*% K[[i]])-1/sigma^2 * t(Y) %*% P %*% K[[i]] %*% P %*% Y;
        expInfor[i,i]=expInfor[i,i]+tr(P %*% K[[i]])-1/sigma^2 * t(Y) %*% P %*% K[[i]] %*% P %*% Y
      }

      exactAI[i,j]=exactAI[j,i]=(obsInfor[i,j]+expInfor[i,j])/2;
      approxAI[i,j]=approxAI[j,i]=(obsInfor[i,j]+expInfor[i,j])/2;
    }
  result=list();
  result$obsInfor=obsInfor;
  result$expInfor=expInfor;
  result$exactAI=exactAI;
  result$approxAI=approxAI;
  result
}
NRpInforsigmaF<-function(parameter,X,Y,K, parameter0, weight, lambda)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  dr0=parameter0[1:length(K)];#dr0=round(dr0,8);
  sigma=parameter[length(K)+1];

  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  obsInfor=matrix(NA,length(K)+1,length(K)+1);
  expInfor=matrix(NA,length(K)+1,length(K)+1);
  exactAI=matrix(NA,length(K)+1,length(K)+1);
  approxAI=matrix(NA,length(K)+1,length(K)+1);
  PY=P %*% Y;
  # For Sigma related terms #
  obsInfor[length(K)+1,length(K)+1]=(-(n-p)/sigma^2+t(Y) %*% PY/sigma^4*3)
  expInfor[length(K)+1,length(K)+1]=(2*(n-p)/sigma^2)
  exactAI[length(K)+1,length(K)+1]=(obsInfor[length(K)+1,length(K)+1]+expInfor[length(K)+1,length(K)+1])/2;
  approxAI[length(K)+1,length(K)+1]=(obsInfor[length(K)+1,length(K)+1]+expInfor[length(K)+1,length(K)+1])/2;
  PK=PKsigma=list();
  for(j in 1:length(K)){
    #PKsigma[[j]]=P %*% (K[[j]]*2*dr[j]);
    PK[[j]]=P %*% K[[j]];
    obsInfor[length(K)+1,j]=obsInfor[j,length(K)+1]=(t(Y) %*% (PK[[j]]*2*dr[j]) %*% PY/2/sigma^4)*(2*sigma);
    expInfor[length(K)+1,j]=expInfor[j,length(K)+1]=(1/2/sigma^2/tr((PK[[j]]*2*dr[j])))*(2*sigma);
    exactAI[length(K)+1,j]=exactAI[j,length(K)+1]=(obsInfor[length(K)+1,j]+expInfor[length(K)+1,j])/2;
    approxAI[length(K)+1,j]=approxAI[j,length(K)+1]=obsInfor[length(K)+1,j]
  }
  for(i in 1:length(K))
    for(j in i:length(K))
    {
      #tmp1=K[[i]] * dr[i]*2;
      #tmp2=K[[j]] * dr[j]*2;
      tmp=PK[[i]] %*% PK[[j]]*2*dr[i]*2*dr[j]
      obsInfor[i,j]=obsInfor[j,i]=-1/2*tr(tmp)+1/sigma^2 * t(Y) %*% tmp %*% PY  ;
      expInfor[i,j]=expInfor[j,i]=-1/2*tr(tmp);
      if(i==j)
      {
        obsInfor[i,i]=obsInfor[i,i]+tr(PK[[i]])-1/sigma^2 * t(Y) %*% PK[[i]] %*% PY;
        expInfor[i,i]=expInfor[i,i]+tr(PK[[i]])-1/sigma^2 * t(Y) %*% PK[[i]] %*% PY;
      }

      exactAI[i,j]=exactAI[j,i]=(obsInfor[i,j]+expInfor[i,j])/2;
      approxAI[i,j]=approxAI[j,i]=(obsInfor[i,j]+expInfor[i,j])/2;
    }
  negsecondpenalty=matrix(0,nrow=nrow(exactAI),ncol=ncol(exactAI))
  diag(negsecondpenalty)[1:length(K)]=lambda*weight/abs(dr0);
  result=list();

  result$obsInfor=obsInfor+negsecondpenalty;
  result$expInfor=expInfor+negsecondpenalty;
  result$exactAI=exactAI+negsecondpenalty;
  result$approxAI=approxAI+negsecondpenalty;
  result
}

NRNegploglikeF<-function(parameter,parameter0,K,Y,X,weight,lambda=0)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)];#dr=round(dr,8);
  dr0=parameter0[1:length(K)];#dr0=round(dr0,8);
  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  tmp1=(n-p)*log(sigma^2)
  tmp2=sum(log(abs(eigen(V0)$value)))
  tmp3=sum(log(abs(eigen(XTV0invX)$value)))
  tmp4= t(Y) %*% P %*% Y/sigma^2;
  like=tmp1+tmp2+tmp3+tmp4;
  plike=-1/2*like-lambda*sum(weight*dr^2/2/abs(dr0));
  #plike=-1/2*like-lambda*sum(weight*abs(dr));
  return(-plike);
}
NRNegpfirstderisigmaF<-function(parameter,parameter0,X,Y,K,weight,lambda=0)
{
  n=nrow(X);p=ncol(X);
  dr=parameter[1:length(K)]; #dr=round(dr,8);
  dr0=parameter0[1:length(K)]; #dr0=round(dr0,8);

  sigma=parameter[length(K)+1];
  V0=diag(n);
  for(i in 1:length(K)) V0=V0+dr[i]^2*K[[i]]
  V0inv=MASS::ginv(V0);
  XTV0invX=t(X) %*% V0inv %*% X
  P=V0inv-V0inv %*% X %*% MASS::ginv(XTV0invX) %*% t(X) %*% V0inv;
  firstderi=NULL;
  for(i in 1:length(K))
  {
    tmp=-1/2*(tr(P %*% K[[i]] *2 *dr[i])-1/sigma^2 * t(Y) %*% P %*% (K[[i]]*2*dr[i]) %*% P %*% Y)
    firstderi=c(firstderi,tmp)
  }
  firstderi=firstderi-lambda*weight*dr/abs(dr0)
  #firstderi=firstderi-lambda*weight*sign(dr);
  firstderi=c(firstderi,(-(n-p)/2/sigma^2+t(Y) %*% P %*% Y/2/sigma^4)*(2*sigma));
  -firstderi;
}

NRoptimStep<-function(parameter,parameter0,Y,K,X,lambda,weight,maxiter=1000, NRtol=1e-6)
{
  run=TRUE; iter=0;
  parameternew=parameter;
  llnew=NRNegploglikeF(parameter=parameter,parameter0=parameter0,K=K,Y=Y,X=X,weight=weight,lambda=lambda);
  while(run & iter< maxiter/100)
  {
    #cat("inner iter", iter)
    #print(parameter)
    iter=iter+1;
    llold=llnew;
    llnew=llold+1;
    parameter=parameternew;
    mm=1;
    la=1;
    firstderi= -NRNegpfirstderisigmaF(parameter=parameter,parameter0=parameter0,K=K,Y=Y,X=X,weight=weight,lambda=lambda)

    infor=NRpInforsigmaF(parameter=parameter,parameter0=parameter0,K=K,Y=Y,X=X,weight=weight,lambda=lambda)$obsInfor
    gm = MASS::ginv(infor)%*%firstderi;
    while(llnew>=llold&&mm<25){
      parameternew=parameter+la*gm;
      llnew <- NRNegploglikeF(parameter=parameternew,parameter0=parameter0,K=K,Y=Y,X=X,weight=weight,lambda=lambda)
      la <- 1/2^mm
      mm <- mm+1
    }
    #cat("dif=",sum((parameternew-parameter)^2),"\n")
    if(mean((parameternew-parameter)^2)<NRtol) run = FALSE
  }
  parameter;
}

NRoptim<-function(parameter,parameter0,Y,K,X,lambda,weight,maxiter=1000,NRtol=1e-6, eps=1e-6)
{
  firstderi= -NRNegpfirstderisigmaF(parameter=parameter,parameter0=parameter0,K=K,Y=Y,X=X,weight=weight,lambda=lambda)
  if(max(abs(firstderi/parameter))>1e3)
  {
    bb4=optim(par=parameter, fn=NRNegploglikeF, gr=NRNegpfirstderisigmaF,K=K,X=X,Y=Y,
              weight=weight,parameter0=parameter0,lambda=lambda,
              control=list(maxit=maxiter), method="L-BFGS-B")
    parameter=parameter0=bb4$par
  }

  parameterfinal=rep(NA,length(parameter));delete=rep(0,length(parameter)-1);
  iter=0;run=TRUE;
  parameterS=parameter;weightnew=weight; Knew=K;
  while(run & iter<maxiter & (sum(delete==0)!=0))
  {
    #print("iter=",iter)
    #print(parameterS)
    iter=iter+1;
    parameterSteps=NRoptimStep(parameter=parameterS,parameter0=parameterS,Y=Y,K=Knew,X=X,lambda=lambda,weight=weightnew,maxiter=maxiter, NRtol=NRtol*10)
    if(sum((parameterSteps-parameterS)^2)<NRtol) run <- FALSE;
    parameterS=parameterSteps;
    drtmp=parameterS[-length(parameterS)]
    if(sum(abs(drtmp)<eps)>0)
    {
      delete[delete==0][abs(drtmp)<eps]=1;
      index=which(abs(drtmp)<eps)

      drtmp=drtmp[!abs(drtmp)<eps];
      weightnew=weight[delete==0];
      index=index[order(index)]
      for(i in 1:length(index))
      {
        Knew[[index[i]]]=NULL;
        index=index-1;
      }
    }
    parameterS=c(drtmp,parameterSteps[length(parameterSteps)])
  }
  delete=c(delete,0)
  parameterfinal[delete==1]=0;
  parameterfinal[delete==0]=parameterS;
  parameterfinal;
}


PLME<-function(K,Y,X=NULL,maxiter=1000,minheri=0.01,lambdarange=c(0,100), tol=1e-6,crit='bic',train.id=NULL, predict=FALSE, weight.fixed=NULL, weight.random=NULL)
{
#  cat("PLME",predict,"-",!is.null(train.id))
  n=length(Y)
  X=cbind(matrix(1,nrow=n,ncol=1),X);
  if(nrow(X)!=n) stop("Error: covariates and the outcomes are not of the same dimension!")
  nntmp=NULL;
  for(i in 1:length(K)) nntmp=c(nntmp,nrow(K[[i]]))
  if(sum(nntmp!=n)!=0) stop("Error: genomic effects (i.e. random effect) and/or the outcomes are not of the same dimension!")
  Ypred=NULL;
  if(!is.null(train.id))
  {
    if(max(train.id)>n | length(train.id) > n) stop("Error: training sample ids are not correct!")
    Ktr=list();for(i in 1:length(K)) Ktr[[i]]=K[[i]][train.id,train.id]
    Xtr=X[train.id,]; Ytr=Y[train.id]
    if(is.null(dim(Xtr))) Xtr=matrix(Xtr,ncol=1)

    fit=PLMEfit(K=Ktr,Y=Ytr,X=Xtr,maxiter=maxiter,minheri=minheri,lambdarange=lambdarange, tol=tol,crit=crit,weight.fixed=weight.fixed,weight.random=weight.random)
    if(predict)
    {
      print("Calculating predicted values......")
      dr=fit$random$dr;
      sigma=fit$random$sigma;
      beta=fit$fixed$fixed;
      Ypred=PLMEpredict(K=K,Y=Y,X=X,dr=dr,sigma=sigma,beta=beta,train.id=train.id)
    }
  }
  if(is.null(train.id))
    fit=PLMEfit(K=K,Y=Y,X=X,maxiter=maxiter,minheri=minheri,lambdarange=lambdarange, tol=tol,crit=crit,weight.fixed=weight.fixed,weight.random=weight.random)
  result=list();
  result$fit=fit;
  result$Ypred=Ypred;
  result

}


PLMEfit<-function(K,Y,X,maxiter=1000,minheri=0.01,lambdarange=c(0,100), tol=1e-6,crit='bic',weight.fixed=NULL, weight.random=NULL)
{
  if(ncol(X)<nrow(X)) random=PLME.Random(K=K,Y=Y,X=X,maxiter=maxiter,minheri=minheri,lambdarange=lambdarange,weight=weight.random)
  if(ncol(X)>=nrow(X))
  {
    X1=X[,-1]
    lambda=glmnet::cv.glmnet(y=Y,x=X1)$lambda.min;
    l1=glmnet::glmnet(y=Y,x=X1,lambda = lambda)
    X1=X[,which(as.matrix(coef(l1))!=0)];
    if(is.null(dim(X1))) X1=X[,1:(nrow(X)/2)]
    if(ncol(X1)>nrow(X1) & !is.null(dim(X1))) X1=X1[,1:(nrow(X1)/2)]
    print("HERE1")

    random=PLME.Random(K=K,Y=Y,X=X1,maxiter=maxiter,minheri=minheri,lambdarange=lambdarange)
  }
  dr=random$dr;
  #print(dr[dr!=0])
  sigma=random$sigma;

  fixed=PLME.Fixed(K=K,Y=Y,X=X,dr=dr,sigma=sigma,weight=weight.fixed);
  result=list();
  result$fixed=fixed;
  result$random=random;
  result;
}




PLME.Fixed<-function(K,Y,X,dr,sigma,maxiter=1000,tol=1e-10,crit='bic',weight=NULL)
{
  set.seed(10);
  print("Getting fixed effects!")
  n=nrow(X);p=ncol(X);
  if(is.null(weight)) weight=rep(1,p);
  if(length(weight)!=p) stop("weight should have the same length as the number of covariates")
  V=diag(n); for(i in 1:length(K)) V=V+dr[i]^2*K[[i]]; V=V*sigma^2;
  Vinv=MASS::ginv(V);

  l=round(KFAS::ldl(Vinv, tol = tol),-log10(tol))
  d <- diag(sqrt(diag(l)))
  diag(l) <- 1
  tmp=t(l %*% d);
  Vinvsqrt=round(tmp,-log10(tol))
  Ynew=Vinvsqrt %*% Y;
  Xnew=Vinvsqrt %*% X;
  if(ncol(Xnew)<nrow(Xnew)) beta0=MASS::ginv(t(Xnew) %*% Xnew) %*% t(Xnew) %*% Ynew
  if(ncol(Xnew)>=nrow(Xnew))
  {
    if(sum(dr==0)==length(dr)){
      lasso_cv <- glmnet::cv.glmnet(x = X[,-1], y = Y,type.measure = "mse",nfold = 10,alpha = 1,intercept=TRUE,standardize=FALSE);
      beta0 <- as.numeric(coef(lasso_cv, s = lasso_cv$lambda.min))
    }
    if(sum(dr==0)!=length(dr)){
      lasso_cv <- glmnet::cv.glmnet(x = Xnew, y = Ynew,type.measure = "mse",nfold = 10,alpha = 1,intercept=FALSE,standardize=FALSE);
      beta0 <- as.numeric(coef(lasso_cv, s = lasso_cv$lambda.min))[-1];
    }
    beta0[abs(beta0)<1e-6]=1e-6
  }

  weight=weight*1/abs(beta0);
  weight[1]=0;
  #print("Getting fixed effects: run!")

  result=list();
  result$lambda=0;
  result$crit=crit;
  result$fixed=beta0;
  if(ncol(Xnew)>1) {
    if(sum(dr==0)==length(dr))
    {
      fit=HDeconometrics::ic.glmnet(x = X[,-1], y = Y,crit=crit,penalty.factor = weight[-1],standardize = FALSE, intercept=TRUE)
      result$fixed=fit$coefficients
    }
    if(sum(dr==0)!=length(dr))
    {
      fit=HDeconometrics::ic.glmnet(x = Xnew, y = Ynew,crit=crit,penalty.factor = weight,standardize = FALSE, intercept=FALSE)
      result$fixed=fit$coefficients[-1];
    }
    result$lambda=fit$lambda
  }
  result
}

PLME.Random<-function(K,Y,X,maxiter=1000,minheri=0.01,lambdarange=c(0,100),weight=NULL)
{
  #print(length(weight))

  set.seed(10);

  # Initilization #
  print("Getting random effects!")
  n=nrow(X); p=ncol(X);
  beta=MASS::ginv(t(X) %*% X) %*% t(X) %*% Y;
  e=Y-X %*% beta;
  sigma=sqrt(t(e) %*% e/(n-p));minherisave=minheri;
  if(0.5/length(K)<=minheri) minheri=0.5/length(K)*0.9;
  tmp=0.5+length(K)*minheri;tmp1=sqrt(1/length(K)*tmp/(1-tmp));
  minheri=minherisave;
  dr=rep(tmp1,length(K)) # genetics total explains 50% heri
  parameter=c(dr,sigma)
  # initial estimates #
  #print("Getting MLEs")
  bb=optim(parameter, fn=negloglikeF, gr=negfirstderisigmaF, K=K,X=X,Y=Y,
           control=list(maxit=maxiter),  method="L-BFGS-B");
  # print(bb$par[1:10]);print(bb$convergence)
  if(bb$convergence!=0) bb=optim(parameter, fn=negloglikeF, gr=negfirstderisigmaF, K=K,X=X,Y=Y,
                                 control=list(maxit=maxiter),  method="L-BFGS-B")
  initialsave=bb;

  parameterfinal=parameter0=bb$par
  if(is.null(weight)) weight=rep(1,length(parameter0[1:length(K)]))

  weight=1/abs(parameter0[1:length(K)])*weight;
  #cat("MLE",tail(parameterfinal),"\n")
  # penalized estimates #
  #print("Getting PLEs:")
  lambdarange=lambdarange[order(lambdarange)]
  seqdif=(lambdarange[2]-lambdarange[1])/10;lambdaall=seq(from=lambdarange[1],to=lambdarange[2],by=seqdif);
  BICs=NULL;
  for(i in 1:length(lambdaall))
  {
    lambda=lambdaall[i];
    exc=parameterfinal[-length(parameterfinal)]^2/(sum(parameterfinal[-length(parameterfinal)]^2)+1) < minheri;
    index=which(!exc);
    Knew=list();
    if(length(index)==0) {parameterfinal[which(exc)]=0;}
    if(length(index)>0)
    {
      parameterfinal[which(exc)]=0
      parameter=c(parameterfinal[index],parameterfinal[length(parameterfinal)]);
      weightnew=weight[index]
      #print(tail(weightnew))

      for(j in 1:length(index)) Knew[[j]]=K[[index[j]]]
      parameter=round(parameter,8);
      #bb1=optim(parameter, fn=ExactNegploglikeF, gr=ExactNegpfirstderisigmaF, K=Knew,X=X,Y=Y,
      #          weight=weightnew,parameter0=parameter,lambda=lambda,
      #          control=list(maxit=maxiter), method="BFGS")
      bb1=NRoptim(parameter=parameter,parameter0=parameter,Y=Y,K=Knew,X=X,lambda=lambda,weight=weightnew,maxiter=maxiter,NRtol=1e-6, eps=1e-6)
      parameterfinalnew=parameterfinal;
      #parameterfinalnew[parameterfinal!=0]=bb1$par
      parameterfinalnew[parameterfinal!=0]=bb1
      parameterfinal=parameterfinalnew;
    }
    bics=2*negloglikeF(parameterfinal,K=K,Y=Y,X=X)+log(length(Y))*(sum(parameterfinal!=0)-1)
    aa=c(lambda,bics,parameterfinal,0)
    BICs=rbind(BICs,aa)
    if(sum(parameterfinal[-length(parameterfinal)]!=0)==0) break;
  }
  colnames(BICs)=c("lambda","BICs",paste("d",1:(length(parameterfinal)-1),sep=""),"sigma","converge")

  #print("Getting PLEs: choosing lambda")
  # finetune lambda #
  select=which(BICs[,"BICs"]==min(BICs[,"BICs"]))[1]
  if(sum(BICs[,"converge"]==0)>0) select=which(BICs[BICs[,"converge"]==0,"BICs"]==min(BICs[BICs[,"converge"]==0,"BICs"]))[1]
  if(nrow(BICs)==1) {BICsfinal=BICs[1,]}
  if(nrow(BICs)>1)
  {
    ## Setting range of new lambdas ##
    {
      if(select==1) lambdaall=seq(from=BICs[select,"lambda"], to=BICs[select+1,"lambda"],by=seqdif/20);
      if(select==nrow(BICs))
      {
        tmp=BICs[select,3:(length(parameterfinal)+1)];
        if(sum(tmp!=0)==0) lambdaall=seq(from=BICs[select-1,"lambda"], to=BICs[select,"lambda"],by=seqdif/20)
        if(sum(tmp!=0)!=0) lambdaall=seq(from=BICs[select-1,"lambda"]+seqdif/2, to=BICs[select,"lambda"]+seqdif/2,by=seqdif/20)
      }
      if(select!=1 & select!=nrow(BICs)) lambdaall=seq(from=BICs[select-1,"lambda"]+seqdif/2,to=BICs[select+1,"lambda"]-seqdif/2,by = seqdif/20)
    }
    ## Rerun under the new lambdas ##
    BICs=NULL;      parameterfinal=parameter0=bb$par
    for(i in 1:length(lambdaall))
    {
      lambda=lambdaall[i];
      exc=parameterfinal[-length(parameterfinal)]^2/(sum(parameterfinal[-length(parameterfinal)]^2)+1) < minheri;
      index=which(!exc);
      Knew=list();
      if(length(index)==0) {parameterfinal[which(exc)]=0;}
      if(length(index)>0)
      {
        parameterfinal[which(exc)]=0
        parameter=c(parameterfinal[index],parameterfinal[length(parameterfinal)]);
        weightnew=weight[index]
        for(j in 1:length(index)) Knew[[j]]=K[[index[j]]]
        parameter=round(parameter,8);
        #bb1=optim(parameter, fn=ExactNegploglikeF, gr=ExactNegpfirstderisigmaF, K=Knew,X=X,Y=Y,
        #          weight=weightnew,parameter0=parameter,lambda=lambda,
        #          control=list(maxit=maxiter), method="BFGS")
        bb1=NRoptim(parameter=parameter,parameter0=parameter,Y=Y,K=Knew,X=X,lambda=lambda,weight=weightnew,maxiter=maxiter,NRtol=1e-6, eps=1e-6)

        parameterfinalnew=parameterfinal;
        #parameterfinalnew[parameterfinal!=0]=bb1$par
        parameterfinalnew[parameterfinal!=0]=bb1
        parameterfinal=parameterfinalnew;
      }
      bics=2*negloglikeF(parameterfinal,K=K,Y=Y,X=X)+log(length(Y))*(sum(parameterfinal!=0)-1)
      aa=c(lambda,bics,parameterfinal,0)
      BICs=rbind(BICs,aa)
      if(sum(parameterfinal[-length(parameterfinal)]!=0)==0) break;
    }
  }

  colnames(BICs)=c("lambda","BICs",paste("d",1:(length(parameterfinal)-1),sep=""),"sigma","converge")
  select=which(BICs[,"BICs"]==min(BICs[,"BICs"]))[1];convs=FALSE;
  if(sum(BICs[,"converge"]==0)>0)
  {
    select=which(BICs[BICs[,"converge"]==0,"BICs"]==min(BICs[BICs[,"converge"]==0,"BICs"]))[1];
    convs=TRUE;
  }
  BICsfinal=BICs[select,];
  names(BICsfinal)=c("lambda","BICs",paste("d",1:(length(parameterfinal)-1),sep=""),"sigma","converge")

  ## Summarize Results ##
  result=list();
  result$lambda=BICsfinal["lambda"];
  result$BIC=BICsfinal["BICs"];
  result$dr=BICsfinal[grep("d",names(BICsfinal))[-1]];
  result$sigma=BICsfinal["sigma"];
  result$e=e
  result
}


PLMEpredict<-function(K,Y,X,dr,sigma,beta,train.id)
{
  nall=nrow(K[[1]]);
  V=diag(nall);for(i in 1:length(K)) V=V+dr[i]^2*K[[i]];
  V=V*sigma^2;
  if(ncol(X)!=length(beta)) stop("the fixed dimension does not match!")
  Ypred=X[-train.id,] %*% beta + V[-train.id,train.id] %*% MASS::ginv(V[train.id,train.id]) %*% (Y[train.id]-X[train.id,]%*%beta)
  Ypred
}


OmicsPLMMPred<-function(Data,predict=FALSE, weight.fixed=NULL,weight.random=NULL,maxiter=1000,minheri=0.01,lambdarange=c(0,100), tol=1e-6,crit='bic',outputall=0)
{
  #print("PRED")
  #print(predict)
  outcome=Data$Y;
  train.id=which(names(Y) %in% Data$trainID);
  Cov=Data$X;
  if(!is.null(Cov))
  {
    if(ncol(Cov)==1) Cov=matrix(scale(Cov),ncol=1);
    if(ncol(Cov)>1) Cov=apply(Cov,2,scale);
  }
  XCov=Cov;
  if(!is.null(XCov))
  {
    if(length(weight.fixed)!=ncol(XCov) & (!is.null(weight.fixed))){
      warning("fixed weight does not have the same length as covaraites, and thus set to 1")
      weight.fixed=NULL;
    }
    Xcovtmp=apply(XCov[train.id,],2,sd)!=0
    if(sum(Xcovtmp)>0) XCov=XCov[,Xcovtmp];
    if(!is.null(weight.fixed) ) weight.fixed=weight.fixed[c(TRUE,Xcovtmp)]
    if(sum(Xcovtmp)==0) {XCov=NULL;weight.fixed=NULL}
  }

  # create similarity and remove redundent ones #
  K=list();start=1;
  for(i in 1:length(Data$KernelOutput))
  {
    if(length(Data$KernelOutput[[i]]$Kinship)>0)
    {
      for(j in 1:length(Data$KernelOutput[[i]]$Kinship))
      {
        if(length(Data$KernelOutput[[i]]$Kinship[[j]])>0){
          add=TRUE;
          if(start==1){
            K[[start]]=Data$KernelOutput[[i]]$Kinship[[j]];
            names(K)[start]=paste(i,Data$KernelOutput[[i]]$IncludeRegions[j],Data$KernelOutput[[i]]$IncludeRegionsNamesAll[j],sep="_")
            start=start+1;
          }
          if(start>1){
            for(k in 1:length(K)){
              if(all.equal(K[[k]],Data$KernelOutput[[i]]$Kinship[[j]])==TRUE) {add=FALSE;break;}
            }
            if(add){
              K[[start]]=Data$KernelOutput[[i]]$Kinship[[j]];
              names(K)[start]=paste(i,Data$KernelOutput[[i]]$IncludeRegions[j],Data$KernelOutput[[i]]$IncludeRegionsNamesAll[j],sep="_")
              start=start+1;
            }
            if(!add){
              names(K)[k]=paste(paste(names(K)[k],":",sep=""),i,Data$KernelOutput[[i]]$IncludeRegions[j],Data$KernelOutput[[i]]$IncludeRegionsNamesAll[j],sep="_");
            }
          }
        }

      }
    }
  }

  print("Fitting models")
  result=PLME(K=K,Y=outcome,X=XCov,maxiter=maxiter,minheri=minheri,lambdarange=lambdarange, tol=tol,crit=crit,weight.fixed=weight.fixed,weight.random=weight.random,train.id=train.id, predict=predict)
  results=list();
  results$Ypred=result$Ypred;
  results$outcome=outcome;
  if(outputall>0) results$fit=result$fit;
  if(outputall==2) results$Similarity=K;
  result=results;
  result;
}

