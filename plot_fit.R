#function to generate predicted response with confidence intervals from a (G)LM(M)
#works with the following model/class: lm, glm, glm.nb, merMod
#this function average over potential covariates
#it also allows for the specification of one or several interacting variables
#these must be factor variables in the model
#for (G)LMM the name of the random terms must be specfied in the RE argument
#for (G)LMM the confidence interval can be either bootstrapped or coming from
#a normal approximation using

#list of arguments:
#@m: a model object, either of class: lm, glm, glm.nb, merMod
#@focal_var: a character, the name of the focal variable that will be on the x-axis
#@inter_var: a character or a character vector, the names(s) of the interacting variables, must be declared as factor variables in the model, default is NULL
#@RE: a charcater or a charcater vector, the name(s) of the random effect variables in the case of a merMod object, default is NULL
#@offset: a character, the name of the offset variable, note that this effect will be averaged out like other continuous covariates, this is maybe not desirable
#@n: an integer, the number of generated prediction points, default is 20
#@n_core: an integer, the number of cores to use in parallel computing for the bootstrapped CI for merMod object, default is 4
#@boot_mer: a logical, whether to use bootstrapped (TRUE) or a normal approximation (FALSE, the default) for the confidence interval in the case of a merMod model 

plot_fit<-function(m,focal_var,inter_var=NULL,RE=NULL,offset=NULL,n=20,n_core=4,boot_mer=FALSE){
  require(arm)  
  dat<-model.frame(m)
  #turn all character variable to factor
  dat<-as.data.frame(lapply(dat,function(x){
    if(is.character(x)){
      as.factor(x)
    }
    else{x}
  }))
  #make a sequence from the focal variable
  x1<-list(seq(min(dat[,focal_var]),max(dat[,focal_var]),length=n))
  #grab the names and unique values of the interacting variables
  isInter<-which(names(dat)%in%inter_var)
  if(length(isInter)==1){
    x2<-list(unique(dat[,isInter]))
    names(x2)<-inter_var
  }
  if(length(isInter)>1){
    x2<-lapply(dat[,isInter],unique)
  }
  if(length(isInter)==0){
    x2<-NULL
  }
  #all_var<-x1
  #add the focal variable to this list
  all_var<-c(x1,x2)
  #expand.grid on it
  names(all_var)[1]<-focal_var
  all_var<-expand.grid(all_var)
  
  #remove varying variables and non-predictors and potentially offset variables
  off_name <- NULL
  if(!is.null(offset)){
    off_name <- grep("^offset",names(dat),value=TRUE)#this is needed because of the weird offset formatting in the model.frame
  }
  dat_red<-dat[,-c(1,which(names(dat)%in%c(focal_var,inter_var,RE,"X.weights.",off_name))),drop=FALSE]
  #if there are no variables left over that need averaging
  if(dim(dat_red)[2]==0){
    new_dat<-all_var
    name_f <- NULL
  }
  else{
    #otherwise add these extra variables, numeric variable will take their mean values
    #and factor variables will take their first level before being averaged out lines 86-87
    fixed<-lapply(dat_red,function(x) if(is.numeric(x)) mean(x) else factor(levels(x)[1],levels = levels(x)))
    #the number of rows in the new_dat frame
    fixed<-lapply(fixed,rep,dim(all_var)[1])
    #create the new_dat frame starting with the varying focal variable and potential interactions
    new_dat<-cbind(all_var,as.data.frame(fixed)) 
    #get the name of the variable to average over
    name_f<-names(dat_red)[sapply(dat_red,function(x) ifelse(is.factor(x),TRUE,FALSE))]
  }  
  #add an offset column set at 0 if needed
  if(!is.null(offset)){
    new_dat[,offset] <- 0
  }
  
  
  #get the predicted values
  cl<-class(m)[1]
  if(cl=="lm"){
    pred<-predict(m,newdata = new_dat,se.fit=TRUE)
  }
  
  if(cl=="glm" | cl=="negbin"){
    #predicted values on the link scale
    pred<-predict(m,newdata=new_dat,type="link",se.fit=TRUE)
  }
  if(cl=="glmerMod" | cl=="lmerMod"){
    pred<-list(fit=predict(m,newdata=new_dat,type="link",re.form=~0))
    #for bootstrapped CI
    new_dat<-cbind(new_dat,rep(0,dim(new_dat)[1]))
    names(new_dat)[dim(new_dat)[2]]<-as.character(formula(m)[[2]])
    mm<-model.matrix(formula(m,fixed.only=TRUE),new_dat)
  }
  #average over potential categorical variables  
  avg_over <- 0 #for cases where no averaging is to be done
  if(length(name_f)>0){
    if(cl=="glmerMod" | cl=="lmerMod"){
      coef_f<-lapply(name_f,function(x) fixef(m)[grep(paste0("^",x),names(fixef(m)))])
    }
    else{
      coef_f<-lapply(name_f,function(x) coef(m)[grep(paste0("^",x),names(coef(m)))])
    }    
    avg_over <- sum(unlist(lapply(coef_f,function(x) mean(c(0,x))))) #averging out all factor effects
    pred$fit<-pred$fit + avg_over
  }
  
  #to get the back-transform values get the inverse link function
  linkinv<-family(m)$linkinv
  
  #get the back transformed prediction together with the 95% CI for LM and GLM
  if(cl=="glm" | cl=="lm" | cl=="negbin"){
    pred$pred<-linkinv(pred$fit)
    pred$LC<-linkinv(pred$fit-1.96*pred$se.fit)
    pred$UC<-linkinv(pred$fit+1.96*pred$se.fit)
  }
  
  #for (G)LMM need to use either bootstrapped CI or use approximate
  #standard error from the variance-covariance matrix
  #see ?predict.merMod and http://glmm.wikidot.com/faq#predconf
  #note that the bootstrapped option is recommended by the lme4 authors
  if(cl=="glmerMod" | cl=="lmerMod"){
    pred$pred<-linkinv(pred$fit)
    if(boot_mer){
      predFun<-function(.) mm%*%fixef(.)+avg_over
      bb<-bootMer(m,FUN=predFun,nsim=200,parallel="multicore",ncpus=n_core) #do this 200 times
      bb$t<-apply(bb$t,1,function(x) linkinv(x))
      #as we did this 200 times the 95% CI will be bordered by the 5th and 195th value
      bb_se<-apply(bb$t,1,function(x) x[order(x)][c(5,195)])
      pred$LC<-bb_se[1,]
      pred$UC<-bb_se[2,] 
    }
    else{
      se <- diag(mm %*% tcrossprod(vcov(m),mm))
      pred$LC <- linkinv(pred$fit - 1.96 * sqrt(se))
      pred$UC <- linkinv(pred$fit + 1.96 * sqrt(se))
    }
  }
  
  #the output
  out<-as.data.frame(cbind(new_dat[,1:(length(inter_var)+1)],pred$LC,pred$pred,pred$UC))
  names(out)<-c(names(new_dat)[1:(length(inter_var)+1)],"LC","Pred","UC")
  return(out)
}