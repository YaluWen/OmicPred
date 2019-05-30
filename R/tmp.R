CalculateSimilarity<-function(OmicsDataEach,OmicsDataEachMap,Annotation,kernels,AllRegions,Genomic="FALSE")
{
  if(nrow(OmicsDataEachMap)!=ncol(OmicsDataEach)) stop("Error: The number of variants included and their annotations do not match")
  Kinship=list();IncludeRegions=NULL; index=1;
  print("Calculating Similarities!")
  if(AllRegions) {Kinship[[1]]=getKernel(Input=OmicsDataEach,kernels = kernels, Genomic = Genomic, AllRegions=AllRegions); IncludeRegions=rbind(IncludeRegions,"chr1:chr22:all:all");index=index+1}
  for(i in 1:nrow(Annotation))
  {
    tmpnames=apply(as.matrix(annotation[,c("gene","chr","start", "end")]),1,paste,collapse=":")[i];
    if(i %% 20 ==1) cat("Calculating Similarities for the ",i, "th region... \n");
    inc=(OmicsDataEachMap$chromosome==Annotation$chr[i]) & (OmicsDataEachMap$position <= Annotation$end[i]) & (OmicsDataEachMap$position >= Annotation$start[i]);
    if(sum(inc)>0)
    {
      input=OmicsDataEach[,inc];
      if(is.null(dim(input))) {input=matrix(input,ncol=1);rownames(input)=rownames(OmicsDataEach);}
      if(ncol(input)>1)
      {
        # impute the missings with the average values #
        nas=is.na(input);
        tmpimpute=apply(input,2,mean,na.rm=T);
        tmpimputegs=matrix(rep(tmpimpute,nrow(input)),nrow=nrow(input),ncol=ncol(input),byrow=T)
        input[nas]=tmpimputegs[nas];
        tmp3=getKernel(Input=input,kernels = kernels, Genomic = Genomic, AllRegions=0)
      }
      if(ncol(input)==1 | is.null(dim(input)))
      {
        tmpimpute=mean(input,na.rm=T);
        input[is.na(input)]=tmpimpute;
        input=matrix(input,ncol=1);rownames(input)=rownames(OmicsDataEach);
        tmp3=getKernel(Input=input,kernels = kernels, Genomic = Genomic, AllRegions=0)
      }
      Kinship[[index]]=tmp3;names(Kinship)[index]=tmpnames;
    }
    if(sum(inc)==0){Kinship[[index]]=list();names(Kinship)[index]=tmpnames;}
    index=index+1;
  }

  # Checking duplicates for similarities #
  IncludeRegionsNames=IncludeRegionsNamesAll=names(Kinship);
  if(length(Kinship)>1)
  {
    for(k in 2:length(Kinship))
    {
      keep=TRUE;current=Kinship[[k]]
      if(!is.na(current))
      {
        for(j in (k-1):1)
        {
          checking=Kinship[[j]]
          if(!is.na(checking)) {if(all.equal(checking,current)==TRUE) {keep=FALSE;break;}}
        }
        if(!keep)
        {
          IncludeRegionsNames[j]=paste(IncludeRegionsNames[j],IncludeRegionsNames[k],sep="=");
          IncludeRegionsNames[k]=NA;
          Kinship[[k]]=list();
        }
      }
    }
  }

  Result=list();
  Result$Kinship=Kinship;
  Result$IncludeRegions=IncludeRegionsNames;
  Result$IncludeRegionsAll=IncludeRegionsNamesAll;
  Result
}

getKernel<-function(Input,kernels,Genomic="FALSE")
{
  IDs=rownames(Input)
  output=NA;
  if(kernels=="linear" & Genomic) output=rrBLUP::A.mat(Input);
  if(kernels=="linear" & !Genomic)
  {
    nas=is.na(Input);
    # impute the missings with the average values #
    tmpimpute=apply(Input,2,mean,na.rm=T);
    tmpimputegs=matrix(rep(tmpimpute,nrow(Input)),nrow=nrow(Input),ncol=ncol(Input),byrow=T)
    Input[nas]=tmpimputegs[nas];
    tmp3=Input %*% t(Input)/ncol(Input); #tmp3=tmp/tmp2;
    output=tmp3
  }
  if(kernels=="IBS" & Genomic){
  output=varComp::IBS(Input)
  }
  if(!is.na(output)) colnames(output)=rownames(output)=IDs
  output;
}





