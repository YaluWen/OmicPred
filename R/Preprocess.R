#' @title Multi-kernel linear mixed model with adaptive lasso for prediction analysis on high-dimensional multi-omics data
#'
#' @description
#' \code{ReadOmicPLMM} prepares the data for analyses
#' @details
#' This function is used to pre-process the data used by the Multi-Kernel Linear Mixed Model with Adaptive Lasso method.
#' @param OmicsData A list with each element represent data from one omic. The rows is observations, and column is the variants
#' @param OmicsDataMap A list with each element represent data annotation from one omic. The rows is variants, and column is the annoation. It should have at least 2 columns ("chromosome","position"). The number of rows should be the same as the number of columns in the corresponding omic data
#' @param phenofile The file containing phenotype information. If provided, will ignore the affected values from the .fam file. The file should contain two columns: pheno: phenotypic values. ID: subject ID for each indidivuals.
#' @param trainID Subject ID for the training individuals.
#' @param annotation A data frame providing information to cut the genomes. It has at least four columns: gene, chr, start, end. The gene column must be unique. The other three columns contain information about chromosome, start and end of each regomic region. If Kinship is not provided, then annotation must be provided.
#' @param Y A vector of phenotypes with each name being subject ID. If provided, will ignore \code{phenofile}.
#' @param X A matrix of demographic variables, should be of the same order as Y (i.e. \code{rownames(X)=names(Y)}). The intercept column is not needed.
#' @param kernelsOmics kernels used for other omic data. Currently, it only takes linear kernel for non-genomic data. For genomic data, it can take linear and IBS.  Its length should be the same to the number of omics data if provided.
#' @param AllRegion A bool to indicate whether GSM should be calculated from the entire genome. Default=0.
#' @return A list that contains Y, Similarities for each omics and trainID, all of which are needed for OmicsPLMMPred algorithm(i.e. the function \code{OmicsPLMMPred}).
#' @examples
#' ## Using the data with pre-calculated GSMs ##
#' load(system.file("extdata", "Data.rda", package = "KLMMAL"))
#'
#' ## read both genotype and phenotype files from the disk, and calculate the GSMs using linear kernels ##
#' trainID=read.table(system.file("extdata", "TrainID.txt", package = "KLMMAL"));
#' annotation=read.table(system.file("extdata", "Annotation.txt", package = "KLMMAL"),header=T);
#' Data=ReadGenomicPLMM(bed=system.file("extdata", "testing.bed", package = "KLMMAL"),trainID=trainID[,1], annotation = annotation, phenofile=system.file("extdata", "pheno.txt", package = "KLMMAL"), kernels = c("linear"), AllRegions = 0);
#'
#' ## read genotype from the disk and calculate the GSMs using linear and poly2 kernels##
#' pheno=read.table(system.file("extdata", "pheno.txt", package = "KLMMAL"),header=T);
#' Y=pheno$pheno;names(Y)=pheno$ID;
#' trainID=read.table(system.file("extdata", "TrainID.txt", package = "KLMMAL"));
#' annotation=read.table(system.file("extdata", "Annotation.txt", package = "KLMMAL"),header=T);
#' Data=ReadGenomicPLMM(bed=system.file("extdata", "testing.bed", package = "KLMMAL"),trainID=trainID[,1], annotation = annotation, Y = Y, kernels = c("linear", "poly2"), AllRegions = 0);
#' @export
ReadOmicPLMM<-function(OmicsData=list(),  OmicsDataMap=list(),phenofile, trainID ,OmicsKernelMatrix=list(), annotation=NULL, Y=NULL, X=NULL, kernelsrOmics=NA, AllRegions=0)
{
  KernerlOutput=list(); kernelstart=1;

  ## import phenotypes ##
  if(is.null(Y) & missing (phenofile) ) stop("Error: You must provide phenotypes for the training data!")
  if(is.null(Y))
  {
    pheno=read.table(phenofile,header = T) # pheno, and ID #
    if(sum(colnames(pheno) %in% c("pheno", "ID"))!=2) stop("Error: Phenotype files have contain pheno and ID columns")
    Y=pheno$pheno;
    names(Y)=pheno$ID;
  }
  if(is.null(names(Y))) stop("Error: There is no IDs for phenotype data!")
  trainID=unique(trainID);
  if(sum(trainID %in% names(Y))!=length(trainID)) stop("Error: The data does not have all the training samples. Check!")
  train.index=which(names(Y) %in% trainID)
  if(sum(!is.na(Y[train.index])<5)) stop("Error: Less than 5 subjects for the training samples. The sample size is too small. ")

  ## Checking Annotations ##
  if(is.null(annotation)) stop("Error: You must provide the start and end of each genomic regions. If you have precalculated all Similarities, then you don't need to run this function")
  if(nrow(annotation)<1 & !AllRegions) stop("Error: Must have at least one genomic regions")
  if(sum(colnames(annotation) %in% c("gene","chr","start", "end"))!=4) stop("Error: the annotation must be a dataframe with at least four columns: gene, chr, start, end (gene/region name(unique), chromosome, start and end of the genomic regions)")
  if(length(unique(annotation$gene))!=nrow(annotation)) stop("Error: gene/region names must be unique")

  ## Checking Omics Data ##
  sn=rank(names(Y));
  if(length(OmicsData)>0)
  {
    for(i in 1:length(OmicsData))
    {
      if(nrow(OmicsData[[i]])!=length(Y)) stop(paste("Error: The", i,"th omics data does not have the same number of observations as the outcome, please check!"))
      if(sum(!(names(Y) %in% rownames(OmicsData[[i]])))!=0) stop(paste("Error: The individuals in the phenotype files and the individuals in the", i, "th omics data do not match"))
      if(sum(!(rownames(OmicsData[[i]]) %in% names(Y)))!=0) stop(paste("Error: The individuals in the phenotype files and the individuals in the", i, "th omics data do not match"))
      if(sum(names(Y)!=rownames(OmicsData[[i]]))!=0) {
        warning(paste("The individual orders in the phenotype files and the individuals in the", i, "th omics data do not match, rearrange the omics data!"))
        or=order(rownames(OmicsData[[i]]));
        OmicsData[[i]]=OmicsData[[i]][or,][sn,]
        if(sum(names(Y)!=rownames(OmicsData[[i]]))!=0) stop(paste("Error: check IDs between", i,"th omics and phenotypes!"))
      }
    }
  }
  if(length(OmicsKernelMatrix)>0)
  {
    for(i in 1:length(OmicsKernelMatrix))
    {
      if(length(OmicsKernelMatrix[[i]])!=nrow(annotation)) OmicsKernelMatrix[[i]]=list();
      if(length(OmicsKernelMatrix[[i]])>0)
      {
        if(length(Y)!=nrow(OmicsKernelMatrix[[i]])) stop(paste("Error: The", i,"th omics data matrix does not have the same number of observations as the outcome, please check!"))
        if(sum(!(names(Y) %in% colnames(OmicsKernelMatrix[[i]])))!=0) stop(paste("Error: The individuals in the phenotype files and the individuals in the", i, "th omics data matrix do not match"))
        if(sum(!(colnames(OmicsKernelMatrix[[i]]) %in% names(Y)))!=0) stop(paste("Error: The individuals in the phenotype files and the individuals in the", i, "th omics data matrix do not match"))
        if(sum(names(Y)!=colnames(OmicsKernelMatrix[[i]]))!=0) {
          warning(paste("The individual orders in the phenotype files and the individuals in the", i, "th omics data matrix do not match, rearrange the omics data!"))
          or=order(colnames(OmicsKernelMatrix[[i]]));
          for(i in 1:length(OmicsKernelMatrix)) OmicsKernelMatrix[[i]]=(OmicsKernelMatrix[[i]][or,or])[sn,sn]
          if(sum(names(Y)!=colnames(OmicsKernelMatrix[[i]]))!=0) stop(paste("Error: check IDs between", i,"th omics matrix and phenotypes!"))
        }
      }
    }
  }

  ## Deal with Omics Data ##
  {
    if(length(OmicsData)>0)
    {
      if(length(OmicsDataMap)!=length(OmicsData)) stop("You must provide the annotation (i.e., chr, position) for each omic data");
      if(!is.na(kernelsOtherOmics)) {
        kernelsOtherOmics=tolower(kernelsOtherOmics);
        if(length(kernelsOtherOmics)!=length(OmicsData)) warning("number of kernels does not match with the data, set to be all linear")
        if(sum(kernelsOtherOmics[toupper(names(OmicsData))!="SNP"]!="linear")!=0) warning("Currently only support linear kernels for non-genomic data")
        kernelsOtherOmics[toupper(names(OmicsData))!="SNP"]=rep("linear",length(OmicsData[toupper(names(OmicsData))!="SNP"]));
        if(length(kernelsOtherOmics[toupper(names(OmicsData))=="SNP"])>0)
        {
          if(!(tolower(kernelsOtherOmics[toupper(names(OmicsData))=="SNP"]) %in% c("linear","ibs")))
            kernelsOtherOmics[toupper(names(OmicsData))=="SNP"]='linear';
        }
      }
      if(is.na(kernelsOtherOmics)) kernelsOtherOmics=rep("linear",length(OmicsData));
      for(i in 1:length(OmicsData))
      {
        KernerlOutput[[kernelstart]]= list();
        OmicsDataEach=OmicsData[[i]];
        OmicsDataEachMap=OmicsDataMap[[i]];
        if(sum(colnames(OmicsDataEachMap) %in% c("position","chromosome"))!=2) stop("Error: the annotation for each omic must be a dataframe with at least four columns: chromosome, position)")
        if(length(OmicsDataEach)>0 & toupper(names(OmicsData))!="SNP") KernerlOutput[[kernelstart]]= CalculateSimilarity(OmicsDataEach,OmicsDataEachMap,Annotation=annotation,kernels=kernelsOtherOmics[i],AllRegions=AllRegions,Genomic=FALSE);
        if(length(OmicsDataEach)>0 & toupper(names(OmicsData))=="SNP") KernerlOutput[[kernelstart]]= CalculateSimilarity(OmicsDataEach,OmicsDataEachMap,Annotation=annotation,kernels=kernelsOtherOmics[i],AllRegions=AllRegions,Genomic=TRUE);
        kernelstart=kernelstart+1;
      }
    }

    if(length(OmicsKernelMatrix)>0)
    {
      if(length(OmicsKernelMatrix)!=length(OmicsData)) warning("Given matrix for omics data not match to the number of other omics, ignore this input")
      if(length(OmicsKernelMatrix)==length(OmicsData))
      {
        for(i in 1:length(OmicsKernelMatrix))
        {
          if(length(OmicsKernelMatrix[[i]])>0) {KernerlOutput[[i]]$Kinship=OmicsKernelMatrix[[i]]}
        }
      }
    }

  }

  ## Compare the consistencies between X and Ys ##
  if(!is.null(X))
  {
    if(rownames(X)!=names(Y)) stop("X and Y should be of the same order")
  }

  Data=list();
  Data$Y=Y;
  Data$trainID=trainID;
  Data$KernerlOutput=KernerlOutput;
  Data$X=X;
  Data
}

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






