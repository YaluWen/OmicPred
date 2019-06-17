#' @title Multi-kernel linear mixed model with adaptive lasso for prediction analysis on high-dimensional multi-omics data
#'
#' @description
#' \code{OmicsPLMM} prepares the data for analyses
#' @details
#' This function is used to pre-process the data used by the Multi-Kernel Linear Mixed Model with Adaptive Lasso method.
#' @param OmicsData A list with each element representing data from one omic. Each row is an observations, and each column is the variants. NA represents missing. For the subject with missing omcis, then the subject should be in the corresponding omic, but with missing entry.
#' @param OmicsDataMap A list with each element represent data annotation from one omic. The row is variants, and column is the annoation. It should have at least 2 columns ("chromosome","position"). The number of rows should be the same as the number of columns in the corresponding omic data
#' @param phenofile The file containing phenotype information. If provided, will ignore the affected values from the .fam file. The file should contain two columns: pheno: phenotypic values. ID: subject ID for each indidivuals.
#' @param trainID Subject ID for the training individuals.
#' @param OmicsKernelMatrix XXX
#' @param annotation A data frame providing information to cut the genomes. It has at least four columns: gene, chr, start, end. The gene column must be unique. The other three columns contain information about chromosome, start and end of each regomic region. #####If Kinship is not provided, then annotation must be provided.
#' @param Y A vector of phenotypes with each name being subject ID. If provided, will ignore \code{phenofile}.
#' @param X A matrix of demographic variables, should be of the same order as Y (i.e. \code{rownames(X)=names(Y)}). The intercept column is not needed.
#' @param kernelsOmics kernels used for other omic data. Currently, it only takes linear kernel for non-genomic data. For genomic data, it can take linear and IBS.  Its length should be the same to the number of omics data if provided.
#' @param AllRegion A bool to indicate whether GSM should be calculated from the entire genome. Default=0.
#' @param stdData A bool to indicate whether data should be centered and scaled with sd=1. Default=TRUE.
#' @param predict A bool indicating if prediction is needed. Default is FALSE.
#' @param weight.fixed Weight that one wants to assign to covariate. Default is NULL.
#' @param weight.random Weight that one wants to assign to random effect. This should be as the same number of random effects included in the model. Suggest to give NULL, and default is NULL.
#' @param maxiter Maximum number of iterations. Default is 1000.
#' @param outputall Controls the details of output. \code{outputall=0}, only the prediction and original outputs are saved.  \code{outputall=1}, in addition to the above, the model fit is also saved. \code{outputall=2}, in addition to the above, the similarity matrix is also save (can be large, not suggested).
#' @return A list that contains original output and the predicted value (if any). The other output is controlled by \code{outputall}.
#' @export
OmicsPLMM<-function(OmicsData=list(),  OmicsDataMap=list(),
                   phenofile, trainID ,OmicsKernelMatrix=list(),
                   annotation=NULL, Y=NULL, X=NULL, kernelsrOmics=NA, AllRegions=0,
                   stdData = TRUE,
                   predict = FALSE, weight.fixed = NULL, weight.random = NULL,
                   maxiter = 1000, outputall = FALSE)
{

  if(stdData & length(OmicsData)>0)
  {
    for(i in 1:length(OmicsData))
    {
      if(length(OmicsData[[i]])>0)
      {
        tmp=OmicsData[[i]];
        if(!is.null(dim(tmp))){
          if(nocl(tmp)!=1) {tmp=apply(tmp,2,scale,na.rm = TRUE);rownames(tmp)=rownames(OmicsData[[i]])}
          if(ncol(tmp)==1) tmp=tmp[,1];
        }
        if(is.null(dim(tmp))){tmp=scale(tmp);tmp=matrix(tmp,ncol=1);rownames(tmp)=rownames(OmicsData[[i]])}
        OmicsData[[i]]=tmp;
      }
    }
  }
  Data=ReadOmicPLMM<-function(OmicsData=OmicsData,  OmicsDataMap=OmicsDataMap,phenofile=phenofile, trainID=trainID ,OmicsKernelMatrix=OmicsKernelMatrix, annotation=annotation, Y=Y, X=X, kernelsrOmics=kernelsrOmics, AllRegions=AllRegions)

  Result=OmicsPLMMPred(Data=Data,predict=predict, weight.fixed=weight.fixed,weight.random=weight.random,maxiter=maxiter,outputall=outputall,minheri=0.01,lambdarange=c(0,100), tol=1e-6,crit='bic')
  Result
}

