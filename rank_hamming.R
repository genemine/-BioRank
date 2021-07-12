library(parallel)
library(doParallel)
library(foreach)
library(fpc)
library(parallel)
library(doParallel)
library(foreach)
library(igraph)
library(e1071)
load("./msigdb_v6.Rdata")

#obtain Biology similarity
BioSimilarity<-function(data,h=1,k=4,size =100){
  if(is.null(rownames(data)))
    stop("data must have rownames(gene names)")
  if(length(which(duplicated(toupper(rownames(data)))))){
    message("gene names have duplicate values after process, these gene will be delete ...")
    a = toupper(rownames(data))
    b = a[duplicated(a)]
    c = c()
    for (i in 1:length(b)) {
      d = which(a == b[i])
      for (j in 1:length(d)) {
        c = c(c , d[j])
      }
    }
    data = data[-c,]
  }
  rownames(data) = toupper(rownames(data))
  
  #choose gene sets
  message("Choosing gene sets...")
  genesetClass=c("c2.all","c5.all","c6.all","c1.all","c3.all","c4.all","c7.all")
  allgeneset = c()
  for(j in 1:length(genesetClass)){
    name = genesetClass[j];
    geneset=msigdb_v6[[`name`]]
    allgeneset=c(allgeneset,geneset)
  }
  badset = c()
  for (i in 1:length(allgeneset)) {
    common_genes = intersect(toupper(allgeneset[[i]]),rownames(data))
    if(length(common_genes)/length(allgeneset[[i]]) < h){
      badset = c(badset,i)
    }
    allgeneset[[i]] = common_genes
  }
  allgeneset = allgeneset[-badset]
  
  if(length(allgeneset)<2000)
    message(paste("The number of annotation gene sets is small:",length(allgeneset)
                  ,", it is recommended to lower the threshold h",sep = ""))
  
  
  message("Clustering by each gene set...")
  #clustering by PAM
  cluster_res = list()
  cl <- makeCluster(k)
  registerDoParallel(cl)
  cluster_res = foreach(gs = 1:length(allgeneset),.combine = "c",.packages = c("fpc","e1071"),.export = "Study") %dopar%
    Study(gs = gs,allgeneset = allgeneset,data = data,size = size)
  stopCluster(cl)
  
  # for (gs in 1:length(allgeneset)) {
  #   cluster_res = c(cluster_res,Study(gs,allgeneset,data))
  # }
  
  message("Coustructing consensus matrix...")
  #coustruct consensus matrix
  nullres = 0
  nsamples = ncol(data)
  global_matrix = matrix(0,nsamples,nsamples)
  rownames(global_matrix) = colnames(data)
  colnames(global_matrix) = colnames(data)
  for (n in 1:length(cluster_res)) {
    if(!is.null(cluster_res[[n]])){
      similar_matrix = matrix(0,nsamples,nsamples)
      for (i in 1:(nsamples-1)) {
        for (j in (i+1):nsamples) {
          if(cluster_res[[n]][[i]] == cluster_res[[n]][[j]]){
            similar_matrix[i,j] = 1
            similar_matrix[j,i] = 1
          }
        }
      }
      diag(similar_matrix) = 1
      global_matrix = global_matrix + similar_matrix
    }else{
      nullres = nullres + 1
    }
  }
  sim = global_matrix/(length(cluster_res)-nullres)
  sim = list(sim)
  usedGeneSetNum = list(length(cluster_res))
  names(usedGeneSetNum) = "usedGeneSetNum"
  names(sim) = "SimilarityMatrix"
  return(c(sim,usedGeneSetNum))
}

#similarity study
Study<-function(gs,allgeneset,data,size){
  print(gs)
  common_genes = intersect(toupper(allgeneset[[gs]]),rownames(data))
  x = data[common_genes,]
  errCol = x[, sapply(x, function(v) var(v, na.rm=TRUE)==0)]
  errCol = as.matrix(errCol)
  if(ncol(errCol)==0){
    if((length(which(x == 0))/ncol(data))<length(common_genes)){
      if(nrow(x)<size){
        a = as.vector(as.data.frame(t(x)))
        codematrix = matrix(0,ncol = ncol(data),nrow = (nrow(x)*(nrow(x)-1))/2)
        # codematrix = as.list(as.data.frame(t(codematrix)))
        row = 1
        for (i in 1:(nrow(x)-1)) {
          for (j in (i+1):nrow(x)) {
            # codematrix[[row]][which(a[[i]]>a[[j]])]=1
            codematrix[row,which(a[[i]]>a[[j]])]=1
            row = row + 1
          }
        }
        # codematrix = as.data.frame(codematrix)
        # codematrix = as.matrix(codematrix)
        # dist = hamming.distance(codematrix)
        dist = hamming.distance(t(codematrix))
        res = pamk(as.dist(dist))
        return(list(res$pamobject$clustering))
      }else{
        res = pamk(as.dist(1-cor(x)))
        return(list(res$pamobject$clustering))
      }
    }else{
      return(list())
    }
  }else{
    return(list())
  }
}



#estimates the optimal number of clusters
# esOptimalNum<-function(sim,range = 2:20){
esOptimalNum<-function(sim,data,kmax = 15){
  K <- 2:kmax
  rst <- sapply(K, function(i){
    print(paste("K=",i))
    result <- pamk(t(data),i)
    stats <- cluster.stats(d = as.dist(1-sim), clustering = result$pamobject$clustering)
    stats$avg.silwidth
  })
  plot(K,rst,type='l',main='轮廓系数与K的关系', ylab='轮廓系数')
  k = c()
  for (i in 2:length(rst)) {
    if(rst[i]>rst[i-1]){
      k = c(k,i)
    }
  }
  k = intersect(order(rst,decreasing = TRUE),k)
  return(paste("K=",K[k[1]]))
}


res = BioSimilarity(data , h = 1)
sim = res$SimilarityMatrix 

