

create.visualization.artifacts <- function(covs) {
    
    ## Check which columns of covs are discrete or continuous
    covsClass <- sapply(as_tibble(covs), typeof)
    covsClass_num <- names(covsClass[which(covsClass %in% c("integer","double"))])
    IntCovsaCols <- apply(covs[,covsClass_num],2,function(x) any(round(as.numeric(x)) != x))
    MaybeCont <- names(IntCovsaCols[which(IntCovsaCols == FALSE)])
    MaybeContData <- apply(covs[,MaybeCont],2,function(x) any(length(levels(as.factor(x)))<(nrow(covs)*0.60)))
    ContinuousData <- names(MaybeContData[which(MaybeContData == FALSE)])
    ContinuousData <- c(ContinuousData,names(IntCovsaCols[which(IntCovsaCols == TRUE)]))
    
    covs.continuous <- covs[,which(colnames(covs) %in% ContinuousData)]
    covs.discrete <- covs[,which(!colnames(covs) %in% ContinuousData)]
    
    ret <- list(covs=covs, discreteCovs=colnames(covs.discrete), continuousCovs=colnames(covs.continuous))
    
  }

write.h5.artifact <- function(seurat, fn, covs, artifactName) {
  require(rhdf5)
  artifacts <- create.visualization.artifacts(covs)
  h5createGroup(fn,paste0("artifacts/",artifactName))
  i <- sapply(artifacts$covs, is.factor)
  artifacts$covs[i] <- lapply(artifacts$covs[i], as.character)
  
  h5write(artifacts$covs, fn, paste0("artifacts/",artifactName, "/covs"))
  h5write(artifacts$discreteCovs, fn, paste0("artifacts/",artifactName, "/discreteCovs"))
  h5write(artifacts$continuousCovs, fn,paste0("artifacts/",artifactName, "/continuousCovs"))
}


write.h5 <- function(seurat.list, fn) {
  unlink(fn)
  seurat = seurat.list$all$seurat  
  require(rhdf5)
  h5createFile(fn)
  
  h5createGroup(fn,"matrix")
  h5write(colnames(seurat), fn, "matrix/barcodes" )
  h5write(seurat@assays$RNA@counts@i, fn, "matrix/indices")
  h5write(seurat@assays$RNA@counts@p, fn, "matrix/indptr")
  h5write(seurat@assays$RNA@counts@x, fn, "matrix/data")
  h5write(seurat@assays$RNA@counts@Dim, fn, "matrix/shape")
  h5write(seurat@assays$RNA@counts@Dimnames[[1]], fn, "matrix/gene_names")
  
  h5createGroup(fn,"artifacts")
  for(group in names(seurat.list)) {
    write.h5.artifact(seurat.list[[group]]$seurat, fn, seurat.list[[group]]$covs, group)
  }
}

