getSelected <- function(data, strSelect) {
  if (!is.null(data)) {
    if (!is.na(strSelect)) {
      name = strSelect %>% strsplit(split = "sep") %>% do.call(c, .)
      data = data[ ,..name]
    } else {
      name = names(data)
    }
    size = ncol(data)
    return (list(data = data, name = name, size = size))
  } else {
    return (list(data = NULL, name = character(0), size = 0))
  }
}

########################### iPat Genomic Selection  #############################
getCovFromGWAS <- function(isGWASAssit, cutoff,
                           sizeN, dataCov,
                           nameProject, nameTrait,
                           rawGenotype, rawMap) {
  ## Config
  cov = dataCov$data
  sizeCov = dataCov$size
  if (isGWASAssist) {
    cat("   Loading QTNs information ...")
    ## Read GWAS result
    tableGWAS = fread(sprintf("iPat_%s_%s_GWAS.txt", nameProject, nameTrait))
    ## Merge GM and p-value
    names(rawMap)[1] = "SNP"
    mapGWAS = data.table(rawMap)
    mapGWAS[ ,P.value := tableGWAS$P.value[match(rawMap$SNP, tableGWAS$SNP)]]
    mapGWAS$P.value[is.na(mapGWAS$P.value)] = 1
    ## Order p-value
    snpOrder = order(mapGWAS$P.value)
    mapGWAS = mapGWAS[snpOrder]
    genotype = rawGenotype[ ,..snpOrder]
    ## Find QTNs
    indexSig = which(mapGWAS$P.value < (cutoff/nrow(tableGWAS)))
    ## Generate a dataframe by number of QTNs
    ### 0 QTNs
    if (length(indexSig) == 0) {
      cGWAS = NULL
      sizeQTN = 0
      ### 1 QTNs
    } else if (length(indexSig) == 1) {
      cGWAS = data.frame(m = genotype[ ,..indexSig])
      sizeQTN = 1
      ### 1+ QTNs
    } else {
      cGWAS = genotype[ ,..indexSig] %>% as.data.frame()
      ## LD Remove
      LD_remain = Blink.LDRemove(cGWAS, .7, indexSig, orientation = "col")
      cGWAS = cGWAS[ ,LD_remain]
      sizeQTN = length(LD_remain)
    }
    cat("Done\n")
  } else {
    cGWAS = NULL
    sizeQTN = 0
  }
  # Prevent c > n
  ## if c + qtn > n
  if (sizeN < sizeCov + sizeQTN) {
    diff = sizeCov + sizeQTN - sizeN
    return (data.frame(cov, cGWAS[ ,1 : (sizeQTN - diff)]))
    ## both cov nor qtn has size greater than 0
  } else if (sizeCov != 0 && sizeQTN != 0) {
    return (data.frame(cov, cGWAS))
    ## if c + qtn <= n
  } else if (sizeCov == 0 && sizeQTN != 0) {
    return (cGWAS)
    ## if only c
  } else if (sizeCov != 0 && sizeQTN == 0) {
    return (cov)
  } else if (sizeCov == 0 && sizeQTN == 0) {
    return (NULL)
  }
}

runCrossValidation <- function(finalP, finalG, finalC, isRaw, countFold, countIter, project, trait) {
  # write header
  nameTableRaw = sprintf("iPat_%s_%s_raw", project, trait)
  nameTableAcc = sprintf("iPat_%s_%s_acc", project, trait)
  if (isRaw)
    cat("obs\tpre\ttrait\titer\tfold\n", file = paste0(nameTableRaw, ".txt"))
  cat("r\trmse\ttrait\titer\tfold\n", file = paste0(nameTableAcc, ".txt"))
  # get sample size
  n = length(finalP)
  # For GBLUP
  if (!is.null(rawKin)) {
    finalK = rawKin
  } else {
    finalK = A.mat(as.matrix(finalG))
  }
  # loop over iterations
  for (iter in 1:countIter) {
    idxFold = rep(1 : countFold, length = n) %>% sample()
    # loop over fold
    for (fold in 1:countFold) {
      # Progress Message
      cat(sprintf("iter : %d/%d, fold : %d/%d\n",
                  iter, countIter, fold, countFold))
      ## Get index for testing set
      idxTest = idxFold == fold
      ## Get phenotype for testing and one with missing value
      yTemp = finalP
      yTemp[idxTest] = NA
      ## Predict
      switch (ANALYSIS,
              "gBLUP" = {
                prediction = getAccByGBLUP(finalG, finalP, finalC, finalK, yTemp, idxTest)
              },
              "rrBLUP" = {
                prediction = getAccByRRBLUP(finalG, finalP, finalC, yTemp, idxTest)
              },
              "BGLR" = {
                prediction = getAccByBGLR(finalG, finalP, finalC, yTemp, idxTest)
              }
      )
      ## Export
      if (isRaw)
        writeRaw(finalP[idxTest], prediction$pre, trait, iter, fold, paste0(nameTableRaw, ".txt"))
      writeAcc(prediction$r, prediction$rmse, trait, iter, fold, paste0(nameTableAcc, ".txt"))
    }
  }
  # Plot
  plotCV(isRaw, nameTableRaw, nameTableAcc)
}

rmse <- function(x, y) {
  sqErr = sum((x - y)^2)
  rmse = sqrt(sqErr/length(y))
  return (rmse)
}

## ---------------------------- gBLUP ---------------------------- ##
runGBLUP <- function(finalP, finalG, finalC, taxa, project, trait) {
  # Write header
  nameTableMarker = sprintf("iPat_%s_%s_GEBV_By_Marker.txt", project, trait)
  nameTableGEBV = sprintf("iPat_%s_%s_GEBV", project, trait)
  nameTableStat = sprintf("iPat_%s_%s_Stat.txt", project, trait)
  cat("Stat\tValue\n", file = nameTableStat)
  # Compute kinship
  if (!is.null(rawKin)) {
    finalK = rawKin
  } else {
    finalK = A.mat(as.matrix(finalG))
  }
  # Run rrBLUP for gBLUP
  gblup = mixed.solve(finalP, X = finalC, K = finalK)
  # Calculate random effect
  Zu = as.matrix(gblup$u)
  # Calculate fixed effect
  etd.beta = as.matrix(gblup$beta)
  # Calculate BLUE and BLUP
  if (!is.null(finalC)) {
    Xb = as.matrix(finalC) %*% etd.beta
  } else {
    Xb = matrix(etd.beta, ncol = 1, nrow = nrow(finalG))
  }
  GEBV = Zu + Xb
  # Export marker effects
  write.table(x = data.frame(u = Zu), file = nameTableMarker, quote = F,
              row.names = T, col.names = T, sep = "\t")
  # Export GEBVs
  write.table(x = data.frame(taxa, GEBV), file = paste0(nameTableGEBV, ".txt"), quote = F,
              row.names = F, col.names = T, sep = "\t")
  # Export Stats
  beta.name = names(gblup$beta)
  Stat = c("Vu", "Ve", paste0("beta.", beta.name), "LL")
  write.table(data.frame(Stat,
                         Value = c(gblup$Vu, gblup$Ve, gblup$beta, gblup$LL)),
              file = nameTableStat, append = TRUE, row.names = F, col.names = F, sep = "\t")
  # Export plot
  plotHistBV(nameTableGEBV)
}

getAccByGBLUP <- function(X, Y, C, K, YTemp, idxTest) {
  gblup = mixed.solve(YTemp, X = C, K = K)
  Zu = as.matrix(gblup$u[idxTest])
  etd.beta = as.matrix(gblup$beta)
  if (!is.null(C))
    Xb = as.matrix(C)[idxTest, ] %*% etd.beta
  else
    Xb = matrix(etd.beta, ncol = 1, nrow = sum(idxTest))
  gebv = Xb + Zu
  acc_r = cor(gebv, Y[idxTest]) %>% c()
  acc_rmse = rmse(gebv, Y[idxTest])
  return (list(r = acc_r, rmse = acc_rmse, pre = gebv))
}

## ---------------------------- rrBLUP ---------------------------- ##
runRRBLUP <- function(finalP, finalG, finalC, taxa, project, trait) {
  # Write header
  nameTableMarker = sprintf("iPat_%s_%s_marker.txt", project, trait)
  nameTableGEBV = sprintf("iPat_%s_%s_GEBV", project, trait)
  nameTableStat = sprintf("iPat_%s_%s_Stat.txt", project, trait)
  cat("Stat\tValue\n", file = nameTableStat)
  # Run rrBLUP
  rrblup = mixed.solve(finalP, X = finalC, Z = finalG)
  # Calculate random effect
  etd.u = as.matrix(rrblup$u)
  # Calculate fixed effect
  etd.beta = as.matrix(rrblup$beta)
  # Calculate BLUE and BLUP
  Zu = as.matrix(finalG) %*% etd.u
  if (!is.null(finalC)) {
    Xb = as.matrix(finalC) %*% etd.beta
  } else {
    Xb = matrix(etd.beta, ncol = 1, nrow = nrow(finalG))
  }
  GEBV = Zu + Xb
  # Export marker effects
  write.table(x = data.frame(u = etd.u), file = nameTableMarker, quote = F,
              row.names = T, col.names = T, sep = "\t")
  # Export GEBVs
  write.table(x = data.frame(taxa, GEBV), file = paste0(nameTableGEBV, ".txt"), quote = F,
              row.names = F, col.names = T, sep = "\t")
  # Export Stats
  beta.name = names(rrblup$beta)
  Stat = c("Vu", "Ve", paste0("beta.", beta.name), "LL")
  write.table(data.frame(Stat,
                         Value = c(rrblup$Vu, rrblup$Ve, rrblup$beta, rrblup$LL)),
              file = nameTableStat, append = TRUE, row.names = F, col.names = F, sep = "\t")
  # Export plot
  plotHistBV(nameTableGEBV)
}

getAccByRRBLUP <- function(X, Y, C, YTemp, idxTest) {
  rrblup = mixed.solve(YTemp, X = C, Z = X)
  etd.u = as.matrix(rrblup$u)
  etd.beta = as.matrix(rrblup$beta)
  Zu = as.matrix(X[idxTest, ]) %*% etd.u
  if (!is.null(C))
    Xb = as.matrix(C)[idxTest, ] %*% etd.beta
  else
    Xb = matrix(etd.beta, ncol = 1, nrow = sum(idxTest))
  gebv = Xb + Zu
  acc_r = cor(gebv, Y[idxTest]) %>% c()
  acc_rmse = rmse(gebv, Y[idxTest])
  return (list(r = acc_r, rmse = acc_rmse, pre = gebv))
}

## ---------------------------- BGLR ---------------------------- ##
runBGLR <- function(finalP, finalG, finalC, taxa, project, trait) {
  # Write header
  nameTableMarker = sprintf("iPat_%s_%s_marker.txt", project, trait)
  nameTableGEBV = sprintf("iPat_%s_%s_GEBV", project, trait)
  nameBGLR = sprintf("iPat_%s_", project)
  # Run BGLR
  ETA = list(list(X = finalG, model = modelBGLR))
  if (!is.null(finalC)) {
    ETA[[2]] = list(X = finalC, model = "FIXED")
  }
  bglr = BGLR(y = finalP, ETA = ETA, verbose = FALSE, saveAt = nameBGLR,
              nIter = nIter, burnIn = burnIn)
  # GEBVs
  GEBV = bglr$yHat
  # Export marker effects
  write.table(x = data.frame(u = bglr$ETA[[1]]$b), file = nameTableMarker, quote = F,
              row.names = T, col.names = T, sep = "\t")
  # Export GEBVs
  write.table(x = data.frame(taxa, GEBV = bglr$yHat), file = paste0(nameTableGEBV, ".txt"), quote = F,
              row.names = F, col.names = T, sep = "\t")
  # Export plot
  plotHistBV(nameTableGEBV)
}

getAccByBGLR <- function(X, Y, C, YTemp, idxTest) {
  nameBGLR = sprintf("iPat_%s_", project)
  ETA = list(list(X = X, model = modelBGLR))
  if (!is.null(C)) {
    ETA[[2]] = list(X = C, model = "FIXED")
  }
  bglr = BGLR(y = YTemp, ETA = ETA, verbose = FALSE, saveAt = nameBGLR,
              nIter = nIter, burnIn = burnIn)
  gebv = bglr$yHat[idxTest]
  acc_r = cor(gebv, Y[idxTest]) %>% c()
  acc_rmse = rmse(gebv, Y[idxTest])
  return (list(r = acc_r, rmse = acc_rmse, pre = gebv))
}


## ---------------------------- Export docs ---------------------------- ##
writeRaw <- function(obs, pre, trait, iter, fold, filename) {
  fwrite(data.table(round(obs, 4), round(pre, 4), trait, iter, fold),
         file = filename,
         row.names = F, col.names = F, sep = "\t", quote = F, append = TRUE)
}

writeAcc <- function(r, rmse, trait, iter, fold, filename) {
  cat(sprintf("%.4f\t%.4f\t%s\t%d\t%d\n",
              r, rmse, trait, iter, fold),
      file = filename, append = TRUE)
}


############################# iPat Specific Ends ###############################

#        LD_remain = Blink.LDRemove(cGWAS, .7, indexSig, orientation = "col")
`Blink.LDRemove`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
  #`Blink.LDRemovebackup`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
  #Objects: LD remove, especially length(Porder)>10000
  #Authors: Yao Zhou
  #Last update: 08/15/2016
  seqQTN=NULL
  is.done=FALSE
  l=1000
  lp=length(Porder)
  tt=1
  while(!is.done){
    tt = tt+1
    Pordern=Blink.LDRemoveDivided(GDneo=GDneo,LD=LD,Porder=Porder,orientation=orientation,model=model,l=lp)
    index=Porder %in% Pordern
    if(orientation=="col"){
      GDneo = GDneo[,index]
    }else{
      GDneo = GDneo[index,]
    }
    ls=length(Pordern)
    if(ls==lp) lp=l*tt
    if(ls<=lp){
      is.done=TRUE
    }
    Porder = Pordern
  }
  if(length(Porder) > 1){
    seqQTN=Blink.LDRemoveBlock(GDneo=GDneo,LD=LD,Porder=Porder,orientation=orientation,model=model)
  }else{
    seqQTN = Porder
  }
  return(seqQTN)
}
