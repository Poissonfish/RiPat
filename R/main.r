runArgs <- function(args) {
  # get package and method
  assign("pkgCalled", args[1], envir = .GlobalEnv)
  if ((pkgCalled == "GAPIT") ||
      (pkgCalled == "FarmCPU") || (pkgCalled == "PLINK"))
    isGWAS = TRUE
  else
    isGWAS = FALSE

  # anything other than PLINK
  if (pkgCalled != "PLINK") {
    for (i in 1:length(args)) {
      switch (
        args[i],
        "-wd" = {
          i = i + 1
          wd = args[i]
          setwd(wd)
        },
        "-project" = {
          i = i + 1
          project = args[i]
        },
        "-phenotype" = {
          cat("   Loading phenotype ...")
          i = i + 1
          rawPhenotype = fread(args[i])
          taxa = rawPhenotype[, 1]
          sizeN = nrow(rawPhenotype)
          cat("Done\n")
        },
        "-pSelect" = {
          i = i + 1
          selectP = args[i]
        },
        "-genotype" = {
          cat("   Loading genotype ...")
          i = i + 1
          rawGenotype = fread(args[i])
          # If have taxa column
          if (is.character(rawGenotype[[1]]))
            rawGenotype = rawGenotype[, -1]
          cat("Done\n")
        },
        "-map" = {
          cat("   Loading map ...")
          i = i + 1
          if (grepl("NA", args[i])) {
            rawMap = NULL
          } else {
            rawMap = fread(args[i])
            names(rawMap) = c("SNP", "Chromosome", "Position")
          }
          cat("Done\n")
        },
        "-cov" = {
          cat("   Checking covariates ...")
          i = i + 1
          if (grepl("NA", args[i])) {
            rawCov = NULL
          } else {
            rawCov = fread(args[i])
          }
          # If have taxa column
          if (is.character(rawCov[[1]]))
            rawCov = rawCov[, -1]
          cat("Done\n")
        },
        "-cSelect" = {
          i = i + 1
          selectC = ifelse(grepl("NA", args[i]),
                           NA, args[i])
        },
        "-kin" = {
          cat("   Checking kinship ...")
          i = i + 1
          if (grepl("NA", args[i])) {
            rawKin = NULL
          } else {
            rawKin = fread(args[i]) %>% as.data.frame()
          }
          # If have taxa column
          if (is.character(rawKin[[1]]))
            rawKin = rawKin[, -1]
          cat("Done\n")
        },
        "-gwas" = {
          i = i + 1
          isGWASAssist = as.logical(args[i])
          cutoff = 0.05
        },
        "-gs" = {
          i = i + 1
          isValid = as.logical(args[i])
          isRaw = TRUE
          i = i + 1
          countFold = as.numeric(args[i])
          i = i + 1
          countIter = as.numeric(args[i])
        },
        "-arg" = {
          if (pkgCalled == "GAPIT") {
            i = i + 1
            assign("model", args[i], envir = .GlobalEnv)
            i = i + 1
            assign("nPC", as.numeric(args[i]), envir = .GlobalEnv)
          } else if (pkgCalled == "FarmCPU") {
            i = i + 1
            assign("method.bin", args[i], envir = .GlobalEnv)
            i = i + 1
            assign("maxLoop", as.numeric(args[i]), envir = .GlobalEnv)            
          } else if (pkgCalled == "gBLUP") {
            # NONE
          } else if (pkgCalled == "rrBLUP") {
            # NONE
          } else if (pkgCalled == "BGLR") {
            i = i + 1
            assign("modelBGLR", args[i], envir = .GlobalEnv)
            i = i + 1
            assign("nIter", as.numeric(args[i]), envir = .GlobalEnv)
            i = i + 1
            assign("burnIn", as.numeric(args[i]), envir = .GlobalEnv)
            assign("rawKin", NULL, envir = .GlobalEnv)
          }
        }
      )
    }

    # Subset Phenotype 
    cat("   Subsetting phenotype ...")
    dataP = getSelected(rawPhenotype, selectP)
    cat("Done\n")

    # Subset Covariates
    cat("   Subsetting covariates ...")
    dataC = getSelected(rawCov, selectC)
    if (is.null(dataC$data)) {
      finalC = NULL
    } else {
      finalC = data.frame(taxa, dataC$data)
    }
    cat("Done\n")

    # iteration
    for (trait in dataP$name) {
      # GWAS (GAPIT or FarmCPU)
      if (isGWAS) {
        if (pkgCalled == "GAPIT") {
          x = GAPIT(
            Y = data.frame(taxa, dataP$data[[trait]]),
            GM = data.frame(rawMap),
            GD = data.frame(taxa, rawGenotype),
            KI = rawKin,
            CV = finalC,
            PCA.total = nPC,
            file.output = F,
            model = model
          )
          dt_gwas = x$GWAS
          dt_gwas = merge(rawMap, dt_gwas[, -c(2, 3)], sort = F) %>% data.table(effect = x$mc)
        } else if (pkgCalled == "FarmCPU") {
          x = FarmCPU(
            Y = data.frame(taxa, dataP$data[[trait]]),
            GM = data.frame(rawMap),
            GD = data.frame(taxa, rawGenotype),
            CV = finalC,
            method.bin = method.bin,
            maxLoop = maxLoop,
            MAF.calculate = TRUE,
            file.output = F
          )
          dt_gwas = x$GWAS
        }

        fwrite(
          x = dt_gwas,
          file = sprintf("iPat_%s_%s_GWAS.txt", project, trait),
          quote = F,
          row.names = F,
          sep = "\t"
        )
        dt_out = dt_gwas[, c(2, 1, 3, 4)] %>% data.frame()
        names(dt_out) = c("CHR", "SNP", "BP", "P")
        iPat.Manhattan(GI.MP = dt_out[, -2],
                       filename = sprintf("iPat_%s_%s", project, trait))
        iPat.QQ(dt_out$P, filename = sprintf("iPat_%s_%s", project, trait))
        iPat.Genotype.View(
          myGD = data.frame(taxa, rawGenotype),
          filename = sprintf("iPat_%s_%s", project, trait)
        )
        iPat.Phenotype.View(
          myY = data.frame(taxa, dataP$data[[trait]]),
          filename = sprintf("iPat_%s_%s", project, trait)
        )

      } else {
        # GS (gBLUP / rrBLUP / BGLR)
        cat(sprintf("   %s is computing for trait %s ...", pkgCalled, trait))
        # Collect covariates
        Cov = getCovFromGWAS(
          isGWASAssist,
          cutoff,
          sizeN = sizeN,
          dataCov = dataC,
          nameProject = project,
          nameTrait = trait,
          rawGenotype,
          rawMap
        )

        # Validation
        if (isValid) {
          # In case have any missing data
          idxNonNA = !is.na(dataP$data[[trait]])
          finalP = dataP$data[[trait]][idxNonNA]
          finalG = rawGenotype[idxNonNA,]
          finalC = Cov[idxNonNA,]
          # Run Validation
          runCrossValidation(finalP,
                             finalG,
                             finalC,
                             isRaw,
                             countFold,
                             countIter,
                             project,
                             trait)

          # No validation
        } else {
          finalP = dataP$data[[trait]]
          finalG = rawGenotype
          finalC = Cov
          if (pkgCalled == "gBLUP") {
            runGBLUP(finalP, finalG, finalC, taxa, project, trait)
          } else if (pkgCalled == "rrBLUP") {
            runRRBLUP(finalP, finalG, finalC, taxa, project, trait)
          } else if (pkgCalled == "BGLR") {
            runBGLR(finalP, finalG, finalC, taxa, project, trait)
          }
        }

        iPat.Genotype.View(
          myGD = data.frame(taxa, rawGenotype),
          filename = sprintf("iPat_%s_%s", project, trait)
        )
        iPat.Phenotype.View(
          myY = data.frame(taxa, dataP$data[[trait]]),
          filename = sprintf("iPat_%s_%s", project, trait)
        )
        cat("Done\n")
      }
    }
  } else {
    # ============= PLINK =============
    for (i in 1:length(args)) {
      print(i)
      switch (
        args[i],
        "-wd" = {
          i = i + 1
          wd = args[i]
          setwd(wd)
        },
        "-project" = {
          i = i + 1
          project = args[i]
        },
        "-phenotype" = {
          cat("   Loading phenotype ...")
          i = i + 1
          if (grepl("NA", args[i]))
            Y.path = "NA"
          else
            Y.path = args[i]
          cat("Done\n")
        },
        "-pSelect" = {
          i = i + 1
          selectP = args[i]
        },
        "-maf" = {
          i = i + 1
          maf = as.numeric(args[i])
        },
        "-ms" = {
          i = i + 1
          ms = as.numeric(args[i])
        },
        "-cSelect" = {
          i = i + 1
          selectC = ifelse(grepl("NA", args[i]),
                           NA, args[i])
        },
        "-genotype" = {
          cat("   Loading genotype ...")
          i = i + 1
          GD.path = args[i]
          cat("Done\n")
        },
        "-map" = {
          cat("   Loading map ...")
          i = i + 1
          GM.path = args[i]
          cat("Done\n")
        },
        "-cov" = {
          cat("   Checking covariates ...")
          i = i + 1
          if (grepl("NA", args[i]))
            C.path = "NA"
          else
            C.path = args[i]
          cat("Done\n")
        },
        "-arg" = {
          i = i + 1
          ci = as.numeric(args[i])
          i = i + 1
          model = args[i]
          i = i + 1
          pathPLINK = args[i]
        }
      )
    }

    # Subset Phenotype
    cat("   Loading phenotype ...")
    # If no phenotype provided, which would be in
    if (Y.path == "NA") {
      Y.data = fread(GD.path, na.strings = c("NA", "NaN")) %>% as.data.frame()
      Y.path = paste0(GD.path %>% substr(1, nchar(.) - 3), "_trait.txt")
      write.table(
        x = data.frame(
          FID = Y.data[, 1],
          IID = Y.data[, 2],
          trait = Y.data[, 6]
        ),
        file = Y.path,
        quote = F,
        row.names = F,
        sep = '\t'
      )
      trait.name = "trait"
      trait_count = 1
      suffix = ".assoc"
      # If phenotype provided
    } else {
      Y.data = fread(Y.path, na.strings = c("NA", "NaN"))
      G.data = fread(GD.path, na.strings = c("NA", "NaN")) %>% as.data.frame()
      # wrong format for PLINK
      FID = G.data[, 1]
      IID = G.data[, 2]
      taxa = IID
      # get selected data
      trait.name = selectP %>% strsplit(split = "sep") %>% do.call(c, .)
      Y.data = data.frame(FID = FID, IID = IID, subset(Y.data, ,trait.name))
      names(Y.data) = c("FID", "IID", trait.name)
      Y.path = paste0(Y.path %>% substr(1, nchar(.) - 4), "_trait.txt")
      trait_count = (names(Y.data) %>% length()) - 2
      suffix = ".qassoc"
    }
    cat("Done\n")

    # Covariate
    if (C.path != "NA") {
      cat("   Loading covariates ...")
      C.data = fread(C.path) %>% as.data.frame()
      if (is.character(C.data[, 1]))
        C.data = C.data[, -1]
      C.name = selectC %>% strsplit(split = "sep") %>% do.call(c, .)
      cat("Done\n")
    } else {
      C.name = character()
    }

    # Basic
    if (model == "GLM") {
      method = "--assoc"
    } else {
      method = "--logistic"
      suffix = paste0(suffix, ".logistic")
    }
    basic = sprintf(
      "%s --ped %s --map %s %s --allow-no-sex --adjust -ci %s --pheno %s --all-pheno --prune -out %s",
      paste0('"', pathPLINK, '"'),
      paste0('"', GD.path, '"'),
      paste0('"', GM.path, '"'),
      method,
      ci,
      paste0('"', Y.path, '"'),
      paste0('"', file.path(wd, project), '"')
    )
    # QC
    if (ms != 1)
      MS = sprintf("--geno %s", ms)
    else
      MS = character()
    if (maf != 0)
      MAF = sprintf("--maf %s", maf)
    else
      MAF = character()

    #Plotting
    setwd(wd)
    for (t in 1:trait_count) {
      ## Rewrite Phenotype
      if (Y.path != "NA") {
        tCol = t + 2
        updateY = Y.data[, c(1:2, tCol)]
        updateY[is.na(updateY[, 3]), 3] = -9
        fwrite(
          x = updateY,
          file = Y.path,
          quote = F,
          row.names = F,
          sep = '\t'
        )
      }

      ## COV and running PLINK
      if (length(C.name) > 0) {
        cov = sprintf(
          "--covar %s --covar-name %s",
          paste0('"', C.path, '"'),
          paste(C.name, collapse = ", ")
        )
        paste(basic, MS, MAF, cov) %>% system(input = "notepad")
      } else{
        paste(basic, MS, MAF) %>% system(input = "notepad")
      }

      #Loading data
      cat(sprintf("   Plotting trait %s ...", t))
      dt_gwas = fread(paste0(project, ".", trait.name[t], suffix), header = T)
      names(dt_gwas) = c("Chromosome",
                          "SNP",
                          "Position",
                          "NMISS",
                          "effect",
                          "SE",
                          "R2",
                          "T",
                          "P.value")
      dt_gwas[, c(2, 1, 3, 9, 5)] %>% fwrite(
        file = sprintf("iPat_%s_%s_GWAS.txt", project, trait.name[t]),
        quote = F,
        row.names = F,
        sep = "\t"
      )
      dt_out = dt_gwas[, c(1:3, 9)] %>% data.frame()
      names(dt_out) = c("CHR", "SNP", "BP", "P")
      iPat.Manhattan(
        GI.MP = dt_out[, -2],
        filename = sprintf("iPat_%s_%s", project, trait.name[t])
      )
      iPat.QQ(dt_out$P,
              filename = sprintf("iPat_%s_%s", project, trait.name[t]))
      iPat.Phenotype.View(
        myY = data.frame(taxa, updateY[, 3]),
        filename = sprintf("iPat_%s_%s", project, trait.name[t])
      )
      cat("Done\n")
    }
  }
}
