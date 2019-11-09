plotHistBV <- function(nameTableGEBV) {
  # Plotting
  dataGEBV = fread(paste0(nameTableGEBV, ".txt"))
  sizeBin = (range(dataGEBV$GEBV) %>% diff())/20
  plotGEBV = ggplot(data = dataGEBV, aes(x = GEBV, fill = ..density..)) +
    labs(title = sprintf("Distribution of GEBVs for %s", trait)) +
    geom_histogram(binwidth = sizeBin, colour = "black") +
    theme(text = element_text(size = 20))
  ggsave(plot = plotGEBV, file = paste0(nameTableGEBV, ".png"), height = 10, width = 14)
}

plotCV <- function(isRaw, nameTableRaw, nameTableAcc) {
  # Plotting
  if (isRaw) {
    dataRaw = fread(paste0(nameTableRaw, ".txt"))
    R2 = cor(dataRaw$obs, dataRaw$pre)^2
    plotRaw = ggplot(data = dataRaw, aes(x = obs, y = pre)) +
      stat_smooth(method = 'lm') +
      geom_point() +
      labs(title = sprintf("Observed Values Against Predictions Over Iterations (r2 = %.2f)", R2),
           x = "Observed Values in Testing Dataset",
           y = "Predicted Values in Testing Dataset") +
      theme(text = element_text(size = 20))
    ggsave(plot = plotRaw, file = paste0(nameTableRaw, ".png"), height = 10, width = 11)
  }
  dataAcc = fread(paste0(nameTableAcc, ".txt"))
  dataAcc2 = data.table(value = c(dataAcc$r, dataAcc$rmse),
                        type = rep(c("Correlation (r)", "RMSE"), each = nrow(dataAcc)))
  plotAcc = ggplot(data = dataAcc2, aes(x = value, fill = type)) +
    labs(title = "Distribution of Prediction Accuracy")  +
    geom_density(alpha = .5) +
    facet_grid(. ~ type, scales = "free") +
    theme(text = element_text(size = 20))
  ggsave(plot = plotAcc, file = paste0(nameTableAcc, ".png"), height = 10, width = 14)
}


`iPat.QQ` <-
  function(P.values, plot.type = "log_P_values", filename = "_", name.of.trait = "Trait", DPP = 50000, plot.style = "rainbow"){
    #Object: Make a QQ-Plot of the P-values
    #Options for plot.type = "log_P_values" and "P_values"
    #Output: A pdf of the QQ-plot
    #Authors: Alex Lipka and Zhiwu Zhang
    # Last update: May 9, 2011
    ##############################################################################################
    # Sort the data by the raw P-values
    P.values=P.values[!is.na(P.values)]
    P.values=P.values[P.values>0]
    P.values=P.values[P.values<=1]

    if(length(P.values[P.values>0])<1) return(NULL)
    N=length(P.values)
    DPP=round(DPP/4) #Reduce to 1/4 for QQ plot
    P.values <- P.values[order(P.values)]

    #Set up the p-value quantiles
    p_value_quantiles <- (1:length(P.values))/(length(P.values)+1)


    if(plot.type == "log_P_values")
    {
      log.P.values <- -log10(P.values)
      log.Quantiles <- -log10(p_value_quantiles)

      index=GAPIT.Pruning(log.P.values,DPP=DPP)
      log.P.values=log.P.values[index ]
      log.Quantiles=log.Quantiles[index]

      if(plot.style=="FarmCPU"){
        pdf(paste0(filename, "_QQ-Plot.pdf"),width = 5,height=5)
        par(mar = c(5,6,5,3))
      }
      if(plot.style=="rainbow"){
        pdf(paste0(filename, "_QQ-Plot.pdf"), width = 5,height=5)
        par(mar = c(5,6,5,3))
      }
      #Add conficence interval
      N1=length(log.Quantiles)
      ## create the confidence intervals
      c95 <- rep(NA,N1)
      c05 <- rep(NA,N1)
      for(j in 1:N1){
        i=ceiling((10^-log.Quantiles[j])*N)
        if(i==0)i=1
        c95[j] <- qbeta(0.95,i,N-i+1)
        c05[j] <- qbeta(0.05,i,N-i+1)
        #print(c(j,i,c95[j],c05[j]))
      }

      #CI Lines
      #plot(log.Quantiles, -log10(c05), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")
      #par(new=T)
      #plot(log.Quantiles, -log10(c95), xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, axes=FALSE, xlab="", ylab="",col="black")

      #CI shade
      plot(NULL, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), type="l",lty=5, lwd = 2, axes=FALSE, xlab="", ylab="",col="gray")
      index=length(c95):1
      polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col='gray',border=NA)

      #Diagonal line
      abline(a = 0, b = 1, col = "red",lwd=2)

      #data
      par(new=T)
      if(plot.style=="FarmCPU"){
        plot(log.Quantiles, log.P.values, cex.axis=1.1, cex.lab=1.3, lty = 1,  lwd = 2, col = "Black" ,bty='l', xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=""),pch=20)
      }
      if(plot.style=="rainbow"){
        plot(log.Quantiles, log.P.values, xlim = c(0,max(log.Quantiles)), ylim = c(0,max(log.P.values)), cex.axis=1.1, cex.lab=1.3, lty = 1,  lwd = 2, col = "Blue" ,xlab =expression(Expected~~-log[10](italic(p))),ylab = expression(Observed~~-log[10](italic(p))), main = paste(name.of.trait,sep=""))
      }

      dev.off()
    }


    if(plot.type == "P_values")
    {
      pdf(paste("QQ-Plot_", name.of.trait,".pdf" ,sep = ""))
      par(mar = c(5,5,5,5))
      qqplot(p_value_quantiles, P.values, xlim = c(0,1),
             ylim = c(0,1), type = "l" , xlab = "Uniform[0,1] Theoretical Quantiles",
             lty = 1, lwd = 1, ylab = "Quantiles of P-values from GWAS", col = "Blue",
             main = paste(name.of.trait,sep=" "))
      abline(a = 0, b = 1, col = "red")
      dev.off()
    }
    #print("GAPIT.QQ  accomplished successfully!")
  }

`iPat.Manhattan` <-
  function(GI.MP = NULL,GD=NULL, filename = "_", name.of.trait = "Trait",plot.type = "Genomewise",
           DPP=50000,cutOff=0.01,band=5,seqQTN=NULL,plot.style="Oceanic",CG=NULL,plot.bin=10^9){
    #Object: Make a Manhattan Plot
    #Options for plot.type = "Separate_Graph_for_Each_Chromosome" and "Same_Graph_for_Each_Chromosome"
    #Output: A pdf of the Manhattan Plot
    #Authors: Alex Lipka, Zhiwu Zhang, Meng Li and Jiabo Wang
    # Last update: Oct 10, 2016
    #Add r2 between candidata SNP and other markers in on choromosome
    ##############################################################################################
    #print("Manhattan ploting...")

    #do nothing if null input
    if(is.null(GI.MP)) return
    #if(is.null(GD)) return
    #print("Dimension of GI.MP")
    #print(dim(GI.MP))
    #print(head(GI.MP))
    #print(tail(GI.MP))
    #print(CG)

    #seqQTN=c(300,1000,2500)
    #Handler of lable paosition only indicated by negatie position
    position.only=F
    if(!is.null(seqQTN)){
      if(seqQTN[1]<0){
        seqQTN=-seqQTN
        position.only=T
      }

    }

    #if(is.null(GD)) print ("GD is not same dim as GM")

    borrowSlot=4
    GI.MP[,borrowSlot]=0 #Inicial as 0
    GI.MP[,5]=1:(nrow(GI.MP))
    if(!is.null(seqQTN))GI.MP[seqQTN,borrowSlot]=1


    GI.MP=matrix(as.numeric(as.matrix(GI.MP) ) ,nrow(GI.MP),ncol(GI.MP))
    GI.MP=GI.MP[order(GI.MP[,2]),]
    GI.MP=GI.MP[order(GI.MP[,1]),]
    #print("@@@@@")
    #print(dim(GD))
    #print(dim(GI.MP))
    if(!is.null(GD))
    {  if(ncol(GD)!=nrow(GI.MP))
    {print("GD does not mach GM in Manhattan !!!")
      return
    }}
    #print("!!")
    #GI.MP[,5]=1:(nrow(GI.MP))
    #print(head(GI.MP,20))
    #Remove all SNPs that do not have a choromosome, bp position and p value(NA)
    GI.MP <- GI.MP[!is.na(GI.MP[,1]),]
    GI.MP <- GI.MP[!is.na(GI.MP[,2]),]

    GI.MP <- GI.MP[!is.na(GI.MP[,3]),]

    #Retain SNPs that have P values between 0 and 1 (not na etc)
    GI.MP <- GI.MP[GI.MP[,3]>0,]
    GI.MP <- GI.MP[GI.MP[,3]<=1,]
    if(!is.null(GD)) GD=GD[,GI.MP[,3]<=1]
    #Remove chr 0 and 99
    GI.MP <- GI.MP[GI.MP[,1]!=0,]
    #GI.MP <- GI.MP[GI.MP[,1]!=99,]
    #print(dim(GI.MP))
    #print("Dimension of GI.MP after QC")
    #print(dim(GI.MP))

    numMarker=nrow(GI.MP)
    bonferroniCutOff=-log10(cutOff/numMarker)

    #Replace P the -log10 of the P-values
    if(!is.null(GD))
    {  if(ncol(GD)!=nrow(GI.MP))
    {print("GD does not mach GM in Manhattan !!!")
      return
    }}
    GI.MP[,3] <-  -log10(GI.MP[,3])
    index_GI=GI.MP[,3]>0
    GI.MP <- GI.MP[index_GI,]
    if(!is.null(GD)) GD=GD[,index_GI]

    GI.MP[,5]=1:(nrow(GI.MP))
    y.lim <- ceiling(max(GI.MP[,3]))
    chm.to.analyze <- unique(GI.MP[,1])
    #print(dim(GI.MP))
    #print(dim(GD))
    #print("name of chromosomes:")
    #print(chm.to.analyze)

    chm.to.analyze=chm.to.analyze[order(chm.to.analyze)]
    numCHR= length(chm.to.analyze)
    #GI.MP[,5]=1:(nrow(GI.MP))
    bin.mp=GI.MP[,1:3]
    bin.mp[,3]=0 # for r2
    bin.mp[,1]=as.numeric(as.vector(GI.MP[,2]))+as.numeric(as.vector(GI.MP[,1]))*(10^(max(GI.MP[,1])+1))


    #as.numeric(as.vector(GP[,3]))+as.numeric(as.vector(GP[,2]))*MaxBP

    bin.mp[,2]=floor(bin.mp[,1]/plot.bin)
    if(!is.null(GD)) X=GD

    #print(head(bin.mp))
    #Chromosomewise plot
    if(plot.type == "Chromosomewise"&!is.null(GD))
    {
      #print("Manhattan ploting Chromosomewise")
      GI.MP=cbind(GI.MP,bin.mp)
      pdf(paste0(filename,"_Manhattan.Plot.Chromosomewise.pdf"), width = 10)
      #par(mar = c(5,5,4,3), lab = c(8,5,7))
      layout(matrix(c(1,1,2,1,1,1,1,1,1),3,3,byrow=TRUE), c(2,1), c(1,1), TRUE)
      for(i in 1:numCHR)
      {
        #Extract SBP on this chromosome
        subset=GI.MP[GI.MP[,1]==chm.to.analyze[i],]
        # print(head(subset))
        subset[,1]=1:(nrow(subset))
        #sub.bin.mp=bin.mp[GI.MP[,1]==chm.to.analyze[i],]
        #subset=cbind(subset,sub.bin.mp)
        sig.mp=subset[subset[,3]>bonferroniCutOff,]
        sig.index=subset[,3]>bonferroniCutOff ### index of significont SNP


        num.row=nrow(sig.mp)
        if(!is.null(dim(sig.mp)))sig.mp=sig.mp[!duplicated(sig.mp[,7]),]
        num.row=nrow(sig.mp)
        if(is.null(dim(sig.mp))) num.row=1
        bin.set=NULL
        r2_color=matrix(0,nrow(subset),2)
        #r2_color
        print(paste("select ",num.row," candidate gene in ",i," chromosome ",sep="") )
        #print(sig.mp)
        if(length(unique(sig.index))==2)
        {
          for(j in 1:num.row)
          {   sig.mp=matrix(sig.mp,num.row,8)

          #print(sig.mp[j,7])
          #print(unique(subset[,7]))
          bin.store=subset[which(subset[,7]==sig.mp[j,7]),]
          if(is.null(dim(bin.store)))
          {subset[which(subset[,7]==sig.mp[j,7]),8]=1
          next
          }
          bin.index=unique(bin.store[,5])
          subGD=X[,bin.store[,5]]
          #print(dim(bin.store))
          if(is.null(CG))candidata=bin.store[bin.store[,3]==max(bin.store[,3]),5]
          if(length(candidata)!=1)candidata=candidata[1]

          for (k in 1:ncol(subGD))
          {
            r2=cor(X[,candidata],subGD[,k])^2
            #print(r2)
            bin.store[k,8]=r2

          }
          #print(bin.store)
          #r2_storage[is.na(r2_storage)]=0
          #print(bin.store)
          subset[bin.store[,1],8]=bin.store[,8]
          #print()
          }###end for each sig.mp
          #sub.bin.mp=bin.mp[subset[,3]>bonferroniCutOff,]
          #print(head(bin.set))

        }###end if empty of sig.mp
        #print("@@@@@@@@@@@@@@@@")
        rm(sig.mp,num.row)
        #print(head(subset))
        #print(head(subset))
        #print(dim(X))
        y.lim <- ceiling(max(subset[,3]))+1  #set upper for each chr
        if(length(subset)>3){
          x <- as.numeric(subset[,2])/10^(6)
          y <- as.numeric(subset[,3])
        }else{
          x <- as.numeric(subset[2])/10^(6)
          y <- as.numeric(subset[3])
        }

        ##print(paste("befor prune: chr: ",i, "length: ",length(x),"max p",max(y), "min p",min(y), "max x",max(x), "Min x",min(x)))
        n_col=10
        r2_color[,2]=subset[,8]
        do_color=colorRampPalette(c("orangeRed", "blue"))(n_col)
        #Prune most non important SNPs off the plots
        order=order(y,decreasing = TRUE)
        y=y[order]
        x=x[order]
        r2_color=r2_color[order,]
        index=GAPIT.Pruning(y,DPP=round(DPP/numCHR))
        x=x[index]
        y=y[index]
        r2_color=r2_color[index,]
        r2_color[which(r2_color[,2]<=0.2),2]=do_color[n_col]
        r2_color[which(r2_color[,2]<=0.4&r2_color[,2]>0.2),2]=do_color[n_col*0.8]
        r2_color[which(r2_color[,2]<=0.6&r2_color[,2]>0.4),2]=do_color[n_col*0.6]
        r2_color[which(r2_color[,2]<=0.8&r2_color[,2]>0.6),2]=do_color[n_col*0.4]
        r2_color[which(r2_color[,2]<=1&r2_color[,2]>0.8),2]=do_color[n_col/n_col]

        #print(unique(r2_color[,2]))

        ##print(paste("after prune: chr: ",i, "length: ",length(x),"max p",max(y), "min p",min(y), "max x",max(x), "Min x",min(x)))

        par(mar=c(0,0,0,0))
        par(mar=c(5,5,2,1),cex=0.8)

        plot(y~x,type="p", ylim=c(0,y.lim), xlim = c(min(x), max(x)),
             col = r2_color[,2], xlab = expression(Base~Pairs~(x10^-6)),
             ylab = "-Log Base 10 p-value", main =       paste("Chromosome",chm.to.analyze[i],sep=" "),
             cex.lab=1.6,pch=21,bg=r2_color[,2])

        abline(h=bonferroniCutOff,col="forestgreen")
        ##print("manhattan plot (chr) finished")
        #layout.show(nf)
        #provcol <-c("darkblue","cyan","green3","brown1","brown1")
        #provcol <-heat.colors(50)
        #par(mar=c(0,0,0,0))
        par(mar=c(15,5,6,5),cex=0.5)

        barplot(matrix(rep(1,times=n_col),n_col,1),beside=T,col=do_color,border=do_color,axes=FALSE,)
        #legend(x=10,y=2,legend=expression(R^"2"),,lty=0,cex=1.3,bty="n",bg=par("bg"))
        axis(3,seq(11,1,by=-2),seq(0,1,by=0.2))

      }# end plot.type == "Chromosomewise"&!is.null(GD)
      dev.off()

      print("manhattan plot on chromosome finished")
    } #Chromosomewise plot


    #Genomewise plot
    if(plot.type == "Genomewise")
    {
      #print("Manhattan ploting Genomewise")
      #Set corlos for chromosomes
      #nchr=max(chm.to.analyze)
      nchr=length(chm.to.analyze)

      #Set color schem
      ncycle=ceiling(nchr/band)
      ncolor=band*ncycle
      #palette(rainbow(ncolor+1))
      cycle1=seq(1,nchr,by= ncycle)
      thecolor=cycle1
      for(i in 2:ncycle){thecolor=c(thecolor,cycle1+(i-1))}
      col.Rainbow=rainbow(ncolor+1)[thecolor]
      col.FarmCPU=rep(c("#CC6600","deepskyblue","orange","forestgreen","indianred3"),ceiling(numCHR/5))
      col.Rushville=rep(c("orangered","navyblue"),ceiling(numCHR/2))
      col.Congress=rep(c("deepskyblue3","firebrick"),ceiling(numCHR/2))
      col.Ocean=rep(c("steelblue4","cyan3"),ceiling(numCHR/2))
      col.PLINK=rep(c("gray10","gray70"),ceiling(numCHR/2))
      col.Beach=rep(c("turquoise4","indianred3","darkolivegreen3","red","aquamarine3","darkgoldenrod"),ceiling(numCHR/5))
      #col.Oceanic=rep(c( '#EC5f67',  '#F99157',  '#FAC863',  '#99C794',  '#5FB3B3',  '#6699CC',  '#C594C5',  '#AB7967'),ceiling(numCHR/8))
      #col.Oceanic=rep(c( '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5',  '#AB7967'),ceiling(numCHR/6))
      col.Oceanic=rep(c(  '#EC5f67',    '#FAC863',  '#99C794',    '#6699CC',  '#C594C5'),ceiling(numCHR/5))
      col.cougars=rep(c(  '#990000',    'dimgray'),ceiling(numCHR/2))

      if(plot.style=="Rainbow")plot.color= col.Rainbow
      if(plot.style =="FarmCPU")plot.color= col.Rainbow
      if(plot.style =="Rushville")plot.color= col.Rushville
      if(plot.style =="Congress")plot.color= col.Congress
      if(plot.style =="Ocean")plot.color= col.Ocean
      if(plot.style =="PLINK")plot.color= col.PLINK
      if(plot.style =="Beach")plot.color= col.Beach
      if(plot.style =="Oceanic")plot.color= col.Oceanic
      if(plot.style =="cougars")plot.color= col.cougars

      #FarmCPU uses filled dots
      mypch=1
      if(plot.style =="FarmCPU")mypch=20

      GI.MP <- GI.MP[order(GI.MP[,2]),]
      GI.MP <- GI.MP[order(GI.MP[,1]),]

      ticks=NULL
      lastbase=0

      #print("Manhattan data sorted")
      #print(chm.to.analyze)

      #change base position to accumulatives (ticks)
      for (i in chm.to.analyze)
      {
        index=(GI.MP[,1]==i)
        ticks <- c(ticks, lastbase+mean(GI.MP[index,2]))
        GI.MP[index,2]=GI.MP[index,2]+lastbase
        lastbase=max(GI.MP[index,2])
      }

      #print("Manhattan chr processed")
      #print(length(index))
      #print(length(ticks))
      #print((ticks))
      #print((lastbase))

      x0 <- as.numeric(GI.MP[,2])
      y0 <- as.numeric(GI.MP[,3])
      z0 <- as.numeric(GI.MP[,1])
      position=order(y0,decreasing = TRUE)
      index0=GAPIT.Pruning(y0[position],DPP=DPP)
      index=position[index0]

      x=x0[index]
      y=y0[index]
      z=z0[index]

      #Extract QTN
      QTN=GI.MP[which(GI.MP[,borrowSlot]==1),]

      #Draw circles with same size and different thikness
      size=1 #1
      ratio=10 #5
      base=1 #1
      themax=ceiling(max(y))
      themin=floor(min(y))
      wd=((y-themin+base)/(themax-themin+base))*size*ratio
      s=size-wd/ratio/2

      #print("Manhattan XY created")
      ####xiaolei update on 2016/01/09
      if(plot.style =="FarmCPU"){
        pdf(paste0(filename,"_Manhattan.Plot.Genomewise.pdf"), width = 13, height=5.75)
      }else{
        pdf(paste0(filename,"_Manhattan.Plot.Genomewise.pdf"), width = 13, height=5.75)
      }
      par(mar = c(3,6,5,1))
      plot(y~x,xlab="",ylab=expression(-log[10](italic(p))) ,
           cex.axis=1.5, cex.lab=2, ,col=plot.color[z],axes=FALSE,type = "p",pch=mypch,lwd=wd,cex=s+.3,main = paste(name.of.trait,sep="      "),cex.main=2.5)

      #Label QTN positions
      if(is.vector(QTN)){
        if(position.only){abline(v=QTN[2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[2], QTN[3], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[2], QTN[3], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
        }
      }else{
        if(position.only){abline(v=QTN[,2], lty = 2, lwd=1.5, col = "grey")}else{
          points(QTN[,2], QTN[,3], type="p",pch=21, cex=2,lwd=1.5,col="dimgrey")
          points(QTN[,2], QTN[,3], type="p",pch=20, cex=1,lwd=1.5,col="dimgrey")
        }
      }

      #Add a horizontal line for bonferroniCutOff
      abline(h=bonferroniCutOff,col="forestgreen")

      #Set axises
      axis(1, at=ticks,cex.axis=1.5,labels=chm.to.analyze,tick=F)
      axis(2, at=1:themax,cex.axis=1.5,labels=1:themax,tick=F)

      box()
      palette("default")
      dev.off()
      #print("Manhattan done Genomewise")

    } #Genomewise plot

    print("GAPIT.Manhattan accomplished successfully!zw")
} #end of GAPIT.Manhattan

`iPat.Genotype.View` <-function(myGD = NULL, filename = "_"){
  # Object: Analysis for Genotype data:Distribution of SNP density,Accumulation,Moving Average of density,result:a pdf of the scree plot
  # myG:Genotype data
  # chr: chromosome value
  # w1_start:Moving Average windows Start Position
  # w1_end:Moving Average windows End Position
  # mav1:Moving Average set value length
  # Authors: You Tang and Zhiwu Zhang
  # Last update: March 11, 2016
  ##############################################################################################

  #heterozygosity of individuals and SNPs (By Zhiwu Zhang)
  #print("Heterozygosity of individuals and SNPs (By Zhiwu Zhang)")
  X=myGD[,-1]
  H=1-abs(X-1)
  het.ind=apply(H,1,mean)
  het.snp=apply(H,2,mean)
  ylab.ind=paste("Frequency (out of ",length(het.ind)," individuals)",sep="")
  ylab.snp=paste("Frequency (out of ",length(het.snp)," markers)",sep="")
  pdf(paste0(filename, "_heterozygosity.pdf"), width = 10, height = 6)
  par(mfrow=c(1,2),mar=c(5,5,1,1)+0.1)
  hist(het.ind,col="gray", main="",ylab=ylab.ind, xlab="Heterozygosity of individuals")
  hist(het.snp,col="gray", main="",ylab=ylab.snp, xlab="Heterozygosity of markers")
  dev.off()
}


`iPat.Phenotype.View` <-function(myY = NULL, filename = "_"){
  # Object: Analysis for Phenotype data:Distribution of density,Accumulation,result:a pdf of the scree plot
  # myY:Phenotype data

  # Authors: You Tang
  # Last update: Sep 7, 2015
  ##############################################################################################
  print("GAPIT.Phenotype.View in press...")
  if(is.null(myY)){stop("Validation Invalid. Please select read valid Phenotype flies  !")}

  y<-myY[!is.na(myY[,2]),2]
  obs<-as.matrix(y)

  traitname=colnames(myY)[2]

  pdf(paste0(filename, "_phnotype_view.pdf"), width = 10, height = 6)
  par(mar = c(5,5,5,5))

  par(mfrow=c(2,2))
  plot(obs,pch=1)
  #hist(obs)
  hist(obs,xlab="Density",main="",breaks=12, cex.axis=1,col = "gray")
  boxplot(obs)
  plot(ecdf(obs),col="red",bg="lightgray",xlab="Density",ylab="Accumulation",main="")

  dev.off()
}

`Blink.LDRemoveDivided`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL,l=NULL){
  #Objects: LD remove, especially length(Porder)>10000
  #Authors: Yao Zhou
  #Last update: 08/15/2016
  seqQTN=NULL
  lp=length(Porder)
  k=ceiling(lp/l)
  GDneo=as.matrix(GDneo)
  if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
  if(orientation=="col"){
    n=nrow(GDneo)
    if(bound){
      GDneo=GDneo[sample(n,200,replace=F),]
    }
  }else{
    n=ncol(GDneo)
    if(bound){
      GDneo=GDneo[,sample(n,200,replace=F)]
    }
    GDneo=t(GDneo)
  }
  for(i in 1:k){
    bottom=(i-1)*l+1
    up=l*i
    if(up>lp) up = lp
    Porderb=Porder[bottom:up]
    GDneob = GDneo[,bottom:up]
    seqQTNs=Blink.LDRemoveBlock(GDneo=GDneob,LD=LD,Porder=Porderb,orientation="col",model=model)
    seqQTN=append(seqQTN,seqQTNs)
  }
  rm(GDneob,seqQTNs,Porderb)
  return(seqQTN)
}

`Blink.LDRemoveBlock`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
  #`Blink.LDRemove`<-function(GDneo=NULL,LD=NULL,Porder=NULL,bound=FALSE,model="A",orientation=NULL){
  #Objects: Calculate LD and remove the correlated SNPs
  #Authors: Yao Zhou
  #Last Update:  03/03/16
  if (model=="D"){
    GDneo=1-abs(GDneo-1)
  }

  GDneo=as.matrix(GDneo)
  if(min(ncol(GDneo),nrow(GDneo))<201) bound=FALSE
  if(orientation=="col"){
    n=nrow(GDneo)
    if(bound){
      GDneo=GDneo[sample(n,200,replace=F),]
    }
  }else{
    n=ncol(GDneo)
    if(bound){
      GDneo=GDneo[,sample(n,200,replace=F)]
    }
    GDneo=t(GDneo)
  }

  corr=cor(GDneo)
  corr[is.na(corr)]=1
  corr[abs(corr)<=LD]=0
  corr[abs(corr)>LD]=1
  Psort=as.numeric(matrix(1,1,ncol(corr)))
  for(i in 2:ncol(corr)){
    p.a=Psort[1:(i-1)]
    p.b=as.numeric(corr[1:(i-1),i])
    index=(p.a==p.b)
    index[(p.a==0)&(p.b==0)]=FALSE
    if(sum(index)!=0) Psort[i]=0
  }
  seqQTN=Porder[Psort==1]
  return(seqQTN)
}

