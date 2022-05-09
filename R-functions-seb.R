### Create a new invisible environment for all the functions to go in so it doesn't clutter your workspace.
.sebenv <- new.env()
################################### 2nd
.sebenv$headm <- function(x, n=6) x[1:n, 1:n]

# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
.sebenv$compareNA <- function(v1,v2) {
    same <- (v1 == v2) | (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}

.sebenv$di <- function() X11.options(display=scan(file="~/.display",what=character()))

.sebenv$label2number <- function(myvec) {
  tmp <- factor(myvec)
  levels(tmp) <- 1:length(levels(tmp))
  return(as.integer(tmp))
}
# label2color(c("a","b","a"))
# [1] 1 2 1


.sebenv$drawhclust2 <- function(pan){
  plot(pan$data,type="n")
  if (!pan$lbvalue %in% names(pan$hlist)) {
    pan$distmat=pan$dists[[pan$lbvalue]]$dist
    tree=hclust(as.dist(1-pan$distmat/max(pan$distmat)),method="single")
    pan$hlist[[pan$lbvalue]]=tree
  }
  if (pan$xrns) points(pan$dists[[pan$lbvalue]]$resampled_points, col=16,pch=4)
  else text(pan$data,col=cutree(pan$hlist[[pan$lbvalue]],pan$slidervalue), labels = as.character(1:nrow(pan$data)))
  return(pan)			# panel has to be returned
}

.sebenv$drawtree <- function(pan){
  text(pan$data,col=cutree(pan$hlist[[pan$lbvalue]],pan$slidervalue), labels = as.character(1:nrow(pan$data)))
  return(pan)			# panel has to be returned
}

.sebenv$oncheck <- function(pan) {
  if (pan$xrns) points(pan$dists[[pan$lbvalue]]$resampled_points, col=16,pch=4)
  else { plot(pan$data,type="n"); drawtree(pan) }
  pan
}

#df$yw}i.sebenv$pi_porm(pi)

.sebenv$plothclust2 <- function(dists,data){
  panel <- rp.control(data=data,dists=dists,slidervalue=1,hlist=list(),distmat=NULL,xrns=FALSE)	# don't take "a"
  rp.listbox(panel, lbvalue,action=drawhclust2, names(dists))
  rp.slider(panel,slidervalue,1,20,drawtree,showvalue=TRUE,res=1)
  rp.checkbox(panel, xrns, action=oncheck)
}

.sebenv$.ls.objects <- function(pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
                                       fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
                             capture.output(print(object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
                      as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[order.by], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
.sebenv$lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

.sebenv$wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
  options(width=as.integer(howWide))
}

.sebenv$sumfilter <- function(x,n=500){filter(x,rep(1,n), sides=2)}

.sebenv$meanfilter <- function(x,n=500){filter(x,rep(1/n,n), sides=2)}

###cross overlapping table types in of g-range data
#given: GRanges object (here genes) with annotation in say type column
#output: table indicating how many of those types overlap
.sebenv$cross_overlap_grange <- function(genes) {
  types <- levels(genes$type)
  m <- matrix(NA,ncol=length(types),nrow=length(types))
  colnames(m) <- rownames(m) <- types
  for (type in types) {
    for (type2 in types) {
      locivsannotation <- findOverlaps(genes[genes$type==type,],genes[genes$type==type2,])
      m[type, type2] <- sum(!duplicated(queryHits(locivsannotation)))
    }
  }
  return(m)
}

.sebenv$makepumatrix <- function(reads,weight=1L) {
  treatments <- paste(rep(21:24,2),rep(c("+","-"),each=4))
  if (class(reads)=="GRanges" & weight==1L) {warning("GRanges-Objects (read in by segmentSeq::readBAM) usually come with weights since they tend to be produced by readBAM, please double check")}
  pumatrix <- matrix(0, ncol=9, nrow=seqlengths(reads))
  pumatrix[,2] <- as.integer(coverage(reads[width(reads)==21 & strand(reads)=="+"],weight=weight)[[1]])
  pumatrix[,3] <- as.integer(coverage(reads[width(reads)==22 & strand(reads)=="+"],weight=weight)[[1]])
  pumatrix[,4] <- as.integer(coverage(reads[width(reads)==23 & strand(reads)=="+"],weight=weight)[[1]])
  pumatrix[,5] <- as.integer(coverage(reads[width(reads)==24 & strand(reads)=="+"],weight=weight)[[1]])
  pumatrix[,6] <- -as.integer(coverage(reads[width(reads)==21 & strand(reads)=="-"],weight=weight)[[1]])
  pumatrix[,7] <- -as.integer(coverage(reads[width(reads)==22 & strand(reads)=="-"],weight=weight)[[1]])
  pumatrix[,8] <- -as.integer(coverage(reads[width(reads)==23 & strand(reads)=="-"],weight=weight)[[1]])
  pumatrix[,9] <- -as.integer(coverage(reads[width(reads)==24 & strand(reads)=="-"],weight=weight)[[1]])
  colnames(pumatrix)<-c("c1",treatments)
  return(pumatrix)
}

.sebenv$makepumatrixgg <- function(aln,weight=1L,sizes=list("21nt"=21,"22nt"=22,"23nt"=23,"24nt"=24),strands=list("+strand"="+","-strand"="-"),samplenames=colnames(aln@data),normalize=FALSE) {
  if (class(aln)=="GRanges" & weight==1L) {warning("GRanges-Objects (read in by segmentSeq::readBAM) usually come with weights since they tend to be produced by readBAM, please double check")}
  pumatrix <- data.frame(coverage = numeric(0), position = numeric(0), chromosome = character(0), strand = character(0), size = numeric(0), sample = character(0))
  if (normalize) {aln@data <- t(t(aln@data)/c(aln@libsizes/mean(aln@libsizes)))}
  for (i in 1:ncol(aln)) {
    reads <- aln@alignments
    reads$counts <- aln@data[,i]
    #tmpalignData@alignments$counts <- aln@data[,i]
    #reads <- tmpalignData@alignments
    for (chr in names(seqlengths(aln@alignments))) {
      for (size in sizes) {
        for (st in strands) {
          cov <- as.integer(coverage(reads[width(reads)%in%size & strand(reads)%in%st],weight=weight)[[chr]])
          pumatrixslice <- data.frame(coverage=cov, position=1:length(cov),chromosome=rep(chr,length(cov)), strand=rep(st,length(cov)),size=rep(size,length(cov)),sample=rep(samplenames[i],length(cov)))
          pumatrix <- rbind(pumatrix,pumatrixslice)
        }
      }
    }
  }
  pumatrix$size <- as.factor(pumatrix$size)
  return(pumatrix)
}
#pmatrix <- makepumatrixgg(readsCYM[exp(CDCYMtmmL@posteriors[,"formaldehydeTreatvsnoForm"]) > 0.9 & CDCYMtmmL@orderings$formaldehydeTreatvsnoForm=="1>2" ,])

#pmatrix[pmatrix$strand=="-",]$coverage <- -pmatrix[pmatrix$strand=="-",]$coverage
#pumatrixgg2 <- pmatrix[pmatrix$coverage!=0,]
#ggplot(pumatrixgg2, aes(x=position,y=coverage,fill=size)) + geom_bar(stat="identity",position = "identity") + facet_grid(sample~chromosome,scales = "free") + theme_bw() + theme(strip.text.y = element_text(size = 7))

# input 2 vectors (of which elementwise M and A is required)
# output: data.frame with m and a column
.sebenv$logma <- function(x,y) {
  a <- 0.5*(log2(x*y))
  m <- log2(x/y)
  cbind(a=a,m=m)
}


.sebenv$findYLim <- function( pumatrix){
  ymax <- 0
  ymin <- 0
  fileMax <- max(pumatrix[,2:length(pumatrix[1,])])
  fileMin <- min(pumatrix[,2:length(pumatrix[1,])])
  if(fileMax > ymax){ymax <- fileMax}
  if(fileMin  < ymin ){ymin <- fileMin}
  c( ymin, ymax)
}

.sebenv$plotGraph <- function( pumatrix,ylims, treatments,title,linecols=colors()[c(554,256,76,28)],legend = c(21:24)){
  ytitle <- "count of sRNAs covering a base"
  xtitle <- paste("distance along ", "sequence", sep = " ")
  #plot.new()
  #PDF
  #pdf(pdffilename)
  colCount <- length(pumatrix[1,]) - 1
  matplot(1:nrow(pumatrix),pumatrix[,2:length(pumatrix[1,])],
          type = "l",
          lty = 1,
          main = title,
          col = linecols,
          ylab = "",
          xlab = "",
          xaxt = "n",
          yaxt = "n",
          #cex.axis = 2.5,
          ylim = ylims,
          cex.main = 0.7,
          cex.lab = 0.5,
          )
  #X
  #xlab <- c("beg","mid","end")
  #xat <- c(1,(chrlength %/% 2),chrlength)
  xat <- pretty(1:nrow(pumatrix), n = 2, min.n = 2, u5.bias = 0.4,high.u.bias = 1)
  xkb <- xat/1000
  xkb <- format(xkb, digits = 2)
  #xlab <- c(as.character(1),as.character(chrlength %/% 2),as.character(chrlength))
  xlab <- as.character(xkb)
  message("chrlen ", nrow(pumatrix), "\n")
  message("xat ",xat, "LABELS ",xkb, "\n")
  #prev
  axis(1, at = xat, labels = xlab, cex.axis = 2.0, padj = 1)
  yat <- pretty(ylims, n = 3,min.n = 3,high.u.bias = 1)
  yneg <- yat * -1
  ynoneg <- ifelse(yat >= 0, yat, yneg)
  ynoneg <- sprintf("%3.0f", ynoneg)
  ynoneg <- format(ynoneg, digits = 2, trim = TRUE, justify = c("left"))
  ylab <- as.character(ynoneg)
  if(length(unique(ylab)) == 1){
    yat <- c(-1,0,1)
    ylab <- c("1","0","1")
  }
  message("ylims ", ylims, "\n")
  message("yat ", yat, "LABELS ",ynoneg, "\n")
  axis(2, at = yat, labels = ylab, cex.axis = 2.0)
  #line
  abline( h= 0, lty = 1, col = "black")
  #don't want this but keep for now to make sure colours right way round
  legend("topright",legend = legend, lty = 1,col =  unique(linecols))
}

#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
.sebenv$multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),byrow = TRUE)
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

###these are some auxiliary functions, needed for later calculation###
# takes "lociData"-object (or a subset thereof) along with a genome (DNAStringSet-object) and extracts the for each loci
# the sequences based on its coordinates, followed by the calculation of sequence similarity for each loci pairing
# the result is saved in the $both slot of the resulting object
# "lociData": contains a list of loci (with coordinates and likelihoods) generated by segmentSeq along with the number of reads mapping to it for each library
# example call:
# unqall_24_3 <- subset2network2(subset=subsetall,cD=classSegLike99,genome=genomemod,ngramwindow=24,maxmismatches=3,tb.end=2)

#!!!##toms new data: /home/tjh48/Code/segmentMap_II/classSegLike.RData
.sebenv$subset2network2 <- function(subset,cD,genome,ngramwindow=18,maxmismatches=2,tb.end=4,mode="both") {
  ## change to v1: in cD object, the annoation slot is no changed to coordinate slot
  ## also the coordinates are now GRanage-objects
  mode <- (match.arg(mode))
  locisubset <- cD@coordinates[subset,]
  #for coloring afterwards, not a necisity
  likelihoodsubset <- cD@locLikelihoods[subset,]
  #naming the loci according to location
  namessubset <- paste(seqnames(locisubset),start(locisubset),end(locisubset),sep="_")
  #extracting sequences information from subset loci
  seqs <- DNAStringSet(vector(mode = "character", length = length(locisubset)))
  for (locus in 1:length(locisubset)) {
    current <- locisubset[locus,]
    seqs[[locus]] <- genome[[as.character(seqnames(current))]][start(current):end(current)] # this might be rewritten more elegantly...
  }
  names(seqs) <- namessubset
  im <- create_interaction_matrix(seqs,ngramwindow=ngramwindow,maxmismatches=maxmismatches,tb.end=tb.end)
  if (mode=="both") im$both <- im$im +im$imr
  im$subset <- subset
  im$seqs <- seqs
  im$like <- likelihoodsubset
  im$locisubset <- locisubset
  im
}

## function to create interaction matrices
## it takes a list of Sequences (BiostringSet-object) and calculates for each pair of this list if they have a common sequence
## based on some criteria (i.e. 24 conscutive matches with at most 3 mismatches).
## 0 in the matrix means not matching this criteria
## 1 or higher means the numerber of stretches matching this criteria
create_interaction_matrix <- function(sequences,sequencesrev=reverseComplement(sequences),ngramwindow=24,maxmismatches=1,tb.end=12) {
  nseqs <- length(sequences)
  interaction_matrix <- matrix(nr=nseqs,ncol=nseqs)
  interaction_matrix_rev <- matrix(nr=nseqs,ncol=nseqs)
  for (sub in 1:nseqs) {
    subject <- sequences[[sub]]
    seqsdict <- PDict(Views(subject, start=1:(length(subject)-ngramwindow+1), end=ngramwindow:length(subject)),tb.start=1,tb.end=tb.end)
    counts <- vcountPDict(seqsdict,sequences,max.mismatch=maxmismatches)
    interaction_matrix[sub,] <- colSums(counts)
    countsrev <- vcountPDict(seqsdict,sequencesrev,max.mismatch=maxmismatches)
    interaction_matrix_rev[sub,] <- colSums(countsrev)
  }
  colnames(interaction_matrix) <- rownames(interaction_matrix) <- colnames(interaction_matrix_rev) <- rownames(interaction_matrix_rev) <- names(sequences)
  list(im=interaction_matrix,imr=interaction_matrix_rev)
}
.sebenv$create_interaction_matrix  <- create_interaction_matrix
rm(create_interaction_matrix )

# compression
lz <-function(seq2, D = c("A", "C", "G", "T")) {
  seq <- strsplit(as.character(seq2),"")[[1]]
  index <- c()
  pos <- 1
  #seq <- sample(c("A", "C", "G", "T"), size = 200, replace = TRUE)
  # algorithm
  while(pos <= length(seq)) {
    maxlen <- max(nchar(D))
    posblock <- sapply(maxlen:1 - 1, function(x) paste(seq[pos:(pos + x)], collapse = ""))
    selblock <- min(which(posblock %in% D))
    index <- c(index, which(D == posblock[selblock]))
    pos <- pos + nchar(posblock[selblock])
    D <- c(D, paste(posblock[selblock], seq[pos], sep = ""))
  }
  #return(list(index,D,length(index)))
  length(index)
}
.sebenv$lz  <- lz
rm(lz )

methoverlaptab <- function(gr, gr2, lib) {
  #this function takes a predifined set of intervals (gr object, such as genes, TEs, promoters or bins)
  #gr needs a following mcols columns: type, ID
  #and overlaps it with a bismarck object containg methyl counts (readsM and readsN).
  #counts for overlapping C are summed up for each element in gr, returning a DataFrame containing
  #the following columns
  #Name,seqname, position, library, context, type
  contexts <- c("CG","CHG","CHH")
  annot <- DataFrame(start=start(gr),seqname=seqnames(gr),mcols(gr)[,c("ID","type")])
  libcolumn <- rep(lib, length(gr))
  methtaball <- NULL
  for (context in contexts) {
    methtab <- annot
    grtmp <-gr2[gr2$context==context]
    overlap <- findOverlaps(gr,grtmp,ignore.strand=TRUE)
    readsM <- readsN <- rep(NA,length(gr))
    contextcolumn <- rep(context,length(gr))
    tmp <- tapply(grtmp[subjectHits(overlap),]$readsM,queryHits(overlap),sum)
    readsM[as.numeric(names(tmp))] <- tmp
    tmp <- tapply(grtmp[subjectHits(overlap),]$readsN,queryHits(overlap),sum)
    readsN[as.numeric(names(tmp))] <- tmp
    methtab$lib <- libcolumn
    methtab$context <- contextcolumn
    methtab$readsM <- readsM
    methtab$readsN <- readsN
    methtab$pctmeth <- readsM/readsN
    methtaball <- rbind(methtaball,methtab)
  }
  methtaball
}
.sebenv$methoverlaptab  <- methoverlaptab
rm(methoverlaptab )

drawhclust <- function(pan){
  require(rpanel)
  plot(pan$data,type="n")
  text(pan$data,col=cutree(pan$h,pan$slidervalue), labels = as.character(1:nrow(pan$data)))
  return(pan)			# panel has to be returned
}
.sebenv$drawhclust  <- drawhclust
rm(drawhclust )

plothclust=function(hobj,data){
  panel <- rp.control(data=data,slidervalue=1,h=hobj)	# don't take "a"
  rp.slider(panel,slidervalue,1,20,drawhclust,showvalue=TRUE,res=1)
}
.sebenv$plothclust <- plothclust
rm(plothclust)

################################### 2nd
drawhclust3 <- function(pan){
  plot(pan$data,type="n")
  if (!pan$lbvalue %in% names(pan$hlist)) {
    pan$distmat=pan$dists[[pan$lbvalue]]$dist
    tree=hclust(as.dist(1-pan$distmat/max(pan$distmat)),method="single")
    pan$hlist[[pan$lbvalue]]=tree
  }
  if (pan$xrns) points(pan$dists[[pan$lbvalue]]$resampled_points, col=16,pch=4)
  else text(pan$data,col=cutree(pan$hlist[[pan$lbvalue]],pan$slidervalue), labels = as.character(1:nrow(pan$data)))
  return(pan)			# panel has to be returned
}
.sebenv$drawhclust3  <- drawhclust3
rm(drawhclust3 )

drawtree <- function(pan){
  text(pan$data,col=cutree(pan$hlist[[pan$lbvalue]],pan$slidervalue), labels = as.character(1:nrow(pan$data)))
  return(pan)			# panel has to be returned
}
.sebenv$drawtree  <- drawtree
rm(drawtree )

.sebenv$oncheck <- function(pan) {
  if (pan$xrns) points(pan$dists[[pan$lbvalue]]$resampled_points, col=16,pch=4)
  else { plot(pan$data,type="n"); drawtree(pan) }
  pan
}

#http://en.wikipedia.org/wiki/Fisher's_method
.sebenv$fisheromnibus <- function(x) {
  number_conditions <- length(x)
  X2 <- -2 * sum(log(x))
  pchisq( X2, 2*number_conditions, lower.tail=FALSE)
}

#fisher omnibus: Fabian: verwende jetzt andere Methode, welche das in eine Normalverteilung umwandelt
.sebenv$stouffer <- function(x) {
  nr_tests <- length(x)
  statistic <- sum(qnorm(x,lower.tail=FALSE))/sqrt(nr_tests)
  pval <- pnorm(statistic,lower.tail=FALSE)
  return(pval)
}
#Stouffer et al 1949 The American Soldier, Vol 1: Adjustment during Army Life, Princeton University Press, Princeton
#
#Vorteile siehst du bei folgenden Fällen
#
#Symmetrie
#Stouffer:c(0.99,0.99) == 1-c(0.01,0.01)
#fisher:c(0.99,0.99) != 1-c(0.01,0.01)
#
#Gewichtung
#stouffer: c(0.999,0.001) = 0.5
#fisher: c(0.999,0.001) = 0.008


.sebenv$plothclust2=function(dists,data){
  panel <- rp.control(data=data,dists=dists,slidervalue=1,hlist=list(),distmat=NULL,xrns=FALSE)	# don't take "a"
  rp.listbox(panel, lbvalue,action=drawhclust3, names(dists))
  rp.slider(panel,slidervalue,1,20,drawtree,showvalue=TRUE,res=1)
  rp.checkbox(panel, xrns, action=oncheck)
}

.sebenv$plotMACDreturn <- function(cD, samplesA, samplesB, normaliseData = TRUE, scale = NULL,
                            xlab = "A", ylab = "M", ...)
{
  if (length(dim(cD)) > 2)
    stop("This function is currently only applicable to 2-dimensional countData objects.")
  if (is.character(samplesA)) {
    Asamps <- which(as.character(cD@replicates) %in% samplesA)
    if (!all(samplesA %in% cD@replicates))
      Asamps <- c(Asamps, which(colnames(cD@data) %in%
                                samplesA[!(samplesA %in% as.character(cD@replicates))]))
    if (!all(samplesA %in% c(colnames(cD@data), as.character(cD@replicates))))
      warning("Some members of 'samplesA' were not found!")
    samplesA <- Asamps
  }
  if (length(samplesA) == 0)
    stop("Can't find any data for sample set A!")
  if (is.character(samplesB)) {
    Bsamps <- which(as.character(cD@replicates) %in% samplesB)
    if (!all(samplesB %in% cD@replicates))
      Bsamps <- c(Bsamps, which(colnames(cD@data) %in%
                                samplesB[!(samplesB %in% as.character(cD@replicates))]))
    if (!all(samplesB %in% c(colnames(cD@data), as.character(cD@replicates))))
      warning("Some members of 'samplesB' were not found!")
    samplesB <- Bsamps
  }
  if (length(samplesB) == 0)
    stop("Can't find any data for sample set B!")
  if (!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  Adata <- cD@data[, samplesA]
  Bdata <- cD@data[, samplesB]
  if (normaliseData) {
    if ("libsizes" %in% names(cD@sampleObservables))
      libsizes <- cD@sampleObservables$libsizes
    else libsizes <- rep(1, ncol(cD))
    Adata <- t(t(Adata)/as.vector(libsizes[samplesA])) *
      mean(libsizes[c(samplesA, samplesB)])
    Bdata <- t(t(Bdata)/as.vector(libsizes[samplesB])) *
      mean(libsizes[c(samplesA, samplesB)])
  }
  if ("seglens" %in% names(cD@rowObservables)) {
    Adata <- Adata/cD@rowObservables$seglens
    Bdata <- Bdata/cD@rowObservables$seglens
  }
  else if ("seglens" %in% names(cD@cellObservables)) {
    Adata <- Adata/cD@cellObservables$seglens
    Bdata <- Bdata/cD@cellObservables$seglens
  }
  Adata <- colSums(t(Adata))/length(samplesA)
  Bdata <- colSums(t(Bdata))/length(samplesB)
  Azeros <- which(Adata == 0)
  Bzeros <- which(Bdata == 0)
  nonzeros <- which(Adata != 0 & Bdata != 0)
  infRatio <- ceiling(max(abs((log2(Adata) - log2(Bdata))[nonzeros]),
                          na.rm = TRUE))
  if (!is.null(scale) && scale > infRatio)
    infRatio <- scale
  M <- log2(Adata) - log2(Bdata)
  M[Azeros] <- -infRatio - 2
  M[Bzeros] <- infRatio + 2
  A <- (log2(Adata) + log2(Bdata))/2
  A[Azeros] <- log2(Bdata[Azeros])
  A[Bzeros] <- log2(Adata[Bzeros])
  plot(y = M, x = A, ylim = c(-infRatio - 3, infRatio + 3),
       axes = FALSE, xlab = xlab, ylab = ylab, ...)
  axis(side = 1)
  maxis <- pretty((-infRatio + 1):(infRatio - 1), min.n = 3,
                  n = length(axTicks(side = 2)))
  maxis <- maxis[maxis < infRatio & maxis > -infRatio]
  axis(side = 2, at = c(-infRatio - 1, maxis, infRatio + 1),
       labels = c(-Inf, maxis, Inf))
  abline(h = c(-1, 1) * (1 + infRatio), col = "orange", lty = 3)
  return(cbind(M,A))
}

## ---- bam2counts
.sebenv$gwcov <- function(bins,mygroupdf,path=".",filecol="File",paired=TRUE,...) {
  require(ShortRead)
  require(reshape2)
  require(GenomicAlignments)
  #tmpdir <- getwd()
  #setwd(path)
  #this functions takes bins (GRanges created by tileGenome) and a data-frame containing metainformation to a list of bam files with the bam-file name as rownames.
  #The bam files are read in and reads are counted for each bin saved in a data.frame in the long and wide format (for ggplot) including the provided annotation in the bin object and mygroupdf
  #bins should ideally have an ID column and a type column
  #e.g.
  # GRanges object with 2 ranges and 3 metadata columns:
  #       seqnames            ranges strand |            ID        type
  #          <Rle>         <IRanges>  <Rle> |   <character> <character>
  #   [1]   I_Chr0 [     1,  500000]      * |      I_Chr0_1      Bin5e5
  #   [2]   I_Chr0 [500001, 1000000]      * | I_Chr0_500001      Bin5e5
  #         genome
  #       <factor>
  #   [1]        I
  #   [2]        I
  #   -------
  #   seqinfo: 30 sequences from an unspecified genome
  #GRanges object with 6 ranges and 2 metadata columns:
  #      seqnames           ranges strand |          ID        type
  #         <Rle>        <IRanges>  <Rle> | <character> <character>
  #  [1]     Chr1 [     1,  50000]      * |      Chr1_1      Bin5e4
  #  [2]     Chr1 [ 50001, 100000]      * |  Chr1_50001      Bin5e4
  #if bins has a genome column, it will be carried over!
  #mygroupdf should looks something like this:
  #   Line Replicate Time  Library       File
  # 1    C         1 wk10 C_1_wk10 file_1.bam
  # 2    C         1  wk4  C_1_wk4 file_2.bam
  myannotcols <- c("seqnames","start","ID","genome")
  bamFls <- mygroupdf[,filecol]
  bin.counts.wide <- matrix(NA, nr=length(bins),nc=length(bamFls))
  colnames(bin.counts.wide) <- mygroupdf[,filecol]
  rownames(bin.counts.wide) <- bins$ID
  message("Paired/Single End: assuming paired=",paired)
  for (bam in bamFls){
    message("reading in ",paste(path,bam,sep="/"))
    if (paired) {
      bin.counts.wide[,bam] <- countOverlaps(bins,readGAlignmentPairs(paste(path,bam,sep="/")),...)
    } else {
      bin.counts.wide[,bam] <- countOverlaps(bins,readGAlignments(paste(path,bam,sep="/")),...)
    }
    #counts_marks[,names(bamFls[i])] <- countOverlaps(gr4, aln[hits==1])
  }
  bin.counts.wide <- cbind(as.data.frame(bin.counts.wide),as.data.frame(bins)[,myannotcols])
  bin.counts.wide$ID  <- as.factor(bin.counts.wide$ID)
  #bin.counts.wide$start <- start(bins)
  #bin.counts.wide$seqnames <- as.character(seqnames(bins))
  #bin.counts.wide$ID <- bins$ID
  #if(!is.null(bins$genome)) {
  #bin.counts.wide$genome <- bins$genome
  #bin.counts.wide$genome[grep(Linepatterns[1],bins$ID)] <- Linecodes[1]
  #bin.counts.wide$genome[grep(Linepatterns[2],bins$ID)] <- Linecodes[2]
  #} else {
  #message("bins object doesnt't have genome column, using NA instead")
  #bin.counts.wide$genome <- rep(NA,length(bins$ID))
  #}
  bin.counts.long <- melt(bin.counts.wide,id.vars=myannotcols)
  #adding mygroupdf information to long formate (note, it is not possible to my knowledge to add it to wide format)
  for (annot in colnames(mygroupdf)){
    bin.counts.long[,annot] <- rep(mygroupdf[,annot],each=length(bins))
  }
  #binannot <- bin.counts.wide[,which(colnames(bin.counts.wide) %in% c("ID","start","seqnames","genome"))]
  return(list(wide=bin.counts.wide,long=bin.counts.long))
}

## ---- bam2counts2
.sebenv$gwcov2 <- function(bins,mygroupdf,path=".",filecol="File") {
  require(ShortRead)
  require(reshape2)
  require(GenomicAlignments)
  myannotcols <- c("seqnames","start","ID","genome")
  bamFls <- mygroupdf[,filecol]
  bin.counts.wide <- matrix(NA, nr=length(bins),nc=length(bamFls))
  colnames(bin.counts.wide) <- mygroupdf[,filecol]
  rownames(bin.counts.wide) <- bins$ID
  for (bam in bamFls){
    message("reading in ",paste(path,bam,sep="/"))
    #     bin.counts.wide[,bam] <- countOverlaps(bins,readGAlignments(paste(path,bam,sep="/")))
    bin.counts.wide[,bam] <- countOverlaps(bins,readGAlignments(paste(path,bam,sep="/")))
    #counts_marks[,names(bamFls[i])] <- countOverlaps(gr4, aln[hits==1])
  }
  bin.counts.wide <- cbind(as.data.frame(bin.counts.wide),as.data.frame(bins)[,myannotcols])
  bin.counts.wide$ID  <- as.factor(bin.counts.wide$ID)
  bin.counts.long <- melt(bin.counts.wide,id.vars=myannotcols)
  #adding mygroupdf information to long formate (note, it is not possible to my knowledge to add it to wide format)
  for (annot in colnames(mygroupdf)){
    bin.counts.long[,annot] <- rep(mygroupdf[,annot],each=length(bins))
  }
  return(list(wide=bin.counts.wide,long=bin.counts.long))
}

gwcov.long <- function(bin.counts.wide,mygroupdf) {
  myannotcols <- c("seqnames","start","ID","genome")
  #	bin.counts.long <- melt(bin.counts.wide,id.vars=myannotcols,measure.vars=mygroupdf[,filecol])
  bin.counts.long <- melt(bin.counts.wide,id.vars=myannotcols)
  #adding mygroupdf information to long formate (note, it is not possible to my knowledge to add it to wide format)
  for (annot in colnames(mygroupdf)){
    bin.counts.long[,annot] <- rep(mygroupdf[,annot],each=nrow(bin.counts.wide))
  }
  return(bin.counts.long)
}
.sebenv$gwcov.long  <- gwcov.long
rm(gwcov.long )

hybrid_normalizing <- function(mycov,mygroupdf,hybridLine="iC",parentLines=c("i","C"),RefGenomes=c("I-Genome","C-Genome")) {
  if (!c("Line") %in% colnames(mygroupdf) ) { stop("mygroupdf needs Line column containing parental and hybrid Line annotation for each library")}
  if ( !all(c(hybridLine,parentLines) %in% levels(mygroupdf$Line)) ) { stop("mygroupdf Line column annotation needs to macht parentLines and hybridLine values.",hybridLine,parentLines," groupdf ",levels(mygroupdf$Line))}
  if ( !all(RefGenomes %in% levels(mycov$genome)) ) { stop("mycov$genome needs to contain RefGenomes")}
  #mycov needs genome column corresponding to RefGenomes parameter values
  #mygroupdf needs Line column, entries have to correspond to parentLines and hybridLine parameters.
  #Refgenome order corresponds to parentLine order! (should ideally be same genomes!)
  #bin_lines <- sapply(strsplit(rownames(mycov),"_"),function(x) x[[1]][1])
  mycovhybrid <- mycov[,which(mygroupdf$Line==hybridLine)]
  libsize.base.hybrid <- mean(colSums(mycovhybrid))
  mycov[,which(mygroupdf$Line==hybridLine)] <- t(t(mycovhybrid * (libsize.base.hybrid)) / colSums(mycovhybrid))
  ##all hybrid libraries are now normalized and can now serve as base levels for the indivdual parents.
  ##The strategy is to determine average libsizes (1 for each parental genome seperatly) by counting normalized hybrid reads mapping for the respective parental genome.
  libsize.base.I <- mean(colSums(mycov[mycov$genome==RefGenomes[1],which(mygroupdf$Line==hybridLine)]))
  libsize.base.C <- mean(colSums(mycov[mycov$genome==RefGenomes[2],which(mygroupdf$Line==hybridLine)]))
  mycovparent <- mycov[,which(mygroupdf$Line==parentLines[1])]
  mycov[,which(mygroupdf$Line==parentLines[1])] <- t(t(mycovparent * (libsize.base.I)) / colSums(mycovparent))
  mycovparent <- mycov[,which(mygroupdf$Line==parentLines[2])]
  mycov[,which(mygroupdf$Line==parentLines[2])] <- t(t(mycovparent * (libsize.base.C)) / colSums(mycovparent))
  return(mycov)
}
.sebenv$hybrid_normalizing <- hybrid_normalizing
rm(hybrid_normalizing)


plotMAgg <- function(x, y=NULL, uselog=TRUE,offset=0,norm=FALSE,method="TMM",dcolor="green", lcolor="red", sp=0.5, ...) {
  if (is.null(y)) {
    if (ncol(x)!=2) {
      stop("if y is omitted, x has to have exactly 2 columns")
    }
    mygroup   <- factor(colnames(x))
    myDGE     <- DGEList(counts=x,group=mygroup)
    ##Note norm factors is only a correction with relation to library size
    # readsCYMlow <-  readsCYM[quantile(rowSums(readsCYM@data), probs = 0.75) > rowSums(readsCYM@data),]
    lib.sizes.effective.tmm <- myDGE$sample$lib.size * calcNormFactors(myDGE,method="TMM")$samples$norm.factors
    lib.sizes.effective.none <- myDGE$sample$lib.size * calcNormFactors(myDGE,method="none")$samples$norm.factors
    lib.sizes.effective.uq <- myDGE$sample$lib.size * calcNormFactors(myDGE,method="upperquartile")$samples$norm.factors
    lib.sizes.effective.uq9 <- myDGE$sample$lib.size * calcNormFactors(myDGE,method="upperquartile",p=0.9)$samples$norm.factors
    tmm <- log2(lib.sizes.effective.tmm[1]/lib.sizes.effective.tmm[2])
    none <- log2(lib.sizes.effective.none[1]/lib.sizes.effective.none[2])
    uq <- log2(lib.sizes.effective.uq[1]/lib.sizes.effective.uq[2])
    uq9 <- log2(lib.sizes.effective.uq9[1]/lib.sizes.effective.uq9[2])
    mydesign <- model.matrix(~mygroup)
    #      plotSmear(myDGEnone)
    y  <- x[,2]
    x  <- x[,1]
  }
  if (uselog) {
    x <- log2(x+offset)
    y <- log2(y+offset)
    message("perfoming log2 on data. Adding offset: ",offset)
  }
  Mean <- (x + y)/2
  Diff <- x - y
  MA.df <- data.frame(Mean,Diff)
  colnames(MA.df) <- c("A","M")
  MA.plot <- ggplot(data=MA.df,aes(x=A,y=M),na.rm=TRUE) +
    geom_point(alpha=.15, color="black", size=1.3) +
    geom_smooth(span=sp,color="red", se = FALSE) +
    labs(subtitle=paste(colnames(x),collapse="\n")) +
    geom_hline(yintercept=tmm,color="blue")+ geom_text(aes( 0, tmm, label = "TMM", size = 3))+
    geom_hline(yintercept=none)+ geom_text(aes( 0, none, label = "none", size = 3)) +
    geom_hline(yintercept=uq,color="orange")+ geom_text(aes( 0, uq, label = "upperquartile", size = 3)) +
    geom_hline(yintercept=uq9,color="green")+ geom_text(aes( 0, uq9, label = "upperquartile p=0.9", size = 3))
  return(MA.plot)
}
.sebenv$plotMAgg <- plotMAgg
rm(plotMAgg)

.sebenv$panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, same.font.size=FALSE,...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, ...)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  if (same.font.size){
    text(0.5, 0.5, txt, cex = cex.cor )
  } else {
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  }
}
#pairs(iris,upper.panel=function(x,y) panel.cor(x,y,method="spearman"))


###
.sebenv$annotate_intervals_arabidopsis <- function(interval) {
  interval$overlapgenes <- rep("none", length(interval))
  ref <- mRNA
  overlap <- findOverlaps(ref,interval)
  uniqueref <- !duplicated(queryHits(overlap))
  interval$overlapgenes[subjectHits(overlap)[uniqueref]] <- as.character(ref[queryHits(overlap)[uniqueref],]$Name)
  #TAIR TEs vs complete set of interval
  interval$overlapTEs <- rep("none", length(interval))
  interval$overlapTEsFamily <- rep("none", length(interval))
  interval$overlapTEsSuperfamily <- rep("none", length(interval))
  ref <- gffTEs
  overlap <- findOverlaps(ref,interval)
  uniqueref <- !duplicated(queryHits(overlap))
  interval$overlapTEs[subjectHits(overlap)[uniqueref]] <- as.character(ref[queryHits(overlap)[uniqueref],]$Name)
  interval$overlapTEsFamily[subjectHits(overlap)[uniqueref]] <- as.character(ref[queryHits(overlap)[uniqueref],]$Family)
  interval$overlapTEsSuperfamily[subjectHits(overlap)[uniqueref]] <- as.character(ref[queryHits(overlap)[uniqueref],]$Super_Family)
  ref <- sidRNAs
  interval$overlapSidRNAs <- rep("none", length(interval))
  overlap <- findOverlaps(ref,interval)
  uniqueref <- !duplicated(queryHits(overlap))
  interval$overlapSidRNAs[subjectHits(overlap)[uniqueref]] <- as.character(ref[queryHits(overlap)[uniqueref],]$Name)
  return(interval)
}

# Input: 2 GRange Objects %>%
# Output: 1 GRange Object (same as first input with additional Column)
# Parameter
#  naming: Name of the new Column
#  uniqueSubject: if subjects overlaps muliple times, should only first be reported? FALSE most sensible in most cases..
.sebenv$overlap_annot <- function(query, subject, naming, ColName="Name", uniqueSubject = TRUE) {
  query@elementMetadata[[naming]] <- rep("none", length(query))
  overlap <- findOverlaps(subject,query)
  if (uniqueSubject) {
    uniquesubject <- !duplicated(queryHits(overlap))
    query@elementMetadata[[naming]][subjectHits(overlap)[uniquesubject]] <-
      as.character(subject[queryHits(overlap)[uniquesubject],]@elementMetadata[[ColName]])
  } else {
    query@elementMetadata[[naming]][subjectHits(overlap)] <-
      as.character(subject[queryHits(overlap),]@elementMetadata[[ColName]])
  }
  return(query)
}
# naming="test"
# ColName="score2"
# df <- data.frame(
#     seqnames = rep(c("chr1", "chr1", "chr1", "chr1"), c(1, 3, 2, 4)),
#     start = c(10, 12, 12, 132, 134, 152, 153, 160, 166, 260),
#     end =   c(14, 15, 17, 132, 155, 154, 159, 166, 171, 290),
#     strand = rep(strand(c("*", "*", "*", "*", "*")), c(1, 2, 2, 3, 2)),
#     row.names = head(letters, 10))
# query <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
# df2 <- data.frame(
#     seqnames = rep(c("chr1", "chr1", "chr1", "chr1"), c(1)),
#     start = c(1, 105, 125,290),
#     end =   c(104, 220, 133, 332),
#     strand = rep(strand(c("*", "*", "*", "*")), c(1)),
#     score2 = 11:14,
#     row.names = head(letters, 4))
# subject <- makeGRangesFromDataFrame(df2, keep.extra.columns=TRUE)
# overlap_annot(query, subject, "test", "score2")
# overlap_annot(query, subject, "test", "score2", FALSE)


.sebenv$headtail  <-  function(df, num = 3) {
  nr <- nrow(df)
  df[c(1:num,(nr-num+1):nr),]
}

.sebenv$venn.pretty <- function (vennlist,filename="venn.pdf",cat.pos=1,cat.cex=1,...) {
  # sample five-set Venn Diagram
  require(VennDiagram)
  require(gplots)
  plotObject <- venn.diagram(
                             x        = vennlist,
                             filename = NULL,
                             col      = "black",
                             fill     = c("dodgerblue", "goldenrod1", "seagreen3", "orchid3", "darkorange1")[1:length(vennlist)],
                             cat.col     = c("dodgerblue", "goldenrod1", "seagreen3", "orchid3", "darkorange1")[1:length(vennlist)],
                             alpha    = 0.50,
                             #                 cex      = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                             #                              1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                             #                 cat.col = c("dodgerblue", "goldenrod0", "darkorange1", "seagreen3", "orchid3")[1:length(vennlist)],
                             cat.cex = cat.cex,,
                             cat.fontface = "bold",
                             margin = 0.05,
                             cat.pos = cat.pos,
                             ext.text = TRUE,
                             ...
                             )
  ggsave(filename=filename,plot=plotObject)
  return(venn(vennlist,show.plot=FALSE))
}



.sebenv$fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
      x <- as.vector(match(elements, x))
      }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

.sebenv$overlapGroups <- function (listInput, sort = TRUE) {
  # listInput could look like this:
  # $one
  # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  # $two
  # [1] "a" "b" "d" "e" "j"
  # $three
  # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
  listInputmat    <- fromList(listInput) == 1
  #     one   two three
  # a  TRUE  TRUE  TRUE
  # b  TRUE  TRUE FALSE
  #...
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
    # attr(,"groups")
    #   one   two three
    # FALSE FALSE  TRUE
    #  f  i
    # 12 13
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

.sebenv$dfinit <- function (col, row, init = NA) {
  # input col, row = vectors with column/row names as stings
  df <- as.data.frame(matrix(init, ncol = length(col), nrow = length(row)))
  rownames(df) <- row
  colnames(df) <- col
  return(df)
}


.sebenv$checkPacks<-function(path){
  require(NCmisc)
  require(stringr)
  require(dplyr)

  ## get all R files in your directory
  ## by the way, extract R code from Rmd: http://felixfan.github.io/extract-r-code/
  files<-list.files(path)[str_detect(list.files(path), file.path(".R$"))]

  ## extract all functions and which package they are from 
  ## using NCmisc::list.functions.in.file
  funs<-unlist(lapply(paste0(path, "/", files), list.functions.in.file))
  packs<-funs %>% names()

  ## "character" functions such as reactive objects in Shiny
  characters<-packs[str_detect(packs, "^character")]

  ## user defined functions in the global environment
  globals<-packs[str_detect(packs, "^.GlobalEnv")]

  ## functions that are in multiple packages' namespaces 
  multipackages<-packs[str_detect(packs, ", ")]

  ## get just the unique package names from multipackages
  mpackages<-multipackages %>%
    str_extract_all(., "[a-zA-Z0-9]+") %>%
    unlist() %>%
    unique()
  mpackages<-mpackages[!mpackages %in% c("c", "package")]

  ## functions that are from single packages
  packages<-packs[str_detect(packs, "package:") & !packs %in% multipackages] %>%
    str_replace(., "[0-9]+$", "") %>%
    str_replace(., "package:", "") 

  ## unique packages
  packages_u<-packages %>%
    unique() %>%
    union(., mpackages)

  return(list(packs=packages_u, tb=table(packages)))

}
# dfinit(c("col1", "col2", "col3"), c("r1", "r2"))
#     a  b
# r1 NA NA
# r2 NA NA

# computes average methylation for each sample in a "methylBase" object:
# by Sebastian Mueller
# e.g. meth:
# methylBase object with 6 rows
# --------------
#   chr start end strand coverage1 numCs1 numTs1 coverage2 numCs2 numTs2 coverage3 numCs3 ...
# 1   1   110 110      *         3      3      0        26     25      1         2      2
# 2   1   115 115      *         4      4      0        50     49      1         3      3
#
# usage: meth_summary(meth)
# output:
# sample1   sample2
# 0.7661584 0.7459081
.sebenv$meth_summary <- function(mymeth) {
  require(dplyr)
  mynames <- mymeth@sample.ids
  mysummary <- as.data.frame(mymeth) %>%
    select(starts_with("num")) %>%
    colSums(na.rm = TRUE) %>%
    tapply(rep(seq_along(mynames), each = 2), function(x) x[1] / sum(x))
  names(mysummary) <- mynames
  return(mysummary)
}

.sebenv$GOfun <- function (targetset, background) {
  interestingGenes <- factor(as.integer( names(background) %in% ((targetset)) ) )
  names(interestingGenes) <- names(background)
  res <- NA
  for (ont in c("BP","MF","CC")) {
    # GOdata            <- new("topGOdata", ontology = ont, allGenes = interestingGenes, annot = annFUN.gene2GO, gene2GO = background)
    GOdata            <- new("topGOdata",
                             ontology = ont,
                             allGenes = interestingGenes,
                             annot = annFUN.gene2GO,
                             # named list of character vectors.  The list names are genes identifiers.  For eachgene the character vector contains the GO identifiers it maps to.  Only the mostspecific annotations are required
                             gene2GO = background)
    resultFisher      <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    resultFisher.Elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    enrichRes         <- GenTable(GOdata, classicFisher = resultFisher, elimFisher = resultFisher.Elim, orderBy = "classicFisher", ranksOf ="elimFisher", topNodes = 20)
    res               <- rbind(res,ont,enrichRes)
  }
  return(res)
}
.sebenv$table2 <- function(...) table(..., useNA = "always")  # counts all three

.sebenv$table2d <- function(...) ftable(addmargins(table(..., useNA = "always")))

.sebenv$table_prop <- function(...) round(100*prop.table(table(...)),1)


# source("/home/sm934/code/R-code-misc-seb/R-functions-seb.r")
attach(.sebenv)

