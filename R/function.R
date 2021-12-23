library(ggplot2)
library(Biostrings)
library(stringdist)
library(rlist)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(parallel)
library(pheatmap)
library(reshape2)
library(ggtree)
library(data.tree)
library(networkD3)
library(stringr)
library(rlist)
library(gtable)
library(tidyverse)
library(ape)
library(ggnewscale)
library(purrr)
library(IRanges)
library(cowplot)
library(S4Vectors)

#' Title
#' @title ReadFasta
#' @description  Function to read fasta file to DNAString object
#' @param filename The input fasta file name
#'
#' @return A DNAString object
#' @export
#' @importFrom Biostrings DNAString
#'
#' @examples
#' fafile=paste0(system.file("extdata",package = 'LinTInd'),"/V3.fasta")
#' ReadFasta(fafile)
#'
ReadFasta = function(filename){
  sv = read.table(filename)
  scarfull = DNAString(sv[2,1])
  return(c("scarfull" = scarfull))
}

#' Title
#' @title ReadCutsite
#' @description Function to create a reference dataframe include each position and its' group
#' @usage ReadCutsite(cutsite,reftype="Accurate")
#' @param segref The cutsite file
#' @param reftype Choose the reference type you want, if reftype="Accurate" (default), there will only the target sites be generated; if reftype="All", each site will be generated
#'
#' @return reference dataframe
#' @export
#'
#' @examples
#' data("example_data",package = "LinTInd")
#' ReadCutsite(cutsite)
#' ReadCutsite(cutsite,reftype="All")
#'
ReadCutsite = function(segref,reftype=NULL){
  colnames(segref) = c("indx","start","end")
  scar = NULL
  type = NULL
  if(is.null(reftype) | reftype=="Accurate"){
    for (i in 2:nrow(segref)) {
      scar = c(scar,segref[i,]$start:segref[i,]$end)
      type = c(type,rep(segref[i,]$indx,(segref[i,]$end-segref[i,]$start)+1))
    }
    scarseg = data.frame("scar" = scar,"type" = as.character(type))
  }else{
    endsite<-NA
    for (i in 2:nrow(segref)) {
      endsite<-c(endsite,segref[["end"]][i]+segref[["start"]][1])
    }
    endsite[nrow(segref)]<-segref[["end"]][1]
    scarseg = data.frame("scar" = c(1:endsite[nrow(segref)]),"type" = NA)
    #endsite<-is.na(endsite)
    for (i in rev(endsite)) {
      if(!is.na(i)){
        scarseg$type[1:i]<-which(endsite==i)-1
      }else{
        break
      }
    }
  }
  return(scarseg)
}

#' Title
#' @title align_to_range
#'
#' @param p A base sequence in character format
#' @param s A base sequence in character format
#' @param cut The distance between the starting sites of two fragments
#'
#' @return A list include two IRanges instances (deletion and insertion)
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom Biostrings DNAString pairwiseAlignment nucleotideSubstitutionMatrix
#' @export
#' @examples
#' align_to_range(p="AAGG---AAATTTCGGAATAAGGAATTT",s="AAGGCCCCAAATTT-CGGAATAAGGAATTT",cut=0)
#'
#'
align_to_range = function(p,s,cut){
  pn <- strsplit(p,"")[[1]]
  sn <- strsplit(s,"")[[1]]
  lenp <- length(pn)
  index <- 1
  i <- 0
  del_flag <- FALSE
  in_flag <- FALSE
  del_start<-c()
  del_end<-c()
  ins_start<-c()
  ins_width<-c()
  ins_editseq<-c()
  while(i < lenp){
    i <- i + 1
    if(sn[[i]] == '-'){
      if(!del_flag){
        del_flag <- TRUE
        del_start<-c(del_start,index)
        width <- 1
      }
      else{
        width <- width + 1
      }
    }
    else{
      if(del_flag){
        del_flag <- FALSE
        del_end<-c(del_end,index)
        #print(paste("del stop width", width))
      }
    }
    if(pn[[i]] == '-'){
      if(!in_flag){
        in_flag <- TRUE
        ins_start<-c(ins_start,index)
        width <- 1

        editseq<-sn[[i]]
      }
      else{
        width <- width + 1

        editseq<-paste0(editseq,sn[[i]])
      }
    }
    else{
      if(in_flag){
        in_flag <- FALSE
        ins_width<-c(ins_width,width)
        ins_editseq<-c(ins_editseq,editseq)
        #print(paste("in stop width", width))
      }
    }
    if(pn[[i]] != '-')
      index <- index + 1
  }
  if(del_flag){
    del_flag <- FALSE
    del_end<-c(del_end,index)
  }
  if(in_flag){
    in_flag <- FALSE
    ins_width<-c(ins_width,width)
    ins_editseq<-c(ins_editseq,editseq)
    #print(paste("in stop width", width))
  }

  ins_start = ins_start-cut
  ins_end = ins_start + ins_width
  del_start = del_start-cut
  del_end = del_end-cut

  ins<-IRanges(start = ins_start,end = ins_end)
  mcols(ins)$seq <- ins_editseq
  del<-IRanges(start = del_start,end = del_end)


  return(list("del" = del,"ins" = ins))

}

#' Title
#' @title FindIndel
#' @description This function can ident indels for each reads in input data, and create IRanges instances for deletion and insertion
#' @param data data frame, include cell barcode, UMI and reads.
#' @param scarfull DNAString of reference sequence
#' @param scar The cutsite data frame
#' @param align_score The minimum alignment score that matched sequence should get, default in this parameter is the score that the reads which all of the target set were cutted got
#' @param type Group name for this data ("None" in default)
#' @param indel.coverage Choose indels selected scope: "Accurate" (default) means only the indles happenned in target site will be idented; "All" means each indel will be detected even they locate on the anchors
#' @param cln The number of threads
#'
#' @return list include IRanges instances (deletion and insertion), a data frame of reads' informations, reference sequenc, dataframe of cut sites
#' @export
#' @importFrom Biostrings DNAString pairwiseAlignment nucleotideSubstitutionMatrix subseq matchPattern pattern alignedPattern alignedSubject
#' @importFrom IRanges IRanges subject
#' @importFrom parallel makeCluster clusterEvalQ parLapply stopCluster clusterExport
#' @importFrom BiocGenerics score
#' @examples
#' data("example_data",package = "LinTInd")
#' scarinfo<-FindIndel(data=data,scarfull=ref,scar=cutsite,indel.coverage="All",type="test",cln=8)
#'
#'
FindIndel = function(data,scarfull,scar,align_score=NULL,type=NULL,indel.coverage=NULL,cln){
  scarref<-ReadCutsite(scar,reftype=indel.coverage)
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  if(is.null(align_score)){
    align_score = length(scarfull)- 2*(scar["end"][1,] - scar["start"][1,])-6
  }else{
    align_score = align_score
  }
  if(is.null(type)){
    type="None"
  }else{
    type=type
  }
  #type="none"
  find_barcode<-function(data){
    #tycpe="none"
    s3<-DNAString(as.character(data))
    alig<-pairwiseAlignment(scarfull,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    if(any(is.null(indel.coverage),indel.coverage=="Accurate")){
      scarshort = subseq(as.character(scarfull[[1]]),scar["start"][1,],scar["end"][1,])
    }else{
      scarshort = as.character(scarfull[[1]])
    }
    if(score(alig)<=align_score){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos = matchPattern(scarshort,as.character(pattern(alig)),with.indels = TRUE,max.mismatch = 10)
      r_read = as.character(subject(alig))
      if(length(scar_pos)!= 0){
        r_scar = testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
      }else{
        r_scar = "unknown"
      }
      stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)

      p <- as.character(alignedPattern(alig)[[1L]])
      s <- as.character(alignedSubject(alig)[[1L]])
      if(any(is.null(indel.coverage),indel.coverage=="Accurate")){
        delins = align_to_range(p,s,scar["start"][1,])
      }else{
        delins = align_to_range(p,s,0)
      }
      del = delins$del
      ins = delins$ins
      if(TRUE %in% (del@start<0) | TRUE %in% (ins@start<0)){
        r_scar = "unknown"
      }
    }
    fin_dat<-data.frame(new.reads=r_read, scar.BC=r_scar,type=type)
    return(list("del" = del,"ins" = ins,fin_dat))
  }

  cl = makeCluster(cln,setup_strategy = "sequential")
  clusterEvalQ(cl,library(Biostrings))
  #environment(align_score) <- .GlobalEnv
  #environment(scarfull) <- .GlobalEnv
  #environment(scar) <- .GlobalEnv
  #environment(data) <- .GlobalEnv
  #environment(mat) <- .GlobalEnv
  #environment(type) <- .GlobalEnv
  #environment(find_barcode) <- .GlobalEnv
  #environment(testFunction) <- .GlobalEnv
  #environment(indel.coverage) <- .GlobalEnv
  clusterExport(cl,c('align_score','scarfull','scar','data','mat','indel.coverage','find_barcode','testFunction',"align_to_range","type"),envir = environment())
  scar_BC = parLapply(cl,data$Read.Seq,find_barcode)
  stopCluster(cl)

  #output
  data_2<-do.call("rbind",sapply(scar_BC,function(x){return(x[3])}))
  data<-cbind(data,data_2)
  scar_BC<-scar_BC[data$scar.BC!="unknown"]
  data<-data[data$scar.BC!="unknown",]
  saveRDS(list(scar_BC,data),"reads_metadata.rds")
  write.table(data,"all_UMI_reads_scar_full.txt",quote=FALSE,sep="\t",row.names=FALSE)
  data$scar_f<-gsub("[-]", "",as.character(data$scar.BC))
  data_v1<-data
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"UMI_reads_scar_full.txt",quote=FALSE,sep="\t",row.names=FALSE)
  saveRDS(scar_BC,"indel.rds")
  return(list("INDEL" = scar_BC,"Scar" = data_v1,"indel.coverage"=indel.coverage,"cutsite"=scar,"ref"=scarfull,"scarref"=scarref))
}


#' Title
#' @title change_form_stat
#' @param indel List include two IRanges instances, contain start and end site of deletions and inserstions
#'
#' @return A scar form string
#'
#' @export
#' @importFrom IRanges IRanges
#' @examples
#' data("example_data",package = "LinTInd")
#' change_form_stat(cellsinfo$indel[[1]])
#'
#'
change_form_stat<-function(indel){
  indel<-indel[c(1,2)]
  indel<-unlist(indel)
  if(length(indel)==0){
    return("unknown")
  }else{
    ins<-data.frame(indel[2])
    #ins$seq=as.character(indel[[2]]@elementMetadata$seq)
    site_ins<-apply(ins,1,function(x){c(x[1]:x[2])})
    if(dim(ins)[1]==1){
      site_ins<-list(as.numeric(site_ins))
    }
    cutsite_ins<-lapply(site_ins,function(x){unique(scarref$type[scarref$scar %in% x])})
    tag_ins<-apply(ins,1,function(x){paste0(x[3],"I+",x[1],x[4])})
    tag_ins<-lapply(tag_ins,function(x){rep(x,length(cutsite_ins[[which(tag_ins==x)]]))})
    tag_ins<-unlist(tag_ins)
    cutsite_ins<-unlist(cutsite_ins)
    del<-data.frame(indel[1])
    site_del<-apply(del,1,function(x){c(x[1]:x[2])})
    if(dim(del)[1]==1){
      site_del<-list(as.numeric( site_del))
    }
    cutsite_del<-lapply(site_del,function(x){unique(scarref$type[scarref$scar %in% x])})
    tag_del<-apply(del,1,function(x){paste0((x[3]-1),"D+",x[1])})
    tag_del<-lapply(tag_del,function(x){rep(x,length(cutsite_del[[which(tag_del==x)]]))})
    tag_del<-unlist(tag_del)
    cutsite_del<-unlist(cutsite_del)
    tag<-c(tag_del,tag_ins)
    cutsite<-c(cutsite_del,cutsite_ins)
    tag_all<-rep("NONE",length(unique(scarref$type)))
    if(length(tag)==0){
      return(paste(tag_all,collapse = "_"))
    }else{
      for(x in c(1:length(tag))){
        if(tag_all[as.numeric(cutsite[x])]=="NONE"){
          tag_all[as.numeric(cutsite[x])]<-tag[x]
        }else{
          tag_all[as.numeric(cutsite[x])]<-paste(tag_all[as.numeric(cutsite[x])],tag[x],sep="&")
        }
      }
    }
    return(paste(tag_all,collapse = "_"))
  }
}
#' Title
#' @title IndelForm
#' @description Generate scar form strings from scarinfo list for each reads
#' @param scarinfo List generate from FindIndel, for more see \code{\link[LinTInd]{FindIndel}}
#' @param cln The number of threads
#'
#' @return A new list of scarinfo, the scarform of each reads will add in the data frame of reads' informations
#' @export
#' @importFrom IRanges IRanges
#' @importFrom parallel makeCluster clusterEvalQ parLapply stopCluster clusterExport
#' @examples
#' data("example_data",package = "LinTInd")
#' IndelForm(scarinfo,cln=4)
#'
IndelForm = function(scarinfo,cln){
  scarref<-scarinfo$scarref
  cl = makeCluster(cln,setup_strategy = "sequential")
  clusterExport(cl,c('scarinfo','change_form_stat',"scarref"), envir = environment())
  scar_form_p<-parLapply(cl,scarinfo$INDEL,change_form_stat)
  stopCluster(cl)
  scar_form<-unlist(scar_form_p)
  scar_form<-gsub(" ","",scar_form)
  data<-scarinfo$Scar
  data$scar_form<-scar_form
  write.csv(data,"indel_pattern.csv",quote=FALSE,row.names = FALSE)
  return(list("INDEL" = scarinfo$INDEL,"Scar" = data,"indel.coverage"=scarinfo$indel.coverage,"cutsite"=scarinfo$cutsite,"ref"=scarinfo$ref,"scarref"=scarinfo$scarref))
}


#' Title
#' @title IndelIdents
#' @description Function to define a scarform for each cell(single cell) or each reads(bulk seq, generate 'cell barcode' for each reads)
#' @param scarinfo List generate from IndelForm, for more see \code{\link[LinTInd]{IndelForm}}
#' @param method.use Select how to determine a scar form string for each cell:
#'                   "reads.num" (default):find the scar with the most reads in the cell；
#'                   "umi.num":find the scar with the most UMIs in the cell；
#'                   "consensus":find the consistent sequences in each cell, and then generate scar form strings from the new reads
#'
#' @param cln The number of threads
#'
#' @return The list generate from FindIndel, but in 'Scar' element a new column contain scar form strings
#' @export
#' @importFrom stringdist stringdistmatrix
#' @importFrom Biostrings consensusString pairwiseAlignment DNAString nucleotideSubstitutionMatrix alignedPattern alignedSubject
#' @importFrom IRanges IRanges
#'
#' @examples
#' data("example_data",package = "LinTInd")
#' IndelIdents(scarinfo,method.use="umi.num",cln=4)
#'
#'
IndelIdents = function(scarinfo,method.use=NULL,cln){
  scarfull<-scarinfo$ref
  scar<-scarinfo$cutsite
  data<-scarinfo$Scar
  scarref<-scarinfo$scarref
  indel.coverage<-scarinfo$indel.coverage
  Cell.BC<-data.frame(table(data$Cell.BC))
  Cell.BC<-Cell.BC[Cell.BC$Freq>1,]
  data_1<-data[data$Cell.BC %in% Cell.BC$Var1,]
  if(all("Cell.BC" %in% names(data),"UMI" %in% names(data))){
    data_1<-data_1[,c("Cell.BC","UMI","scar_f","scar_form")]
  }else if(any("Cell.BC" %in% names(data),"UMI" %in% names(data))){
    names(data_1)[which(names(data_1) %in% c("Cell.BC","UMI"))]<-"Cell.BC"
    data_1<-data_1[,c("Cell.BC","scar_f","scar_form")]
  }else{
    data_1<-data_1[,c("scar_f","scar_form")]
  }
  if(length(names(data_1))==2){
    pattern<-data_1$scar_form
    data_con = data.frame("Cell.BC" = c(1:length(pattern)),
                          "pattern" = pattern,
                          "reads_num" = 1,
                          "reads_pro" = 1,
                          "umi_num" = 1,
                          "umi_pro" = 1,
                          stringsAsFactors = FALSE)
    INDEL_ranges <- scarinfo$INDEL
    INDEL_ranges <- INDEL_ranges[data_con$pattern!="unkown"]
    #INDEL_ranges <-INDEL_ranges[data$Cell.BC %in% Cell.BC$Var1]
    INDEL_ranges_man<-list()
    for(index in seq(length(INDEL_ranges))){
      INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
    }

    data_con<-data_con[data_con$pattern!="unkown",]
    write.csv(data_con,"final_scarform.csv",quote=FALSE,row.names = FALSE)
  }else{
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  # x = as.character(Cell.BC$Var1)[2]c
  # dat =data_1
  max_reads_stat = function(x,dat,method.use=NULL){
    temreads = dat[dat$Cell.BC==x,]
    read_data = data.frame(table(as.character(temreads$scar_f)))
    read_data = read_data[order(-read_data$Freq),]
    scar_data = data.frame(table(as.character(temreads$scar_form)))
    scar_data = scar_data[order(-scar_data$Freq),]
    if(scar_data$Freq[1]==1){
      #can't define
      pattern="unkown"
      reads_num=0
      reads_pro=0
      umi_num=0
      umi_pro=0
      #UMI=0
      del=NA
      ins=NA
    }else{
      if(method.use=="consensus"){
        scarstrdist = stringdistmatrix(as.character(read_data$Var1),as.character(read_data$Var1))
        scarindex = which(apply(scarstrdist,1,function(x){sum(read_data$Freq[which(x<9)])>(sum(read_data$Freq)/3)}))
        if(length(scarindex) > 0){
          fin_read = consensusString(temreads$scar_f[temreads$scar_f %in% as.character(read_data$Var1)[scarindex]])
          reads_pro = round(length(which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex]))/sum(read_data$Freq),4)
          reads_num = length(which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex]))
          umi_num = length(unique(temreads$UMI[which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex])]))
          umi_pro = round(umi_num/length(unique(temreads$UMI)),4)
          fin_read_cons = gsub("\\?","",fin_read)
          s1 = DNAString(fin_read_cons)
          aligc = pairwiseAlignment(scarfull,s1,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension = 1)
          stopifnot(is(aligc, "PairwiseAlignments"), length(x) == 1L)
          p <- as.character(alignedPattern(aligc)[[1L]])
          s.cons <- as.character(alignedSubject(aligc)[[1L]])
          if(any(is.null(indel.coverage),indel.coverage=="Accurate")){
            indel = align_to_range(p,s.cons,scar["start"][1,])
          }else{
            indel = align_to_range(p,s.cons,0)
          }
          pattern =  change_form_stat(indel)
        }else{
          pattern = "unknown"
          reads_num=0
          reads_pro=0
          umi_num=0
          umi_pro=0
          #UMI=0
          #del=NA
          #ins=NA
          indel=NULL
        }
        #print("consensus")
      }else if(method.use=="umi.num"){
        #method.use=="umi.num"
        read_data_umi = unique(temreads)
        read_data_umi = data.frame(table(read_data_umi$scar_form))
        read_data_umi = read_data_umi[order(-read_data_umi$Freq),]
        pattern = as.character(read_data_umi$Var1[1])
        umi_num = read_data_umi$Freq[1]
        umi_pro = round(umi_num/length(unique(temreads)$UMI),4)
        reads_num = scar_data[which(scar_data$Var1 == pattern),"Freq"]
        reads_pro = round(reads_num/sum(scar_data$Freq),4)
        #print("umi.num")
      }else{
        pattern = as.character(scar_data$Var1[1])
        reads_num = scar_data$Freq[1]
        reads_pro = round(reads_num/sum(scar_data$Freq),4)
        umi_num = length(unique(temreads$UMI[which(temreads$scar_form %in% pattern)]))
        umi_pro = round(umi_num/length(unique(temreads$UMI)),4)
        #print("reads.num")
      }
    }
    fin_line = data.frame("Cell.BC" = x,
                          "pattern" = pattern,
                          "reads_num" = reads_num,
                          "reads_pro" = reads_pro,
                          "umi_num" = umi_num,
                          "umi_pro" = umi_pro,
                          stringsAsFactors = FALSE)
    if(method.use=="consensus"){
      return(list(indel,fin_line))
    }else{
      return(fin_line)
    }
  }
  cl = makeCluster(cln,setup_strategy = "sequential")
  clusterEvalQ(cl,library(Biostrings))
  clusterEvalQ(cl,library(stringdist))
  #environment(data_1) <- .GlobalEnv
  #environment(max_reads_stat) <- .GlobalEnv
  #environment(Cell.BC) <- .GlobalEnv
  #environment(scarref) <- .GlobalEnv
  #environment(align_to_range) <- .GlobalEnv
  #environment(change_form_stat) <- .GlobalEnv
  #environment(scarfull) <- .GlobalEnv
  #environment(method.use) <- .GlobalEnv
  #environment(indel.coverage) <- .GlobalEnv
  clusterExport(cl,c('data_1','indel.coverage','method.use','max_reads_stat','Cell.BC',"align_to_range","scarref","change_form_stat","scarfull","mat","scar"), envir = environment())
  data_con<-parLapply(cl,as.character(Cell.BC$Var1),function(x)tryCatch(max_reads_stat(x,dat=data_1,method.use=method.use),error=function(e) NULL))
  stopCluster(cl)
  data_con<-data_con[!sapply(data_con,is.null)]

  if(method.use=="consensus"){
    INDEL_ranges_man <-list()
    for(i in seq(length(data_con))){
      INDEL_ranges_man <- c(INDEL_ranges_man,list(data_con[[i]][[1]]))
    }
    data_con_sub<-data.frame(Cell.BC=NA,pattern=NA,reads_num=NA, reads_pro=NA, umi_num=NA, umi_pro=NA)
    data_con_sub<-data_con_sub[-1,]
    for(i in seq(length(data_con))){
      data_con_sub<-rbind(data_con_sub,data_con[[i]][[2]])
    }
    data_con<-data_con_sub
    INDEL_ranges_man<-INDEL_ranges_man[data_con$pattern!="unkown"]
    data_con<-data_con[data_con$pattern!="unkown",]
    write.csv(data_con,"final_scarform.csv",quote=FALSE,row.names = FALSE)
    saveRDS(INDEL_ranges_man,"indel_ident.rds")
  }else{
    data_con<-do.call("rbind",data_con)
    data_con<-data_con[data_con$pattern!="unkown",]
    write.csv(data_con,"final_scarform.csv",quote=FALSE,row.names = FALSE)

    INDEL_ranges <- scarinfo$INDEL
    #INDEL_ranges <-INDEL_ranges[data$Cell.BC %in% Cell.BC$Var1]
    INDEL_ranges_man<-list()
    for(scarform in data_con$pattern){
      index=which(data$scar_form==scarform)[1]
      INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
    }
    saveRDS(INDEL_ranges_man,"indel_ident.rds")
  }
  }
  return(list("indel"=INDEL_ranges_man,"info"=data_con,"indel.coverage"=indel.coverage,"cutsite"=scarinfo$cutsite,"ref"=scarinfo$ref,"scarref"=scarinfo$scarref))
}


#' Title
#' @title TagProcess
#' @description Split each indel from scar form strings and map indel information to cell barcodes
#' @param data List generate from IndelIdents, for more see \code{\link[LinTInd]{IndelIdents}}
#' @param Cells (optional) Dataframe of cells' annotation, with two columns: "Cell.BC" and "Cell.type"
#' @param prefix (optional) Indels' prefix
#'
#' @return List with two dataframes: Indels for each cell barcode and cells' annotation
#' @export
#' @importFrom S4Vectors na.omit
#'
#' @examples
#' data("example_data",package = "LinTInd")
#' TagProcess(cellsinfo$info,Cells=celltype)
#'
TagProcess = function(data,Cells=NULL,prefix=NULL){

  TagStat = function(x) {
    x = as.character(x)
    umi = x[5]
    CB = x[1]
    x = x[2]
    x = unlist(strsplit(x,"_|&"))
    x = x[!x%in%c("NONE")]
    x = unique(x)
    if(length(x) == 0){
      return(NA)
    }else{
      return(data.frame(Cell.BC = CB,Reads_num = umi,Tag = x))
    }
  }

  #tag = NULL
  #common.CB = NULL
  #if((dim(data)[1]>1)){
  # for (i in 1:(dim(data)[1]-1)) {
  #  common.CB = c(common.CB,intersect(data$Cell.BC[i], data$Cell.BC[i+1]))
  #}
  #}else{
  #common.CB = data$Cell.BC[1]
  #}

  #for (i in 1:dim(data)[1]) {
  if(!is.null(Cells)){
    data=data[data$Cell.BC %in% Cells$Cell.BC,]
    data$Cell.type<-Cells$Cell.type[match(data$Cell.BC,Cells$Cell.BC)]
  }
  #data[[i]]=data[[i]][data[[i]]$Cell.BC %in% common.CB,]
  tagi = apply(data,1,TagStat)
  tagi = do.call("rbind",tagi)
  tagi = na.omit(tagi)
  #tagi = data.frame(table(tagi$Tag)/length(as.character(unique(tagi$Cell.BC))))
  #black list filter
  if(!is.null(Cells)){
    tagi$Cell.type<-Cells$Cell.type[match(tagi$Cell.BC,Cells$Cell.BC)]
  }
  if(!is.null(prefix)){
    tagi$Tag = paste(prefix[i], tagi$Tag, sep = "")
  }else{
    tagi$Tag=tagi$Tag
  }
  #tag = rbind(tag,tagi)
  #}
  return(list("tag"=tagi,"celltype"=Cells))
}


#' Title
#' @title BuildTree
#' @description Generate an array generant tree of a data.tree data structure and save it
#' @param tag List generate from TagProcess, for more see \code{\link[LinTInd]{TagProcess}}
#'
#' @return list with two elements, a data.tree data structure and a dataframe of array information for each cell barcode
#' @export
#' @importFrom networkD3 saveNetwork diagonalNetwork
#' @importFrom data.tree ToListExplicit as.Node
#' @examples
#' data("example_data",package = "LinTInd")
#' treeinfo<-BuildTree(tag)
#'
BuildTree = function(tag){
  #tag stat
  Cells=tag$celltype
  tag=tag$tag
  Tag = data.frame(table(tag$Tag))
  tag$Cell.num = Tag$Freq[match(tag$Tag,Tag$Var1)]
  tag_tab = acast(tag,Cell.BC~Tag)
  tag_tab[!is.na(tag_tab)] = 1
  tag_tab[is.na(tag_tab)] = 0
  #tag integrate
  cell_tab = data.frame(table(tag$Cell.BC))
  cell_tab = cell_tab[order(-cell_tab$Freq),]
  tags_all = lapply(as.character(cell_tab$Var1),function(x){sort(as.character(tag$Tag[tag$Cell.BC == x]))})
  tags_paste = sapply(tags_all,function(x){paste(x,collapse = "_")})
  tags_tab = data.frame(table(tags_paste))
  tags_tab$num = unlist(lapply(as.character(tags_tab$tags_paste),function(x){length(strsplit(x,split = "_")[[1]])}))
  tags_tab = tags_tab[order(-tags_tab$num,-tags_tab$Freq),]
  tags_uni = lapply(as.character(tags_tab$tags_paste),function(x){strsplit(x,split = "_")[[1]]})
  Tag_1 = Tag[Tag$Var1 %in% unlist(tags_uni[tags_tab$num == 1]),]
  Tag_1 = Tag_1[order(-Tag_1$Freq),]
  tags_uni[tags_tab$num == 1] = as.list(as.character(Tag_1$Var1))
  tags_tab[tags_tab$num == 1,] = tags_tab[tags_tab$num == 1,][match(as.character(Tag_1$Var1),
                                                                    as.character(tags_tab$tags_paste[tags_tab$num==1])),]

  #node build
  cluster_stat = function(i){
    x = tags_uni[[i]]
    n = tags_tab$num[i]
    tags_belone = NA
    for(y_ind in which(tags_tab$num<n)){
      y=tags_uni[[y_ind]]
      if(length(intersect(x,y))>0 & length(setdiff(y,x))==0){
        tags_belone = y_ind
        break
      }else{
        next
      }
    }
    return(tags_belone)
  }

  belons = sapply(which(tags_tab$num > 1),cluster_stat)
  belons[tags_tab$num == 1]= NA

  Tags = as.list(which(tags_tab$num == 1))
  names(Tags) = which(tags_tab$num == 1)
  nodes = as.list(c(1:length(tags_tab$num)))

  for(i in c(1:length(tags_tab$num))){
    n = belons[[i]]
    while(!is.na(n)){
      nodes[[i]] = c(n,nodes[[i]])
      n = belons[[n]]
    }
  }
  nodes_len = sapply(nodes,length)

  for(i in c(1:length(tags_tab$num))){
    nodes[[i]] = as.character(tags_tab$tags_paste)[nodes[[i]]]
    if(length(nodes[[i]])<max(nodes_len)){
      nodes[[i]][(length(nodes[[i]])+1):max(nodes_len)] = NA
    }else{
      next
    }
  }
  nodes=data.frame(do.call("rbind",nodes))
  names(nodes) = paste("N",as.character(1:ncol(nodes)),sep = "")
  nodes$pathString = do.call(paste, c("N0",nodes, sep="/"))
  nodes$pathString = gsub("/NA","",nodes$pathString)

  #save tree figure and rds
  population = as.Node(nodes)
  saveNetwork(diagonalNetwork(ToListExplicit(population, unname = TRUE),
                              margin = 10,fontSize = 8,width=15*dim(nodes)[2] ,height = 30*dim(nodes)[1]),
              file = "tree.html")

  saveRDS(population,"tree.rds")

  #save celltype tab
  cell_tab$tags = tags_paste
  if(!is.null(Cells)){
    cell_tab$celltype = Cells$Cell.type[match(as.character(cell_tab$Var1),Cells$Cell.BC)]
  }
  write.csv(cell_tab,"cell_tab.csv",row.names = FALSE,quote = FALSE)
  return(list(tree=population,info=cell_tab))

}


#' Title
#' @title PlotTree
#' @description Function to visualise the array generate tree
#' @param treeinfo List generate from BuildTree, for more see \code{\link[LinTInd]{BuildTree}}
#' @param data.extract (optional) If "FALSE" (default), will not return the indel's information, if it's "TRUE", the opposite
#' @param annotation (optional) If "TRUE" (default), the annotation of each cell barcodes have to be provided before, and a heatmap of cells' distribution for each array will be return
#' @param prefix (optional) Indels' prefix
#'
#' @return A list include a ggplot object, a dataframe show the distribution of each array contained in each group of cells (optional), and a dataframe to create the ggplot object
#' @export
#' @importFrom data.tree ToDataFrameNetwork FromDataFrameNetwork ToNewick
#' @importFrom reshape2 dcast melt
#' @importFrom purrr map_chr %>%
#' @importFrom ape rotateConstr
#' @importFrom ggtree ggtree geom_tiplab xlim_expand geom_facet read.tree
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 theme aes ggplot geom_ribbon geom_line scale_x_continuous xlab ylab theme_bw geom_segment geom_tile
#' @importFrom cowplot plot_grid
#' @importFrom BiocGenerics start end width type
#' @importFrom ggnewscale new_scale_color
#' @importFrom stringr str_extract str_extract_all
#' @importFrom dplyr summarise group_by
#' @importFrom rlist list.insert list.order
#' @examples
#' data("example_data",package = "LinTInd")
#' plotinfo<-PlotTree(treeinfo = treeinfo,data.extract = "TRUE",annotation = "TRUE")
#' plotinfo<-PlotTree(treeinfo = treeinfo,data.extract = "TRUE",annotation = "FALSE")
#'
PlotTree = function(treeinfo,data.extract=NULL,annotation=NULL,prefix=NULL){
  tree=treeinfo$tree
  celltab=treeinfo$info

  #1.trans data.tree to newick with internal node
  td = ToDataFrameNetwork(tree)
  from_node = unique(td$from)
  td[which(td$to %in% from_node),"to"] = paste("node.",td[which(td$to %in% from_node),"to"],sep = "")
  td = rbind(data.frame("from" = from_node, "to" = from_node),td)
  td$from = paste("node.",td$from,sep = "")
  td = td[which(td$to != "N0"),]
  td_node = FromDataFrameNetwork(td)
  tree_nwk = ToNewick(td_node)
  write(tree_nwk,"tree.nwk")


  #2.data build
  #(1)tree_nwk for plotting tree
  tree_nwk = read.tree("tree.nwk")
  tl = tree_nwk$tip.label

  #(2)tree_cell_plot for plotting cell composition
  tree_cell = celltab %>% group_by(tags, celltype) %>% summarise(count = sum(Freq))
  tree_cell = dcast(tree_cell,tags~celltype,value.var = "count")
  tree_cell[is.na(tree_cell)] = 0
  #normlize
  tree_cell[,-1] = log(tree_cell[,-1]+1)
  tree_cell[,-1] = t(apply(tree_cell[,-1], 1, function(x) {x/sum(x)}))
  tree_cell_l = melt(tree_cell)
  #tree_cell_plot = tree_cell_l

  #option for cell composition bar plot
  # tree_cell_l = melt(tree_cell)
  # tree_cell_plot = left_join(data.frame(id=tl),tree_cell_l,by=c("id" = "tags"))

  #option for
  # tree_cell_plot$variable = as.character(tree_cell_plot$variable)
  # tree_cell_plot[which(tree_cell_plot$value>0),3] = tree_cell_plot[which(tree_cell_plot$value>0),2]
  # tree_cell_plot[which(tree_cell_plot$value==0),3] = NA

  #(3)indel for plotting pattern
  #library(stringr)

  indel = data.frame(id=NA,start=NA,end=NA,width=NA,type=NA)
  indel = indel[-1,]
  #build indel data frame
  for(i in c(1:length(tl))){
    x = tl[i]
    x_str = unlist(strsplit(x,split = "_"))
    del = NULL
    ins = NULL
    if(!is.null(prefix)){
      for (k in 1:length(prefix)) {
        x_pre = x_str[[1]][grepl(prefix[k], x_str[[1]])]
        if(length(x_pre)>0){
          x_pre_del = x_pre[grep("D",x_pre)]
          if(length(x_pre_del) > 0){
            x_pre_del = str_extract_all(x_pre_del, "[0-9]+")
            for (j in 1:length(x_pre_del)) {
              line = as.numeric(x_pre_del[[j]])
              del = rbind(del, data.frame(id=x,start=line[length(line)] + (k-1)*300,
                                          end=line[length(line)]+line[length(line)-1] + (k-1)*300,width=line[length(line)-1],
                                          type="deletion"))
            }

          }
          x_pre_ins = x_pre[grep("I",x_pre)]
          if(length(x_pre_ins) > 0){
            x_pre_ins = str_extract_all(x_pre_ins, "[0-9]+")
            for (j in 1:length(x_pre_ins)) {
              line = as.numeric(x_pre_ins[[j]])
              ins = rbind(ins, data.frame(id=x,start=line[length(line)] + (k-1)*300,
                                          end=line[length(line)]+line[length(line)-1] + (k-1)*300,width=line[length(line)-1],
                                          type="insertion"))
            }
          }

        }
      }
    }else{
      x_pre_del = x_str[grep("D",x_str)]
      if(length(x_pre_del) > 0){
        x_pre_del = str_extract_all(x_pre_del, "[0-9]+")
        for (j in 1:length(x_pre_del)) {
          line = as.numeric(x_pre_del[[j]])
          #del = rbind(del, data.frame(id=x,start=line[length(line)] + (k-1)*300,
                                      #end=line[length(line)]+line[length(line)-1] + (k-1)*300,width=line[length(line)-1],
                                      #type="deletion"))
          del = rbind(del, data.frame(id=x,start=line[length(line)],
                                      end=line[length(line)]+line[length(line)-1],width=line[length(line)-1],
                                      type="deletion"))
        }

      }
      x_pre_ins = x_str[grep("I",x_str)]
      if(length(x_pre_ins) > 0){
        x_pre_ins = str_extract_all(x_pre_ins, "[0-9]+")
        for (j in 1:length(x_pre_ins)) {
          line = as.numeric(x_pre_ins[[j]])
          #ins = rbind(ins, data.frame(id=x,start=line[length(line)] + (k-1)*300,
                                      #end=line[length(line)]+line[length(line)-1] + (k-1)*300,width=line[length(line)-1],
                                      #type="insertion"))
          ins = rbind(ins, data.frame(id=x,start=line[length(line)],
                                      end=line[length(line)]+line[length(line)-1],width=line[length(line)-1],
                                      type="insertion"))
        }
      }
    }
    indel = rbind(indel,del,ins)
  }


  #(4)reorder tree
  SortTree = function(label){

    #construct data structure
    taglst = list()
    for (i in 1:length(label)) {
      tchr = label[i]
      tline = strsplit(tchr,"_")
      element = list("labels" = tchr,"tags" = tline[[1]], "tag_number" = length(tline[[1]]))
      taglst[[i]] = element
    }
    taglst.sort = taglst[list.order(taglst,-tag_number)]
    for (i in 2:(length(taglst.sort) - 1)) {
      point = taglst.sort[[i]]
      score = 0
      for (j in (i+1):length(taglst.sort)) {
        queue = taglst.sort[[j]]
        qscore = length(intersect(point$tags,queue$tags))

        if(qscore > score){
          score = qscore
          taglst.sort[[j]] = NULL
          taglst.sort = list.insert(taglst.sort,i+1,queue)
        }
      }
    }


    #result
    label.sort = map_chr(taglst.sort, 1)
  }
  label.sort =  SortTree(tree_nwk$tip.label)
  tree_data_sort = ape::rotateConstr(tree_nwk, label.sort)

  #3.plot integrated tree
  #set expand can change width of the tree
  p = ggtree(tree_data_sort,size = 0.1, ladderize=FALSE) + geom_tiplab(size = 0.3) +
    xlim_expand(c(0, 200),panel = "Tree")

  if(annotation=="FALSE"){
    p1=p + geom_facet(panel = "indel pattern", data = indel,
                      geom = geom_segment,
                      mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
                      size = 0.3)
  }else{
    p1 = p + geom_facet(panel = "cells",data = tree_cell_l,
                        geom =  geom_tile,
                        mapping = aes(x = as.numeric(as.factor(variable)),fill = value,color = variable)) +
      #    scale_fill_viridis_d(option="D", name="discrete\nvalue") +
      scale_fill_gradient(low = "white",high = "#440130")  +
      new_scale_color() +
      geom_facet(panel = "indel pattern", data = indel,
                 geom = geom_segment,
                 mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
                 size = 0.3)
  }
  if(data.extract=="TRUE" & annotation=="TRUE"){
    return(list(p=p1,tagsinfo=tree_cell,indelinfo=indel))
  }else if(data.extract=="TRUE" & is.null(annotation)){
    return(list(p=p1,tagsinfo=tree_cell,indelinfo=indel))
  }else if(data.extract=="TRUE" & annotation=="FALSE"){
    return(list(p=p1,indelinfo=indel))
  }else{
    return(p1)
  }

  # p1 = facet_plot(p,panel = "cell composition",data = tree_cell_l,
  #                 geom = geom_barh,
  #                 mapping = aes(x = value,fill = variable),
  #                 stat="identity") %>%
  #   facet_plot(panel = "indel pattern", data = indel,
  #              geom = geom_segment,
  #              mapping = aes(x = start, xend = end, y =y, yend = y,color=type),
  #              size = 1)


  #change facet grid
  #gt = ggplot_gtable(ggplot_build(p1))
  #gt$widths[7] = 0.5*gt$widths[7] # in this case it was colmun 7 - reduce the width by a half
  #gt$widths[9] = 0.5*gt$widths[9] # in this case it was colmun 7 - reduce the width by a half
  #grid.draw(gt) # plot with grid draw

  #ggsave(gt,filename = paste(outname,"pdf",sep = "."),height = 100,width = 15,units = "cm")
}


#' Title
#' @title IndelPlot
#' @description Return 2 line charts, show the probability of insertion and deletion at each site
#' @param cellsinfo List generate from IndelIdents, for more see \code{\link[LinTInd]{IndelIdents}}
#'
#' @return 2 line charts
#' @export
#'
#' @importFrom IRanges IRanges
#' @importFrom ggplot2 theme aes ggplot geom_ribbon geom_line scale_x_continuous xlab ylab theme_bw
#' @importFrom cowplot plot_grid
#' @importFrom BiocGenerics start end width type
#'
#' @examples
#' data("example_data",package = "LinTInd")
#' IndelPlot(cellsinfo = cellsinfo)
#'

IndelPlot<-function(cellsinfo){
  scar=cellsinfo$cutsite
  indel.coverage<-cellsinfo$indel.coverage
  INDEL_ranges<-cellsinfo$indel
  del_ranges<-unlist(lapply(INDEL_ranges,function(x){x[1]}))
  ins_ranges<-unlist(lapply(INDEL_ranges,function(x){x[2]}))
  del_r<-del_ranges[[1]]
  for(i in 2:length(del_ranges)){
    del_r<-c(del_r,del_ranges[[i]])
  }
  ins_r<-ins_ranges[[1]]
  for(i in 2:length(ins_ranges)){
    ins_r<-c(ins_r,ins_ranges[[i]])
  }
  if(any(is.null(indel.coverage),indel.coverage=="Accurate")){
    scarref<-data.frame(scar=NA,type=NA)
    scarref<-scarref[-1,]
    for(i in c(2:dim(scar)[1])){
      scarref_sub<-data.frame(scar=c(scar$start[i]:(scar$end[i])),type=scar$indx[i])
      scarref<-rbind(scarref,scarref_sub)
    }
  }else{
    scarref<-data.frame(scar=NA,type=NA)
    scarref<-scarref[-1,]
    startP<-scar$start[1]
    for(i in c(2:dim(scar)[1])){
      scarref_sub<-data.frame(scar=c(scar$start[i]:(scar$end[i]))+startP,type=scar$indx[i])
      scarref<-rbind(scarref,scarref_sub)
    }
  }
  scarref$type<-as.character(scarref$type)
  all_site_stat<-function(i,dat){
    return(c(start(dat[i]):end(dat[i])))
  }
  #all_site_per_stat<-function(i,dat){
  #	return(length(which(dat==i))/length(INDEL_ranges))
  #}

  del_allsite<-unlist(lapply(seq(del_r),all_site_stat,dat=del_r))
  del_allsite_per<-data.frame(table(del_allsite))
  del_allsite_per$Freq<-del_allsite_per$Freq/length(INDEL_ranges)
  names(del_allsite_per)[1]<-"Site"
  del_allsite_per<-del_allsite_per[del_allsite_per$Site %in% c(1:scar$end[1]),]
  del_allsite_per$Site<-as.numeric(as.character(del_allsite_per$Site))
  if(length(setdiff(c(1:scar$end[1]),del_allsite_per$Site))>0){
    del_allsite_per<-rbind(del_allsite_per,data.frame(Site=setdiff(c(1:scar$end[1]),del_allsite_per$Site),Freq=0))
  }else{
    del_allsite_per<-del_allsite_per
  }
  if(indel.coverage=="Accurate"){
    del_allsite_per<-del_allsite_per[del_allsite_per$Site<=max(scarref$scar+1),]
  }
  p_del<-ggplot()+
    geom_ribbon(data=scarref,aes(x=scar,ymin=0,ymax=max(del_allsite_per$Freq),fill=type),alpha=0.1)+
    geom_line(data=del_allsite_per,aes(x=Site,y=Freq),size=1,colour="lightskyblue")+theme(legend.position="none")+
    scale_x_continuous(breaks=c())+xlab("")+ylab("Deletion Frequency")+theme_bw()

  ins_allsite<-unlist(lapply(seq(ins_r),all_site_stat,dat=ins_r))
  ins_allsite_per<-data.frame(table(ins_allsite))
  ins_allsite_per$Freq<-ins_allsite_per$Freq/length(INDEL_ranges)
  names(ins_allsite_per)[1]<-"Site"
  ins_allsite_per<-ins_allsite_per[ins_allsite_per$Site %in% c(1:scar$end[1]),]
  ins_allsite_per$Site<-as.numeric(as.character(ins_allsite_per$Site))
  if(length(setdiff(c(1:scar$end[1]),ins_allsite_per$Site))>0){
    ins_allsite_per<-rbind(ins_allsite_per,data.frame(Site=setdiff(c(1:scar$end[1]),ins_allsite_per$Site),Freq=0))
  }else{
    ins_allsite_per<-ins_allsite_per
  }
  if(indel.coverage=="Accurate"){
    ins_allsite_per<-ins_allsite_per[ins_allsite_per$Site<=max(scarref$scar+1),]
  }
  p_ins<-ggplot()+
    geom_ribbon(data=scarref,aes(x=scar,ymin=0,ymax=max(ins_allsite_per$Freq),fill=type),alpha=0.1)+
    geom_line(data=ins_allsite_per,aes(x=Site,y=Freq),size=1,colour="indianred1")+theme(legend.position="none")+
    scale_x_continuous(breaks=c())+xlab("")+ylab("Insertion Frequency")+theme_bw()
  return(plot_grid(p_del,p_ins,nrow = 2))
}


#' Title
#' @title TagDist
#' @description If the cell barcode and the anntation of each cell are provided, this function can calculate the relationship between each cell type in three way
#' @param tag List generate from TagProcess, for more see \code{\link[LinTInd]{TagProcess}}
#' @param method Denote which method to use:
#'
#' - "Jaccard"(default): calculate the weighted jaccard similarity of indels between each pair of groups;
#' - "P": right-tailed test, compare the Indels intersection level with the hypothetical  result generated from random editing, and the former is expected to be significantly higher than the latter;
#' - "spearman": Spearman correlation of indels between each pair of groups
#'
#' @return 2 figures are saved to show the distribution of INDEL and the relationship between groups respectively, the matrix of the relationship between groups is returned
#' @export
#'
#' @importFrom stats pnorm sd cor
#' @importFrom pheatmap pheatmap
#' @importFrom reshape2 acast
#'
#' @examples
#' tag_dist=TagDist(tag,method = "spearman")
#' tag_dist=TagDist(tag)
#' tag_dist=TagDist(tag,method = "Jaccard")
#' tag_dist=TagDist(tag,method = "P")
#' tag_dist=TagDist(tag,method = "spearman")
#'
TagDist<-function(tag,method=NULL){

  one_jac_stat2<-function(l,VBC1,clone_stat_data){
    VBC2<-as.numeric(clone_stat_data[,l])
    num <- sum(sapply(1:length(VBC1), function(x)(min(VBC1[x],VBC2[x]))))
    den <- sum(sapply(1:length(VBC1), function(x)(max(VBC1[x],VBC2[x]))))
    return(num/den)
  }
  one_jac_stat1<-function(i,clone_stat_data){
    VBC1<-as.numeric(clone_stat_data[,i])
    return(unlist(lapply(clu,one_jac_stat2,VBC1=VBC1,clone_stat_data=clone_stat_data)))
  }
  one_op_stat2<-function(y,x,jac,jac_sample){
    jac_real<-jac[x,y]
    jac_pred<-as.numeric(unlist(lapply(jac_sample,function(z){z[x,y]})))
    ob<-jac_real/mean(jac_pred)
    zscore<-(jac_real-mean(jac_pred))/sd(jac_pred)
    p<-pnorm(zscore,lower.tail = FALSE)
    return(c(ob,p))
  }
  one_op_stat1<-function(x,jac,jac_sample,clu){
    return(lapply(clu,one_op_stat2,x=x,jac=jac,jac_sample=jac_sample))
  }
  tag<-tag$tag
  clone_tab<-acast(tag,Tag~Cell.type)
  clu<-colnames(clone_tab)
  clone_tab<-clone_tab[apply(clone_tab,1,sum)>=2,]
  annotation_col = data.frame(Group=factor(unlist(str_extract_all(row.names(clone_tab), "[D|I]+"))))
  row.names(annotation_col)<-rownames(clone_tab)
  pheatmap(t(clone_tab),cluster_cols = FALSE,cluster_rows = FALSE,show_colnames = FALSE,border=FALSE,annotation_col = annotation_col,scale = "row",filename = "tag_heatmap_scale.pdf",width = 5,height = 3)

  if(any(is.null(method),method=="Jaccard")){
    all_jac<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab)),ncol=length(clu)))
    names(all_jac)<-clu
    row.names(all_jac)<-clu
    pheatmap(all_jac,file="group_dist.pdf",width = 4.9,height = 4.5)
    return(all_jac)
  }else if(method=="P"){
    all_jac<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab)),ncol=length(clu)))
    names(all_jac)<-clu
    row.names(all_jac)<-clu
    jac_sample<-list()
    for(t in c(1:500)){
      tag$tags_sample<-sample(tag$Tag,length(tag$Tag))
      clone_tab_sample_sub<-acast(tag,tags_sample~Cell.type)
      jac_sample_sub<-data.frame(matrix(unlist(lapply(clu,one_jac_stat1,clone_stat_data=clone_tab_sample_sub)),ncol=length(clu)))
      names(jac_sample_sub)<-clu
      row.names(jac_sample_sub)<-clu
      jac_sample<-c(jac_sample,list(jac_sample_sub))
    }

    #plot pvalue
    all_ob_p<-lapply(clu,one_op_stat1,jac=all_jac,jac_sample=jac_sample,clu=clu)
    all_ob<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==1]]}))
    all_ob<-data.frame(all_ob)
    all_p<-do.call("rbind",lapply(all_ob_p,function(x){unlist(x)[c(1:length(unlist(x)))[c(1:length(unlist(x)))%%2==0]]}))
    all_p<-data.frame(all_p)
    names(all_ob)<-clu
    row.names(all_ob)<-clu
    names(all_p)<-clu
    row.names(all_p)<-clu

    all_p[is.na(all_p)] = 0
    pheatmap(-log2(all_p+min(all_p[all_p!=0])),main = "log2(P.val)",file="group_dist.pdf",width = 4.9,height = 4.7)
    return(all_p)
  }else if(method=="spearman"){
    clone_tab_new<-t(apply(clone_tab,1,function(x){x/sum(x)}))
    clone_tab_new<-apply(clone_tab_new,2,function(x){x/sum(x)})
    cor_s<-cor(clone_tab_new,method = "spearman")
    pheatmap(cor_s,file="group_dist.pdf",width = 4.9,height = 4.5)
    return(cor_s)
  }
}
