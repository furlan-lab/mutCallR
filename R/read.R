#' @export
#' @title Glimpse into mutCaller counts.txt.gz file for the names of variants and number of reads
#' @description See what variants are in a counts.txt.gz output file from mutCaller
#' @returns a table of counts for each variant contained in the mutCaller output file
#' @param file path to file name (only 1)

glimpse_countsfile<-function(file){
  d<-readr::read_delim(file, col_names = c("cb", "umi", "seq", "loc", "name", "call", "count"))
  return(table(d$name))
}




#' @export
#' @title Add mutcaller data to a Seurat or Monocle3 object
#' @description This function aims to import mutCaller data into Seurat/Monocle3 objects.  Output files for mutCaller will include a counts.txt.gz file.
#' This (These) file location(s) should be provided using the 'counts_file' argument.  If specific prefixes have been cellular barcode (CB) in your seurat/monocle
#' object, you can provide them uses the 'prefixes' for each mutCaller counts.txt.gz file (optional).  The 'name' of the variant to be extracted from the counts.txt.gz file
#' is a required argument and this function will need to be invoked for each variant contained in the counts.txt.gz file.
#'
#' Use the glimpse_countsfile function in this package, to quickly see all the names of variants contained in a mutCaller counts.txt.gz file and the number of unique CB/UMI combinations
#' associated with each variant.  After loading counts data, read_mutcaller will classify each cell as 'Ref', 'Query', or 'Other'.  Note that calls for 'Other' (that is variants that
#' neither fit the Ref or the Query) are omitted by default, but can be retained by setting the use_other argument to TRUE.  The read_thresh argument is used
#' to omit calls with substantiating reads fewer than this threshold.  Note that in cases where multiple reads with the same CB/UMI have discrepent calls,
#' read_mutcaller will prioritize the call with the highest readcount supporting the variant.  This calculation includes a consideration of cases in which
#' reads that have the same CB, but UMIs that differ by only 1 bp in their UMIs.  In these situations, read_mutcaller calls will lump these likely errors in UMI sequence with their
#' closely related neighbors.
#'
#' @returns a seurat/monocle3 object with a new column in the metadata( or colData) with cell classificatiuoon for a given variant (previously counted using mutcaller)
#' @param obj input Seurat or Monocole3 object (optional).  If not provided this function will process all mutCaller counts and return cell calls
#' @param count_files a vector of one or more paths to counts.txt.gz file output by mutcaller; e.g. '/home/user/mutcaller_run/out_star/counts.txt.gz'
#' @param name the name of the variant to be added to the object (or returned as raw counts).  This is found in the 5th column of the counts.txt.gz file
#' @param prefixes the optional name of prefixes to be added to cell barcodes.  The length of this cvector of prefixes should be the same as the length of count_files
#' @param cores num of cores to use for processing counts
#' @param read_thresh, the minimum number of reads required to support a molecule call
#' @param use_other, retain information regarding 'Other' cell calls in output; default is F
#' @param query_top, option for factoring the cell call output such that the levels of the factor enable plotting query cells on top; default = T
#' @importFrom pbmcapply pbmclapply
#' @importFrom readr read_delim
#'
read_mutcaller<-function(obj=NULL, count_files, name, prefixes=NULL, cores=1, read_thresh=5,  use_other=F, query_top=T){
  if (length(count_files)!=length(prefixes)){stop("Error: number of count files and number of prefixes do not match!")}
  #count_file<-count_files[1]
  datlist<-pbmclapply(count_files, function(count_file) {
    # d<-fread(count_file)
    d<-readr::read_delim(count_file, col_names = c("cb", "umi", "seq", "loc", "name", "call", "count"))
    #colnames(d)<-c("cb", "umi", "seq", "loc", "name", "call", "count")
    d$cb<-substr(d$cb, 1, 16)
    d<-d[d$count>read_thresh,]
    d<-d[d$name==name,]
    datlist<-lapply(unique(d$cb), function(ucb){
      ##i<-1
      # i=i+1
      # ucb<-unique(d$cb)[i+1]
      slice<-d[d$cb %in% ucb,]
      # slice
      if(nrow(slice)>1){
        slice<-LDfunc(slice)
        countSlice(slice)
      }else{
        if(slice$call=="Ref"){
          data.frame(cb=slice$cb, name=slice$name, Ref=1, Query=0, Other=0)
        }else if(slice$call=="Query"){
          data.frame(cb=slice$cb, name=slice$name, Ref=0, Query=1, Other=0)
        } else {
          data.frame(cb=slice$cb, name=slice$name, Ref=0, Query=0, Other=1)
        }
      }
    })
    do.call(rbind, datlist)
  }, mc.cores = cores)
  if(!is.null(prefixes)){
    datlist<-lapply(1:length(datlist), function (n){
      if(is.null(datlist[[n]])){
        datlist[[n]]
      }else{
        datlist[[n]]$cb<-paste0(prefixes[n], datlist[[n]]$cb)
        datlist[[n]]
      }
    })
  }
  alldat<-do.call(rbind, datlist)
  if(is.null(obj)){
    return(alldat)
  }
  alldat<-call_cell(alldat, use_other)
  if(class(obj)=="Seurat"){
    clean_cb<-clean_cb_format(Cells(obj))
  } else if(class(obj)=="cell_data_set"){
    clean_cb<-colnames(obj)
  } else {
    stop("Object class is not known")
  }
  final_call<-alldat$call[match(clean_cb, alldat$cb)]
  if(query_top){
    final_call<- factor(final_call, levels=c(NA, "Ref", "Ref+Query", "Query"))
  }
  obj[[name]]<-final_call
  obj
}

countSlice<-function(slice){
  outslice<-data.frame(cb=slice$cb[1], name = slice$name[1], Ref=0, Query=0, Other=0)
  wd<-table(LDfunc(slice)$call)
  for(i in 1:length(wd)){
    outslice[[names(wd)[i]]]<-as.numeric(wd)[i]
  }
  outslice
}

#' @importFrom stringdist stringdistmatrix
LDfunc<-function(slice, ldthresh=2){
  #this function takes a set of data with column umi and column count and tests for similarity between the umis using levenshtein distance (ld).  For rows of data with ld similarity < ldthresh, this function will keep only the row with the highest count value return the original data minus the eliminated row.
  mdist<-stringdistmatrix(slice$umi, slice$umi, method = "lv")
  losers<-vector()
  dim<-length(slice$umi)
  for(i in 1:dim){
    for(j in 1:dim){
      if(i > j){
        if(mdist[i,j]<ldthresh){
          #print(paste0("i = ", i, "j = ", j))
          loser<-which.min(c(slice$count[i], slice$count[j]))
          losers<-c(losers, c(i,j)[loser])
        }
      }
    }
  }
  if(length(losers)>0){
    slice[-unique(losers),]
  }else{
    slice
  }
}



call_cell<-function(data, use_other=F){
  if(!all(colnames(data)==c("cb", "name","Ref","Query","Other"))){stop("Data does not appear to be in the correct format")}

  if(use_other){
    collapsed_by_cell<-apply(data[,3:5], 1, function(x) {
      names(which(x>0))
    })
    return (data.frame(cb=data$cb, call=sapply(collapsed_by_cell, function(entry) paste0(entry, collapse="+"))))
  } else {
    collapsed_by_cell<-apply(data[,3:4], 1, function(x) {
      names(which(x>0))
    })
    df<-data.frame(cb=data$cb, call=sapply(collapsed_by_cell, function(entry) paste0(entry, collapse="+")))

    return (df[df$call!="",])

  }

}

#' @importFrom stringr str_sub
#' @importFrom stringr str_length
clean_cb_format<-function(cb){
  #check for "-N"
  if(names(table(str_sub(cb, str_length(cb)-1, str_length(cb))))=="-1"){
    message("'-1' suffix found at the end of all cell barcodes; temporarily removing these for matching")
    return(str_sub(cb, 1, str_length(cb)-2))
  } else if(names(table(str_sub(cb, str_length(cb)-1, str_length(cb))))){

  } else if (all(grepl("-", t))){
    message("'-N' suffix found at the end of all cell barcodes; temporarily removing these for matching")
    return(str_sub(cb, 1, str_length(cb)-2))
  } else if (all(grepl("-", t2))){
    message("No suffixes found at the end of all cell barcodes; matching as provided")
    return(cb)
  }
}
