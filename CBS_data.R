library(DNAcopy)
"CBS_data"<-function(){
  args<-commandArgs(T)
  name = paste(args[1],"txt",sep=".")
  data = read.table(name)
  #data = data[-1, -(1:2)]
  #data = log2(data/2)
  head = matrix(0, nrow(data), 3)
  head[,1] = 1
  head[,2] = 1:nrow(data)
  head[,3] = 1:nrow(data)

  chrom <- rep(1, nrow(data))
  maploc <- 1:nrow(data)
  seg.file_g = matrix(0,1,6)
  seg.file_g_one = matrix(0,1,6)
  seg.file = matrix(0, nrow(data), 1)

  stac_amp = matrix(0, 1, nrow(data))
  stac_amp[1,] = 1:nrow(data)
  stac_amp_one = matrix(0, 1, nrow(data))

  stac_del = matrix(0, 1, nrow(data))
  stac_del[1,] = 1:nrow(data)
  stac_del_one = matrix(0, 1, nrow(data))


  for (j in 1:ncol(data)){

    smoothed.data <- smooth.CNA(CNA(data[,j],chrom,maploc))
    seg<- segment(smoothed.data)
    for (k in 1:length(seg$output$loc.start)){
      seg.file_g_one[1,1]=j
      seg.file_g_one[1,2]=1
      seg.file_g_one[1,3]=seg$output$loc.start[k]
      seg.file_g_one[1,4]=seg$output$loc.end[k]
      seg.file_g_one[1,5]=seg$output$num.mark[k]
      seg.file_g_one[1,6]=seg$output$seg.mean[k]
      seg.file_g=rbind(seg.file_g,seg.file_g_one)
      seg.file_g_one=matrix(0,1,6)

    }


  }
  seg.file_g = seg.file_g[-1,]

  out.file=getwd()
  out.file=paste(out.file,'seg',sep="/")

  write.table(seg.file_g,file=out.file,row.names=F,col.names=F,quote=F,sep="\t")
  
}
CBS_data()




  


