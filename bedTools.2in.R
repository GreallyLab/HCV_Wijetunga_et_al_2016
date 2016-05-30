# Function to run bedTools in R from:
# http://zvfak.blogspot.com/2011/02/calling-bedtools-from-r.html

bedTools.2in<-function(functionstring="intersectBed",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=try(read.table(out,header=F),silent=T)
  
  unlink(a.file);unlink(b.file);unlink(out)
  if (is(res, "try-error")) return(mat.or.vec(0,3)) else return(res)

}