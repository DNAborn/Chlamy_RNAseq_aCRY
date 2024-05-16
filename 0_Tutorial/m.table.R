samplelist <- {}
mappingrates <- {}
for (i in list.files(path = quantdir)){
  print(i)
  si <- paste("",str_sub(i,11,-7),sep = "")
  si
  samplelist <- c(samplelist,si)
  f <- readLines(paste(quantdir,i,"logs/salmon_quant.log", sep="/"))
  line <- grep("Mapping rate = ",f,value=TRUE)
  sl <- str_length(line)
  sl
  notime <- substring(line,30,sl)
  notime
  manual <- substring(line,sl-7,sl-1)
  val <- as.numeric(str_extract(notime,"[0-9.]+"))
  val
  valr<-round(val, digits=2)
  print(paste("Mapping rate of ",si," is: ",valr," %"))
  mappingrates <- c(mappingrates,valr)
}
sample.table$mappingrates <- mappingrates
# Make table

m.table <- data.frame(sample.table,mappingrates)