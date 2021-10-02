setwd("/Users/yuwenyi/SP-workshop")
a <- scan("1581-0.txt",what="character",skip=156)
n <- length(a)
a <- a[-((n-2909):n)] ## strip license

split_punct<- function(p) {
ip<-grep(p, a, fixed=TRUE)
xa<-gsub(p,"",a,fixed=TRUE)
xxa<-rep("",length(xa)+length(ip))
iia<-ip+1:length(ip)
xxa[iia]<-paste(p)
xxa[-iia]<-xa
return(xxa)
}
punc<-c(",", ".", ";", "!", ":", "?")
for (p in punc) {
  a<-split_punct(p)
}