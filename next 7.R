setwd("/Users/yuwenyi/SP-workshop")
a <- scan("1581-0.txt",what="character",skip=156)
n <- length(a)
a <- a[-((n-2909):n)] ## strip license

#step 4
split_punct<- function(p) {
ip<-grep(p, a, fixed=TRUE)
xa<-gsub(p,"",a,fixed=TRUE)
xxa<-rep("",length(xa)+length(ip))
iia<-ip+1:length(ip)
xxa[iia]<-paste(p)
xxa[-iia]<-xa
return(xxa)
}

#step 5
punc<-c(",", ".", ";", "!", ":", "?")
for (p in punc) {
  a<-split_punct(p)
}

#step 6
a<-tolower(a)
b<-unique(a)#finds the vector of unique words
c<-match(a,b)#finds the vector of indicies
freq<-tabulate(c)#indicates the number of times each unique word occurs in the text
oder_index<-order(freq,decreasing=TRUE)[1:1000]#？If we need to certify why 1000 is the most common?not 1001 
b_index<-sort(oder_index)
b<-b[b_index]
b

