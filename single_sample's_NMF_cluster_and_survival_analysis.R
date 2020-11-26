library("NMF")
source("./mutil_sur.R")
#use Rscript to get outer args
#usage Rscript single samples's NMF cluster and survival.R filename
args <- commandArgs(T)
datasname<-args[1]

filename=strsplit(datasname,"\\.")[[1]][1]
#This function is used to create a directory
create.in<-function(dir.name,change.dir=T){
  if(dir.exists(dir.name)){
    print(paste(dir.name,"exists remove it"))
    unlink(dir.name,recursive = T)
    print(paste("create",dir.name))
    dir.create(paste0("./",dir.name))
  }else{
    print(paste("create",dir.name))
    dir.create(paste0("./",dir.name))
  }
  if(change.dir==T){
    print(paste("change workdir to",dir.name))
    setwd(paste0("./",dir.name))
  }
}
#read in expressiondata
datas<-read.csv(file = datasname,row.names = 1)
#read in clinical data
clincal<-read.csv("clinical_time.csv")

#outprefix=expdata_mads_k_method
single.NMF.and.survival<-function(exp_data,k,method,out.prefix){
  out.prefix = paste0(out.prefix,"_k=",k,"_",method)
  create.in(paste0("./",out.prefix),change.dir = F)
  nmf.fit<-nmf(exp_data,rank = k,method,nrun=100,seed=123456,.opt="v")
  groups<-data.frame(group=predict(nmf.fit))
  groups$samples<-row.names(groups)
  #将分组信息输出
  write.csv(groups,paste0("./",out.prefix,"/groups.csv"))
  figure.name=(paste0("./",out.prefix,"/",out.prefix))
  try(nmf.survival.result<-mutil_sur(group.data = groups,clinical.data = clincal,
                                     figure_name = figure.name))
  save(file = paste0("./",out.prefix,"/",out.prefix,".Rdata"),list = "nmf.fit")
  
  png(paste0("./",out.prefix,"/",out.prefix,"_consensusMap.png"),width = 1062,height = 802,pointsize = 72)
  consensusmap(nmf.fit)
  dev.off()
  
}
create.in("mads_file",change.dir = F)
mad_data=data.frame(mads=apply(datas,1,mad),row.names = row.names(datas))
write.csv(mad_data,"./mads_file/mads_file.csv")
methods=c("nsNMF","offset",
         "brunet","KL","lee","snmf/r","snmf/l",
         "Frobenius")
ks=c(2,3,4,5)
for(percent in seq(0.3,1,0.05)){
  mad_data_mad <- datas[row.names(mad_data)[order(mad_data$mads,decreasing = T)[1:(dim(datas)[1]*percent)]],]
  #将mads文件保存到mads_file里面
  write.csv(mad_data_mad,paste0("./","mads_file/",percent,".csv"))
  for(method in methods){
    lapply(ks,function(x) single.NMF.and.survival(exp_data = mad_data_mad,method,k = x,
                                                  out.prefix =paste(filename,"mads",percent,sep = "_")))
  }

}
