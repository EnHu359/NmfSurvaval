library(survminer)
library(survival)
#clinicaldata samples OS.time OS.state DFS.time DFS.state
#groupfile samples group
#save files in current directiory
mutil_sur<-function(group.data,clinical.data,figure_name,figur_any=F){
  rootpath=getwd()
  #读取分组信息
  groups<<-group.data
  groups$samples<<-gsub("^X|T$|N$","",groups$samples)
  #print(groups)
  roman<-paste("Subtype",c("I","II","III","IV","V","VI","VII","VIII","IX","X"),sep = " ")
  groups$group<-as.numeric(groups$group)
  for(i in levels(as.factor(groups$group))){
    groups$group[groups$group==i]<-paste0(roman[as.numeric(i)]," (n=",sum(groups$group==i),")")
  }
  
  #创建OS数据
  clinical.data<<-clinical.data[which(clinical.data$OS.state!=2),]
  #合并分组信息和生存信息
  clinical_groups<<-merge(clinical.data,groups,by="samples")
  print(clinical_groups)
  #构建生存对象
  formula_OS<-as.formula(paste0('Surv(time=OS.time,event=OS.state)~group'))
  fit_OS<-surv_fit(formula_OS,data = clinical_groups)
  #logrnak检验
  setff_os<-survdiff(formula_OS,data=clinical_groups)
  #计算卡方P值
  P.val.os<- 1 - pchisq(setff_os$chisq,length(setff_os$n) - 1)
  
  #创建DFS数据
  clinical.data<-clinical.data[which(clinical.data$DFS.state!=2),]
  clinical_groups<-merge(clinical.data,groups,by="samples")
  formula_DFS<-as.formula(paste0('Surv(time=DFS.time,event=DFS.state)~group'))
  fit_DFS<-surv_fit(formula_DFS,data = clinical_groups)
  setff_DFS<-survdiff(formula_DFS,data=clinical_groups)
  P.val.DFS<-1-pchisq(setff_DFS$chisq,length(setff_DFS$n) - 1)
  #print(paste0("OS:",P.val.os," DFS:",P.val.DFS,figurename))
  #只有当两个p值都小于0.05时才绘图
  if((P.val.os<=0.05 & P.val.DFS <=0.05) | figur_any){
    #图片标题设置未什么好呢
    return("T")
    p1<<-ggsurvplot(fit_OS, data = clinical_groups,
                    surv.median.line = "hv",
                    #conf.int = TRUE,
                    risk.table = TRUE,
                    pval = TRUE,
                    palette = "aaas",
                    legend="top",
                    #legend = c(0.8,0.9),
                    legend.title="Group",
                    break.x.by = 300,
                    xlab = "OS time",
                    font.x="bold",
                    font.y="bold",
                    ggtheme = theme_bw(),
                    
                    
                    
    )
    
    p2<<-ggsurvplot(fit_DFS, data = clinical_groups,
                    surv.median.line = "hv",
                    #conf.int = TRUE,
                    risk.table = TRUE,
                    pval = TRUE,
                    palette = "aaas",
                    legend="top",
                    #legend = c(0.8,0.9),
                    legend.title="Group",
                    break.x.by = 300,
                    xlab = "DFS time",
                    font.x="bold",
                    font.y="bold",
                    ggtheme = theme_bw()
                    
    )
    png(paste0(figure_name,"_OS_",round(P.val.os),".png"),width = 1062,height = 802,pointsize = 72)
    print(p1)
    dev.off()
    png(paste0(figure_name,"_DFS_",round(P.val.DFS),".png"),width = 1062,height = 802,pointsize = 72)
    print(p2)
    dev.off()
  }
}