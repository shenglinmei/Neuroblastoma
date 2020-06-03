source('/home/meisl/bin/FunctionLib/Lib/pagodaLib.r')

p2=readRDS('~/Workplace/neuroblastoma//allp2.rds')

load('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/Neuron/Neu.Mar18.RData')


conM=readRDS('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/MES/Fib_conos.rds')

conM2=readRDS('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/MES/Fib.noCycle_conos.rds')

c1=conM$clusters$leiden$groups
c2=conM2$clusters$leiden$groups

c3=readRDS('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/Fib.cluster.rds')

table(annot[names(c3)])


con$plotGraph(groups=c1,plot.na=F)
con$plotGraph(groups=c2,plot.na=F)
con$plotGraph(groups=c3,plot.na=F)



load('F3new.RData')

table(neu4)


neu4=Toch(neu4)
neu4[neu4=='SOX11/WNT']='Immature'
table(neu4)

neu5=Toch(neu5)
neu5[neu5=='SOX11/WNT']='Immature'

neu5=ordered(as.factor(neu5),levels=c("Mature" ,"Immature","Proliferating"))
table(neu5)

names(neu.pal)[4]="Immature"

neu.palf <- function(n) return(neu.pal)













con1=readRDS('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/Neuron//raw.neuron2_noCycle_new_conos.rds')

con1$plotGraph(groups=neu4)

#con1$findCommunities(resolution=0.65)

#con1$plotGraph()

# c1=con1$clusters$leiden$groups
# nn=names(c1[c1==1])
# neu4[names(neu4) %in% nn]='Neuronal'


#inter=intersect(names(neu4),rownames(con1$embedding))
#DEheatmap(p2,con1$clusters$leiden$groups[inter],'test',)

neu4=readRDS('neu4.cell.ano.rds')
neu4[neu4=="SOX11_WNT" ]="SOX11/WNT" 
neu4[neu4=="Neuronal" ]="Mature" 
neu4[neu4=="noAmp Noradrenergic" ]="NonAmplified" 



neu5=neu4[!neu4 %in% c('NonAmplified')]

neu4=as.factor(neu4)


#neu.pal <- setNames(sample(rainbow(length(levels(neu4)))),levels(neu4));

neu.pal=readRDS('neu.pal.rds')
names(neu.pal)=c("Mature","NonAmplified" ,'Proliferating',"SOX11/WNT" )

neu.palf <- function(n) return(neu.pal)


setwd('/home/meisl/Workplace/neuroblastoma/Figures/Fig3_new')
load('F3.RData')
load('F3new.RData')


a2=con1$plotGraph(groups=neu5,alpha=0.2,size=0.2,plot.na=F,font.size = c(5, 6),palette=neu.palf)
a2
ggsave('F3A.png',a2,height=2.6,width=3.5)




neu5=ordered(as.factor(neu5),levels=c("Mature" ,"SOX11/WNT","Proliferating"))


con2=readRDS('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/Neuron/neuronal_Mar17_all_conos.rds')
con22=readRDS('/home/meisl/Workplace/neuroblastoma/Figures/Fig3/Neuron/neuronal_Mar17_all2_conos.rds')



save(con2,de2,neu5,neu.pal,file='F3.heatmap.RData')




de1=con2$getDifferentialGenes(groups=neu4,z.threshold = 2,append.auc = TRUE)

de2=con2$getDifferentialGenes(groups=neu5,z.threshold = 2,append.auc = TRUE)
names(de2)
names(de1)




pp <- plotDEheatmap(ncdcon,sannot,sannot.de,n.genes.per.cluster = 20 ,show.gene.clusters=T,column.metadata=list(samples=ncdcon$getDatasetPerCell()),order.clusters = T, column.metadata.colors = list(clusters=sannot.pal,samples=sample.pal))


pdf(file='sannot.heatmap.pdf',width=10,height=30); print(pp); dev.off();


source("/home/meisl/bin/conos/R/plot.R")

#options(error=dump.frames())
tmpf=conos::plotDEheatmap(con2,neu5,de2,n.genes.per.cluster = 25 ,show.gene.clusters=T,column.metadata.colors = list(clusters=neu.pal), order.clusters = T, use_raster=F,min.auc = 0.55, row.label.font.size = 3)
tmpf
pdf(file='sannot.small.pdf',width=5,height=4); print(tmpf); dev.off();

png(file='sannot.small.png')
tmpf
dev.off()




de3=p2$getDifferentialGenes(groups=neu5,z.threshold = 2)



allg1=lapply(de1,function(x) {
  x=x[x$Z>0,]
  
  gs=rownames(x[order(x$Z,decreasing=T),])[1:150]
  gs[!is.na(gs)]
})

allg2=lapply(de2,function(x) {
  x=x[x$Z>0,]
  x=x[x$AUC>0.55,]
  
  gs=rownames(x[order(x$AUC,decreasing=T),])[1:150]
  gs[!is.na(gs)]
})


allg3=lapply(de3,function(x) {
  x=x[x$Z>0,]
  gs=rownames(x[order(x$Z,decreasing=T),])[1:150]
  gs[!is.na(gs)]
})





unlist(lapply(allg1,length))

unlist(lapply(allg2,length))

exp=p2$counts

t=DEheatmap2(allg2,exp,neu5,'Neu5',removeGene=NULL,num=25)

DEheatmap2(allg2,exp,neu4,'Neu4',removeGene=NULL,num=25)


allg3$Neuronal=c(allg3$Neuronal,c('PHOX2B','TH','NPY'))


allg4=allg3

tmp=DEheatmap2(allg3,exp,neu5,'FigD',removeGene=NULL,num=23)

allg4[["Neuronal" ]]=c('EEF1A1','EIF4A2','PRPH','DBH','SCG5',"STMN2","STMN4", 'MEG3','RPL10','RBP1','EEF2','GAP43')
allg4[["SOX11" ]]=c("HNRNPH1" , "SOX11" , "MAP1B" ,"SET" , "KIF5A","CTNNB1" , "MARCKS",'WTAP')
allg4[["Proliferating" ]]=c("CDK1","EZH2","FOXM1","TK1",'TOP2A','TMPO','GTSE1','MKI67','NUF2')






tmp=DEheatmap2(allg4,exp,neu5,'FigD',removeGene=NULL,num=23)


aka3 = list('cellType' = neu.pal[as.character(unique(neu5))])
colano <- data.frame('cellType'=tmp$colano$CellType,row.names = rownames(tmp$colano))
head(colano)

rowano=data.frame('cellType'=tmp$rowano$group,row.names = rownames(tmp$rowano))


fout='F3b.png'
x=tmp$x
rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
heat=pheatmap(x,cluster_cols=FALSE,annotation_col = colano,show_colnames = F, #annotation_row = rowano,
              annotation_legend = TRUE, annotation_colors = aka3[1],
              cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =7,width=5.5,height=6,   #3.5*0.02*length(markers),
              breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
heat






fout='F3b2.png'
x=tmp$x
rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
heat=pheatmap(x,cluster_cols=FALSE,annotation_col = colano,show_colnames = F,annotation_row = rowano,
              annotation_legend = TRUE, annotation_colors = aka3[1],
              cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =7,width=5.5,height=6.6,   #3.5*0.02*length(markers),
              breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
heat



texp=tmp$backx

texp=t(texp[colnames(x),rownames(x)])
expression.quantile=0.99
# transform expression values
x1 <- t(apply(as.matrix(texp), 1, function(xp) {
  qs <- quantile(xp,c(1-expression.quantile,expression.quantile))
  xp[xp<qs[1]] <- qs[1]
  xp[xp>qs[2]] <- qs[2]
  xp-min(xp);
  xpr <- diff(range(xp));
  if(xpr>0) xp <- xp/xpr;
  xp
}))


dim(x1)





#gm2=x
gm2=x
require(ComplexHeatmap)
ha = HeatmapAnnotation(
  cells=neu5[colnames(gm2)],
 # annotation_name_side = "bottom",
  #patient=ncdcon$getDatasetPerCell()[colnames(gm2)],
  # specify colors
  col = list(cells=neu.pal
             #pseudotime=circlize::colorRamp2(c(-1, 0, 1), c('darkgreen','grey90','orange'))
  ),
 annotation_legend_param = list(
  # title_position = "lefttop",
   cells = list(nrow = 1)
  
   ),

  border = T
)


col=circlize::colorRamp2(c(-2, 0, 2), c("blue","white","red"))

#col=circlize::colorRamp2(c(0, 0.5, 1), c("blue","white","red"))

hm2 <- Heatmap(gm2,col=col, cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, top_annotation = ha, name='expression', show_heatmap_legend = F, show_row_dend = F, row_names_gp = grid::gpar(fontsize = 8),use_raster=T, raster_device = "CairoPNG")
hm2



pdf(file='F3.B2.pdf',width=4,height=4.9);

draw(hm2, annotation_legend_side = "bottom",column_title=NULL)

dev.off()











gm2=x
require(ComplexHeatmap)
ha = HeatmapAnnotation(
  cells=neu5[colnames(gm2)],
  #patient=ncdcon$getDatasetPerCell()[colnames(gm2)],
  # specify colors
  col = list(cells=neu.pal
             #pseudotime=circlize::colorRamp2(c(-1, 0, 1), c('darkgreen','grey90','orange'))
  ),
  border = T
)

hm2 <- Heatmap(gm2 cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, top_annotation = ha, name='expression', show_heatmap_legend = F, show_row_dend = F, row_names_gp = grid::gpar(fontsize = 8),use_raster=T, raster_device = "CairoPNG")
hm2

pdf(file='F3.heatmap.pdf',width=6,height=7); hm2; dev.off()













allg4 <- list()
allg4[["Neuronal" ]]=c('PHOX2B','TH','NPY','EEF1A1','EIF4A2','PRPH',"STMN2","STMN4", 'MEG3','RTN1','EEF2')
allg4[["SOX11"]]=c("HNRNPH1" ,"CTNNB1" ,    "MARCKS" ,   "SOX11"  ,"MAP1B" , "SET" ,"KIF5A" ,"ATP1B1" , 'WTAP',   "TRA2A")
allg4[["Proliferating" ]]=c('CKS1B','MKI67',"SMC4","TK1","CDK1","TK1" ,"EZH2","FOXM1",'H2AFZ')

tmp=DEheatmap2(allg4,exp,neu5,'Fig3b3',removeGene=NULL,num=30)




gm2=tmp$x
require(ComplexHeatmap)
ha = HeatmapAnnotation(
  cells=neu5[colnames(gm2)],
  #patient=ncdcon$getDatasetPerCell()[colnames(gm2)],
  # specify colors
  col = list(cells=neu.pal
  ),
  border = T
)


hm2 <- Heatmap(gm2, cluster_columns = F, cluster_rows = F,show_column_names = F,border=T, top_annotation = ha, name='expression', show_heatmap_legend = F, show_row_dend = F, row_names_gp = grid::gpar(fontsize = 8),use_raster=T, raster_device = "CairoPNG")
hm2

pdf(file='F3.heatmap32.pdf',width=4,height=5); hm2; dev.off()








###




library(scProcess)

lis=list()

gs=c('ISL1','PNMT','GAP43','PRPH','NTRK1','SOX11','EZH2','MKI67')  #,'GAP43'


for (x in gs){
  dis=(-3.+nchar(x)*0.2)
  if (x %in% c('PHOX2B')){
    dis=-2.8+nchar(x)*0.2
    
  }
  tmp=con1$plotGraph(colors=exp[name2,x],plot.na=F,size=0.1,alpha=0.7)+annotate("text", x = dis, y=-3.45, label = x,size=5.4) +
  scale_x_continuous(limits = range(emb[,1]), expand = c(0, 0)) +
    scale_y_continuous(limits = range(emb[,2]), expand = c(0, 0)) 
  
  lis[[x]]=tmp
}


names(lis)

lrow=3
lcol=3

#lis[['ano']]=con1$plotGraph(groups=neu4,plot.na=F,size=0.1,alpha=0.3)


b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
fout='F3C4.png'
ggsave(fout,b,width = 2.3*lcol,height=1.6*lrow)





lrow=4
lcol=2
b=  cowplot::plot_grid(plotlist=lis, ncol=lcol, nrow=lrow)
fout='F3C6.png'
ggsave(fout,b,width = 2.3*lcol,height=1.6*lrow)



fout='F3C6.pdf'
ggsave(fout,b,width = 2.3*lcol,height=1.6*lrow)





# velocity 


spliced=readRDS('/home/meisl/Workplace/neuroblastoma/velocity/spliced.rds')
unspliced=readRDS('/home/meisl/Workplace/neuroblastoma/velocity/unspliced.rds')

output=function(appname,group4,emb,spliced,unspliced){
  
  obs=data.frame('name'=names(group4),'clusters'=group4)
  cell=names(group4)
  
  splicedf=spliced[,cell]
  unsplicedf=unspliced[,cell]
  gs=rownames(splicedf)
  
  Matrix::writeMM(splicedf,paste(appname,'.splicedf.mtx',sep=''))
  Matrix::writeMM(unsplicedf,paste(appname,'.unsplicedf.mtx',sep=''))
  write.csv(obs,paste(appname,'.obs.csv',sep=''),row.names=F)
  write.csv(emb,paste(appname,'.emb.csv',sep=''),row.names=F,col.names=F)
  write.csv(gs,paste(appname,'.gs.csv',sep=''),row.names=F,col.names=T)
}



samp=con1$getDatasetPerCell()
neu=neu5
name2=intersect(rownames(con1$embedding),names(neu5))

table(samp[name2],neu[name2])

name2=intersect(name2,colnames(spliced))

#table(samp[setdiff(colnames(spliced),name2)])



a2=con1$plotGraph(groups=neu5[name2],alpha=0.3,size=0.25,plot.na=F,font.size = c(4, 5),palette=neu.palf)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#a2=a2+xlim(c(-3.01,3.4))+ylim(c(-4,2.8))
ggsave('F3A3.png',a2,height=2.2,width=3.4)

#a2+theme(panel.border=element_blank())
a2=a2+scale_x_continuous(limits = range(emb[,1]), expand = c(0, 0)) +
  scale_y_continuous(limits = range(emb[,2]), expand = c(0, 0)) 


a2

ggsave('F3A2.png',a2,height=2.2,width=3.4)





name3=intersect(name2,rownames(con1$embedding))
emb=con1$embedding[name2,]

length(name2)
length(name3)
length(neu5)

plot(emb)

range(emb[,1])
range(emb[,2])


emb=con$embedding[names(neu4),]
plot(emb)

emb=emb[(emb[,1]>(-2) & emb[,1] <30),]
name2=rownames(emb)
a2=con$plotGraph(groups=neu4[name2],alpha=0.08,size=0.2,plot.na=F,font.size = c(4, 5),palette=neu.palf)
a2
ggsave('S3a2.png',a2,height=2.3,width=3.3)


x='STMN2'
con$plotGraph(colors=exp[name2,x],plot.na=F,size=0.1,alpha=0.7,title=x) #+annotate("text", x = -3.2+nchar(i)*0.2, y=-3.6, label = x,size=4.2)



con$plotGraph()

library(scProcess)

lis=list()

gs=c('NPY','TH','DBH','STMN2')

lis=lapply(scProcess:::sn(gs),function(x) 
  con$plotGraph(colors=exp[name2,x],plot.na=F,size=0.1,alpha=0.7,title=x) #+annotate("text", x = -3.2+nchar(i)*0.2, y=-3.6, label = x,size=4.2)
)

names(lis)

lrow=2
lcol=2

#lis[['ano']]=con1$plotGraph(groups=neu4,plot.na=F,size=0.1,alpha=0.3)


b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
fout='S3C2.png'
ggsave(fout,b,width = 2.3*lcol,height=1.9*lrow)


















#name2=name2[name2 %in% CNV.cells]

con1$plotGraph(groups=neu[name2])

select=c('NB26', 'NB09' ,'NB12', 'NB13', 'NB24')

select=c('NB26', 'NB09' ,'NB12', 'NB13', 'NB15')

#select=c('NB12','NB13','NB24')
for (i in select){
  tnname=name2[grepl(i,name2)]
  print(table(neu[tnname]))
  output(paste('dat.neuronal.',i,sep=''),neu[tnname],con1$embedding[tnname,],spliced,unspliced)
}

tnname=name2
output(paste('dat.neuronal.','all',sep=''),neu[tnname],con1$embedding[tnname,],spliced,unspliced)






##




library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
cc.genes <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), filters = 'go', values = c('GO:0045202'), mart = ensembl)

ggs=list('synaptic machinery'=as.character(cc.genes[,2]))

length(cc.genes[,2])

cc.genes2 <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), filters = 'go', values = c('GO:0042551'), mart = ensembl)

ggs=list('neuron maturation'=as.character(cc.genes2[,2]))









neu6=Toch(neu4[neu4 %in% c("SOX11/WNT","Neuronal")])

x <- as.matrix(exp[names(neu6),intersect(as.character(cc.genes[,2]),colnames(exp))])

x=t(x)

library(preprocessCore)
exp2=normalize.quantiles(x)
colnames(exp2)=colnames(x)
rownames(exp2)=rownames(x)
x=exp2


score=scale(rowMeans(t(x)))


boxplot(score~neu6,las=2,outline=F)

dat=data.frame('score'=score,'group'=neu6)

p <- ggplot(dat, aes(x=group, y=score)) + 
  geom_boxplot()
p



library(ggpubr)



  tmp=dat
  sig=compare_means(score ~ group,  data = tmp)
  sig
  
 
  sig=sig[sig$p.signif!='ns',]
  sig=sig[order(sig$p.adj),][1:3,]
  sig=sig[!is.na(sig[,1]),]
  siglis=split(sig, seq(nrow(sig)))
  pair=lapply(siglis,function(x) as.character(x[,2:3]))
  pair

  limHeight=1.2
  p1=ggboxplot(dat, x = 'group', y = "score",fill ="group",xlab = "",ylab='synaptic machinery signature',width=0.6,notch = FALSE,
               color = "black",outlier.shape = NA)+
    stat_compare_means(comparisons = pair,label = "p.signif",hide.ns=TRUE,size=3,tip.lengt=0.01,bracket.size =0.3,label.y.npc = "bottom")  #+ # Add pairwise comparisons p-value
 # p1=p1+ylim(c(min(dat$score)*.7,max(dat$score)*limHeight))
  #p1=p1+ geom_point(data = dat,aes(fill = group), size = 0.5, shape = 21, position = position_jitterdodge()) +scale_fill_manual(values=neu.pal)
  p1=p1+ theme(legend.position="none")
  p1=p1+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p1
  
  
  ggsave(paste('synaptic machinery','.boxplot2.pdf',sep=''),p1,height=3.5,width=3)
  
  ggsave(paste('neuron maturation','.boxplot2.pdf',sep=''),p1,height=3.5,width=3)
  
  
  
  names(score)=colnames(x)
  con1$plotGraph(colors=score,plot.na=F)
  
  
  
  table(neu5)
  
  exp[1:4,1:4]
  
  library(scProcess)
  
  gs=intersect(as.character(cc.genes[,2]),colnames(exp))
  cname=names(neu5)
  samp=Toch(con1$getDatasetPerCell())
  
  cname=cname[!grepl('NB16|NB17|NB18|NB22',cname)]
  
  table(samp[cname])
  
  dfTIM=getScore(t(exp),gs,neu5[cname],samp[cname],neu5[cname],magnitude.normal=NULL,rscore=NULL)
  
  drawfig(dfTIM,'neuron maturation') 
  
  drawfig(dfTIM,'synaptic machinery') 
  
  
  
  
###  velocity of  bridge and nearby genes 
  
  
  
  #head(dfM2)
  ylab='synaptic machinery score'
  tmp=dfTIM
  limHeight=1.3
  sig=compare_means(score ~ cell,  data = tmp)
  sig=sig[sig$p.signif!='ns',]
  sig=sig[order(sig$p.adj),][1:3,]
  sig=sig[!is.na(sig[,1]),]
  siglis=split(sig, seq(nrow(sig)))
  pair=lapply(siglis,function(x) as.character(x[,2:3]))
  pair
  
  library(ggpubr)

  p1=ggplot(na.omit(tmp),aes(x=cell,y=score,fill=cell))+geom_boxplot(notch=FALSE,outlier.shape=NA)+theme_classic()
  p1=p1+ geom_point(data = na.omit(tmp),shape=21, size = 1) + #, position = position_jitterdodge()
    stat_compare_means(comparisons = pair,label = "p.signif",hide.ns=TRUE,size=3,tip.lengt=0.01,bracket.size =0.3,label.y.npc = "bottom")  #+ # Add pairwise comparisons p-value
  
  p1=p1+ theme(legend.position="none")+ylab(ylab)
  p1=p1+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')
  p1=p1+ylim(c(min(tmp$score),max(tmp$score)*limHeight))+scale_fill_manual(values=neu.pal)
  p1
  
  ggsave(paste(ylab,'2.pdf',sep=''),p1,height=3.2,width=2)
  
  
  #aes(shape = sample)
  
  
  p1=ggplot(na.omit(tmp),aes(x=cell,y=score,dodge=cell,fill=cell))+geom_boxplot(notch=FALSE,outlier.shape=NA)+theme_classic()
  p1=p1+ geom_point(data = na.omit(tmp),aes(shape = sample), size = 1) + #, position = position_jitterdodge()
    stat_compare_means(comparisons = pair,label = "p.signif",hide.ns=TRUE,size=3,tip.lengt=0.01,bracket.size =0.3,label.y.npc = "bottom")  #+ # Add pairwise comparisons p-value
  
  p1=p1+ theme(legend.position="right")+ylab(ylab)
  p1=p1+theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')
  p1=p1+ylim(c(min(tmp$score),max(tmp$score)*limHeight))+scale_fill_manual(values=neu.pal)
  p1
  
  ggsave(paste(ylab,'3.pdf',sep=''),p1,height=3.2,width=3)
  
  
  
  
  
  
  
  
  
  
  
  tmp=annot[!annot %in% c("Noradrenergic")]
  
  anonew=c(Toch(tmp),Toch(neu4))
  
  table(anonew)
  
  
  samp=Toch(con$getDatasetPerCell())
  cname=names(anonew)
  ano2=data.frame('Cell'=anonew[cname],'Sample'=samp[cname])
  
  # Annotation vs sample
  tmp2 <- acast(ano2, Cell ~ Sample, fun.aggregate=length)
  head(tmp2)
  # Normalise for the number of cells in each library
  tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN='/'))

  
  tab=table(samp[names(neu4)])
  tab=names(tab[tab>100])
  tab
  
  dim(tmp3)
  
  plot(tmp3['Immature',tab],tmp3['Proliferating',tab],xlab='SOX11 ratio',ylab='Proliferatingratio')
  
  
  df <- data.frame(patient=tab,x=tmp3['Immature',tab],y=tmp3['Proliferating',tab])
  
  library(ggrepel)
  
  ct <- cor.test(df$x,df$y,alternative='greater')
  p1 <- ggplot(df, aes(x=x,y=y,label=patient)) + geom_point(aes(color=patient))  + 
    geom_smooth(method='lm',color='gray40',linetype='dashed',alpha=0.15,size=0.5) + theme_bw() + 
    theme(axis.text.y=element_text(angle=90)) + geom_text_repel(aes(color=patient)) + guides(color=F) + 
    xlab("Immature abundance") + ylab("Proliferating abundance") +# xlim(range(df$x))+ylim(c(min(df$y),0.16))+
    geom_text(x=-Inf,y=Inf,label=paste('R=',round(ct$estimate,2),' ','p=',round(ct$p.value,2),sep=''),hjust=-0.05,vjust=1.2,size=3.5)
  
  p1
  
  pdf(file='fig4a.assoc2.pdf',width=2.4,height=2.4); print(p1); dev.off();
  
  
  
  
  
  cor(tmp3['SOX11/WNT',tab],tmp3['Proliferating',tab])
  
  
  pdf1=pheatmap(cor(t(tmp3[,tab]),method='pearson'))
  
  pdf('heatmap.Corr.pdf',height=5,width=6)
  print(pdf1) 
  dev.off()
  
  
  cor(t(tmp3[,tab]),method='spearman')[as.character(unique(neu4)),as.character(unique(neu4))]
  
  
  
  
  scon=readRDS('/home/meisl/Workplace/neuroblastoma/Mar2020/all.con_conos.rds')
  
  
  sn=function(x) { names(x) <- x; return(x); }
  
  
  cm <- scon$getClusterCountMatrices(groups=anonew)
  samples=scon$samples
  
  
  
  ctdm <- lapply(cm,function(x) {apply(x,2,function(x) log10((x/pmax(1,sum(x)))*1e3+1 ) )} )
  
  
  ctdm2 <- lapply(sn(colnames(ctdm[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(ctdm,function(x) x[,ct]))
    tcm
  })
  
  
  
  ctdm2=readRDS('ctdm2.rds')
  
  d1=ctdm2$`SOX11/WNT`
  d2=ctdm2$Proliferating
  
  g1='SOX11'
  g2='MKI67'
  
  #for (g2 in c('MKI67','EZH2','TOP2A','CDK1')){
  
  dat=data.frame('x'=d1[tab,g1],'y'=d2[tab,g2],'patient'=tab)
  cr=round(cor(dat$x,dat$y),3)
  
  p=ggplot(dat, aes(x=x, y=y)) +xlab(g1)+ylab(g2)+
    geom_point(size=2, shape=23)+ggtitle(paste('cor',cr))
    ggtitle(title)
  p
  
  
  
  
  ggsave(paste(g2,'.cor.withSox11.pdf',sep=''),p,height=3,width=3)
  
  cor.test(dat$x,dat$y)
  
}
  


df=dat
ct <- cor.test(df$x,df$y,alternative='greater')
p1 <- ggplot(df, aes(x=x,y=y,label=patient)) + geom_point(aes(color=patient))  + 
  geom_smooth(method='lm',color='gray40',linetype='dashed',alpha=0.15,size=0.5) + theme_bw() + 
  theme(axis.text.y=element_text(angle=90)) + geom_text_repel(aes(color=patient)) + guides(color=F) + 
  xlab("SOX11 exp in Immature") + ylab("MKI67 exp in Proliferating") + 
  geom_text(x=-Inf,y=Inf,label=paste('R=',round(ct$estimate,2),' ','p=',round(ct$p.value,2),sep=''),hjust=-0.05,vjust=1.2,size=3.5)

p1

pdf(file='figExp2.pdf',width=2.4,height=2.4); print(p1); dev.off();




  rexp=rowMeans(as.matrix(exp[names(neu5[neu5=='Proliferating']),]))

  tmp=as.matrix(exp[names(neu5[neu5=='Proliferating']),])
  tlen=apply(tmp, 2, function(x) length(x[x>0]))
  tlen[1:4]
  ratio=tlen/nrow(tmp)
  
  ratio2=ratio[ratio>0.25]


  
  corr=apply(d2[tab,],2,function(x) cor(d1[tab,g1],x))
  
  corr2=corr[names(ratio2)]
  corr2=sort(corr2,decreasing = TRUE)
  corr2[1:12]
  
  
 tt= corr[c(allg4$Proliferating,'PRR11','PHF19')]
 
 pdf('cor.bar.pdf',height=3.4,width=3.7)
 barplot(sort(tt,decreasing=T),las=2,ylab='correlation with SOX11',col=rainbow(12))
 dev.off()
 
  
 
 
  corr['TMEM263']
  ratio[allg4$Proliferating]
  
  
  
  
  tail(sort(corr),n=33)
  
  
  con1$plotGraph(colors=allp2$counts[names(neu4),'EXOC6'],plot.na=F)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  exp=t(allp2$counts)
  
  cname=names(neu4)
  #cname=intersect(cname,colnames(exp))
  ano2=Toch(neu4[cname])
  gs=c('SOX11','MKI67','TK1','FOXM1','NUF2','TOP2A')
  
  tmp=exp[gs,cname]

  
  samp=apply(data.frame(names(ano2)),1, function(x) strsplit(x,'_')[[1]][1])
  
  ttype=paste(ano2,samp,sep='|')
  names(ttype)=names(ano2)
  
  sel= table(ttype) %>% .[.>5] %>% names()
  ttype=ttype[ttype %in% sel]
  
  tmp=tmp[,names(ttype)]
  tm=apply(tmp,1,function(x) tapply(x,ttype,mean))
  
  
  ano=lapply(sn(rownames(tm)),function(x) strsplit(x,'[|]')[[1]])
  sample=unlist(lapply(ano,function(x) x[2]))
  faction=unlist(lapply(ano,function(x) x[1]))
  
  
  index1=( (sample %in% unique(sample)) & faction=="SOX11/WNT")
  tm[index1,'SOX11']

  index2=( (sample %in% unique(sample)) & faction=="Proliferating")
  tm[index2,'MKI67']
  

  plot(tm[index1,'SOX11'],tm[index2,'FOXM1'])
  cor(tm[index1,'SOX11'],tm[index2,'FOXM1'])
  
  
  
  
  
  
  dat=readRDS('/home/meisl/Workplace/neuroblastoma/raw.1210.rds')
  dat1=readRDS('/home/meisl/Workplace/neuroblastoma/raw.rds')
  
  dim(dat[[1]])
  dim(dat1[[1]])
  
  
  dat2=lapply(dat,function(x) t(x))
  aggr2=lapply(dat2[1:7],collapseCellsByType2, groups = neu4, min.cell.count = 0)
  
  aggr2=list()
  for ( i in tab){
    aggr2[[i]]=t(collapseCellsByType2(dat2[[i]],groups = neu4, min.cell.count = 0))
  }
  
  ctdm <- lapply(aggr2,function(x) {apply(x,2,function(x) log10((x/pmax(1,sum(x)))*1e3+1 ) )} )
  
  ctdm2 <- lapply(sn(colnames(ctdm[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(ctdm,function(x) x[,ct]))
    tcm
  })
  
  
  
  
  
  
  
  
  
  
  
  library(dplyr)
  
  sn=function(x) { names(x) <- x; return(x); }
  
  
  collapseCellsByType2=function (cm, groups, min.cell.count)
  {
    #cm=dat2[[1]]
    #groups=neu4
    g1 <- groups[intersect(names(groups), rownames(cm))]
    t1 <- as.numeric(table(g1))
    names(t1) <- levels(g1)
    droplevels <- names(t1)[t1 < min.cell.count]
    g1.n <- names(g1)
    g1 <- as.character(g1)
    names(g1) <- g1.n
    g1[g1 %in% droplevels] <- NA
    g1 <- as.factor(g1)
    aggr <- Matrix.utils::aggregate.Matrix(cm[names(g1),], g1)
    aggr <- aggr[rownames(aggr) != "NA", ]
    return(aggr)
  }
  
  
  
  
getClusterCountMatrices2=function (samples,clustering = NULL, groups = NULL, common.genes = TRUE, 
            omit.na.cells = TRUE) 
  {
    "Estimate per-cluster molecule count matrix by summing up the molecules of each gene for all of the cells in each cluster.\n\n\n       Params:\n\n       - clustering: the name of the clustering that should be used\n       - groups: explicitly provided cell grouping\n\n       - common.genes: bring individual sample matrices to a common gene list\n       - omit.na.cells: if set to FALSE, the resulting matrices will include a first column named 'NA' that will report total molecule counts for all of the cells that were not covered by the provided factor.\n       \n\n       Return: a list of per-sample uniform dense matrices with rows being genes, and columns being clusters\n      "
    if (is.null(groups)) {
      groups <- getClusteringGroups(clusters, clustering)
    }
    groups <- as.factor(groups)
    matl <- lapply(samples, function(s) {
      m <- t(s)
      cl <- factor(groups[match(rownames(m), names(groups))], 
                   levels = levels(groups))
      tc <- conos:::colSumByFactor(m, cl)
      if (omit.na.cells) {
        tc <- tc[-1, , drop = F]
      }
      t(tc)
    })
    if (common.genes) {
      gs <- unique(unlist(lapply(matl, rownames)))
      matl <- lapply(matl, function(m) {
        nm <- matrix(0, nrow = length(gs), ncol = ncol(m))
        colnames(nm) <- colnames(m)
        rownames(nm) <- gs
        mi <- match(rownames(m), gs)
        nm[mi, ] <- m
        nm
      })
    }
    matl
}



cm <- getClusterCountMatrices2(dat1[1:17],groups=anonew)









neu6=neu5

glist=scProcess:::getMarkers()


lis=list()
for (iterm in c("G1.S","S","G2.S","M","M.G1" )){

ggs=glist[[iterm]]
x <- as.matrix(exp[names(neu6),intersect(as.character(ggs),colnames(exp))])

x=t(x)

library(preprocessCore)
exp2=normalize.quantiles(x)
colnames(exp2)=colnames(x)
rownames(exp2)=rownames(x)
x=exp2

score=scale(rowMeans(t(x)))
names(score)=colnames(exp2)

lis[[iterm]]=con1$plotGraph(colors=score,plot.na=F,size=0.1,alpha=0.7,title=iterm)
}


b=  cowplot::plot_grid(plotlist=lis, ncol=2, nrow=3)
fout='cycleScore.png'
ggsave(fout,b,width = 2.3*2,height=2.3*3)








lis=lapply(sn(c('CCNA2','CCNB2','CCNE2','CCND2')), function(x) con1$plotGraph(colors=exp[name2,x],plot.na=F,title=x,size=0.1,alpha=0.7))
b=  cowplot::plot_grid(plotlist=lis, ncol=2, nrow=2)
fout='cycle.png'
ggsave(fout,b,width = 2.3*2,height=2.1*2)





