
library(scProcess)

setwd('/home/meisl/Workplace/neuroblastoma/Figures/Fig2')

con=readRDS('/home/meisl/Workplace/neuroblastoma/ncdcon.rds')

con=Conos$new(con)
con$plotGraph(gene='MYCN',size=0.3,alpha=0.9,title='MYCN')






annot <- readRDS("/d0-mendel/home/meisl/Workplace/neuroblastoma/cell.annotation.Jan2020.rds")$cellano
annot <- setNames(as.character(annot),names(annot))
annot[annot=='unknown'] <- 'Th' # Ca2+ mediated apoptosis, involivng LCK, CD3D, RAC2, etc. Lumping in with Th 
annot[annot=='Plasmacytoid'] <- 'pDC'
annot[grep('Cytotoxic',annot)] <- 'Tcyto'
annot[grep('T helper',annot)] <- 'Th'
annot[grep('NK',annot)] <- 'NK'
annot[grep('Bcell',annot)] <- 'B'
annot[grep('PlasmaCell',annot)] <- 'Plasma'
annot[grep('Mast',annot)] <- 'Mast'
annot[grep("Bridge",annot)] <- 'SCP-like'


dannot <- as.factor(annot)

annot[grep("SOX11|Stress|euronal|Prolif",annot)] <- 'Noradrenergic'
annot[grep("Bridge",annot)] <- 'SCP-like'
annot <- as.factor(annot);




# Load ggplot2
library(ggplot2)
sample=readRDS("/d0-mendel/home/meisl/Workplace/neuroblastoma/cell.annotation.Jan2020.rds")$sample
sample=sample[names(annot)]


# pie chart for cell frequnce among samles 

set.seed(3)
annot.pal <- setNames(sample(rainbow(length(levels(annot)),v=0.95,s=0.7)),levels(annot));
annot.palf <- function(n) return(annot.pal)
dannot.pal <- setNames(sample(rainbow(length(levels(dannot)),v=0.85,s=0.9)),levels(dannot));
x <- which(levels(dannot) %in% levels(annot)); dannot.pal[x] <- annot.pal[match(levels(dannot)[x],levels(annot))]
dannot.palf <- function(n) return(dannot.pal)



# for samples
sample.pal <- setNames(sample(rainbow(length(con$samples),v=0.7,s=0.9)),names(con$samples))
sample.pal <- sample.pal[c(sort(names(sample.pal))[-1],'Adr')]
sample.palf <- function(n) return(sample.pal[1:n])
sample.palf <- function(n) { browser(); return(sample.pal[1:n])}




##

annot=Toch(annot)
annot[annot=='Noradrenergic'] <- 'Adrenergic'

names(annot.pal)[11]='Adrenergic'








cname=names(annot)

# expression matrix 
exp=mergeDat2(lapply(con$samples,function(x) t(x$counts)))



gs=c('TH','ACTA2','DCN','VWF','PLVAP','CD3D','CD4','CD8A','XCL1','IGHG1','MS4A1','TPSAB1','LILRA4','LYZ',
     'S100B','PDGFRB')       

length(gs)

lcol=4
lrow=4
lis=list()
for (gene in gs){
  print(gene)
  alpha =0.2
  if (gene %in% c('IL10','IL1B')){
    alpha =0.5
  }
 # a=con$plotGraph(colors =exp[gene,cname],title=gene,alpha =alpha,plot.na=F,size=0.2)+ theme(plot.title = element_text(size=19,hjust = 0.5))
  a=con$plotGraph(colors =exp[gene,cname],alpha =alpha,plot.na=F,size=0.1)+annotate("text", x = -33+nchar(gene), y=46, label = gene,size=3.9)
  lis[[gene]]=a
}
b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
ggsave('S1.markerGene.pdf',b,width = 2*lrow,height=2*lcol)




###
annot=as.factor(annot)
sannot <- setNames(as.character(annot),names(annot));
sannot[grep("^T|ILC3|NK",sannot)] <- 'T cells'
sannot[grep("^B$|Plasma|pDC",sannot)] <- 'B cells'
sannot[grep("Macrop|Mono|mDC|Mast",sannot)] <- 'Myeloid'
sannot <- as.factor(sannot)
sannot=ordered(sannot,levels=c("Adrenergic"  ,   "B cells"   , "Pericytes"    ,    "Endothelial" ,  "Myeloid"  ,  "Mesenchymal"     ,    "Myofibroblasts",
                                 "SCP-like"   ,    "T cells" ))

sannot.pal <- setNames(annot.pal[match(levels(sannot),names(annot.pal))],levels(sannot))
sannot.pal['T cells'] <- annot.pal['Tcyto']
sannot.pal['B cells'] <- annot.pal['B']
sannot.pal['Myeloid'] <- annot.pal['Macrophages']




library(reshape2)

ano2=data.frame('Cell'=sannot[cname],'SampleType'=sample[cname])

# Annotation vs sample
tmp2 <- acast(ano2, Cell ~ SampleType, fun.aggregate=length)
head(tmp2)
# Normalise for the number of cells in each library
tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN='/'))
tmp4 <- melt(tmp3)
head(tmp4)
names(tmp4) <- c('cell', 'sample','pc.of.sample')
head(tmp4)


p=ggplot(tmp4, aes(x=sample, fill=cell, y = pc.of.sample)) +theme_bw()+
  geom_bar(stat='identity', position='fill') + 
  scale_fill_manual(values=sannot.pal)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y=element_blank(),panel.border = element_blank(),
        axis.title.x=element_blank())

p  
  
ggsave(paste('F1','.porpotion2.png',sep=''),plot=p,height=4.2,width=4.8)


ggsave(paste('F1','.porpotion2.pdf',sep=''),plot=p,height=4.2,width=4.8)





# pie chart
ratio=(table(sample)/length(sample))*100

# Create Data
data <- data.frame(
  cell=paste(names(ratio),': ',round(ratio,2),'%',sep=''),
  value=as.numeric(ratio)
)

# Basic piechart
a1=ggplot(data, aes(x="", y=value, fill=cell)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=as.character(sample.pal[names(ratio)]))+
  theme_void()
a1



a1=ggplot(data, aes(x="", y=value, fill=cell)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=as.character(annot.palf(length(unique(sample.palf)))))+
  theme_void()



ggsave('S1a.pdf',a1,height=5,width=5)




nname <- readPagoda2SelectionFile('/home/meisl/Workplace/neuroblastoma/Figures/Fig2/bridge2')[['bridge2']]$cells

nname <- readPagoda2SelectionFile('NBF')[['NBF']]$cells


nname <- readPagoda2SelectionFile('bridge.flank5')[['bridge.flank5']]$cells



#nname <- readPagoda2SelectionFile('bridge.flank')[['bridge.flank']]$cells


emb=con$embedding
emb=emb[(emb[,1]>(-13) & emb[,1] <(10)),]
nname=rownames(emb)

plot(emb)

bridge.ano=Toch(annot[nname])
bridge.ano=bridge.ano[bridge.ano %in% c('Mesenchymal','SCP-like','Noradrenergic')]

#bridge.ano=bridge.ano[bridge.ano %in% c('SCP-like')]

nname=names(bridge.ano)

emb=emb[nname,]
plot(emb)
abline(h=(-19))

emb=emb[emb[,2]<(-19),]


nname=rownames(emb)

dannot.palf2 <- function(n) return(annot.pal[unique(Toch(bridge.ano))])

bridge.ano=bridge.ano[nname]
a2=con$plotGraph(groups=bridge.ano,alpha=0.3,size=0.5,plot.na=F,font.size = c(5, 6),palette=dannot.palf2)
a2

ggsave('F2a.pdf',a2,height=3.5,width=3.5)




gs=c('TH','PHOX2B','NPY','ERBB3','PLP1','SOX10','PRRX1','DCN','COL1A1')

lcol=3
lrow=3
lis=list()
for (gene in gs){
  print(gene)
  alpha =0.6
  if (gene %in% c('IL10','IL1B')){
    alpha =0.5
  }
  # a=con$plotGraph(colors =exp[gene,cname],title=gene,alpha =alpha,plot.na=F,size=0.2)+ theme(plot.title = element_text(size=19,hjust = 0.5))
  a=con$plotGraph(colors =exp[gene,nname],alpha =alpha,plot.na=F,size=0.3)+annotate("text", x = -11+nchar(gene)*0.5, y=-20.9, label = gene,size=3.9)
  lis[[gene]]=a
}
b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
ggsave('S2a.png',b,width = 2*lrow,height=2*lcol)

plot(con$embedding[nname,])




nname=names(annot[annot %in% c('SCP-like')])

tab=sort(table(sample[nname]),decreasing=T)
pdf('barplot.pdf',height=3.7,width=4.3)
barplot(tab,las=2,col=sample.pal,ylab='Frequency of Bridge cells')
dev.off()







cname=names(annot[annot %in% c('SCP-like','Mesenchymal','Noradrenergic')])
ano2=data.frame('Cell'=Toch(annot)[cname],'Sample'=sample[cname])

table(annot[cname])


# Annotation vs sample
tmp2 <- acast(ano2, Cell ~ Sample, fun.aggregate=length)
head(tmp2)
# Normalise for the number of cells in each library
tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN='/'))

tab=tmp3['SCP-like',]

tab=sort(tab,decreasing=T)
tab





# barplot with error bar 

dat=data.frame('SCP'=tmp2['SCP-like',],'all'=colSums(tmp2),'sample'=colnames(tmp2))
dat$ratio=dat$SCP/dat$all


getSE=function(p_hat,n){
  z=1.96
  se=z*sqrt(p_hat*(1-p_hat)/n)
  return(se)
}


dat$se=apply(dat,1,function(x) getSE(as.numeric(x[4]),as.numeric(x[2])))

ratio=dat$ratio 
names(ratio)=dat$sample
levels=names(sort(ratio,decreasing = TRUE))

dat$sample=ordered(as.factor(dat$sample),levels=levels)


# Use 95% confidence intervals instead of SEM
p=ggplot(dat, aes(x=sample, y=ratio, fill=sample))   +theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")  + ylab('Fraction of Bridge cells')+
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))
p=p+scale_fill_manual(values=sample.pal)
p=p+theme(legend.position="none",legend.title = element_blank(),legend.key.size = unit(0.8, "lines"),
          legend.text = element_text(size=8))

p=p+geom_text(data = dat, aes(x = sample, y = ratio+se + 0.03, label = SCP),size=2.9)

p

ggsave(paste('S1','.boxplot4.pdf',sep=''),p,height=2.7,width=3.6)







###

load('/home/meisl/Workplace/neuroblastoma/Figures/CNV_fig/S2/FDR5.RData')

load('/home/meisl/Workplace/neuroblastoma/Figures/CNV/CNV.0211.RData')

amplif=Toch(annot)
amplif[ names(amplif) %in% all.CNV2]='Amplified'
amplif[ !(names(amplif) %in% all.CNV2)]='NonAmplified'



cname=names(annot[annot %in% c('SCP-like')])
ano2=data.frame('Cell'=amplif[cname],'Sample'=sample[cname])

# Annotation vs sample
tmp2 <- acast(ano2, Cell ~ Sample, fun.aggregate=length)
tmp2

dat=data.frame('SCP'=tmp2['Amplified',],'all'=colSums(tmp2),'sample'=colnames(tmp2))
dat$ratio=dat$SCP/dat$all












#
# velocity  (prepare data )


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

# NB12 NB13 NB15 NB17 NB18 NB20 NB23 NB26

nname=intersect(colnames(spliced),nname)

table(sample[nname],Toch(bridge.ano[nname]))




select=c('NB12', 'NB13' ,'NB15', 'NB17', 'NB18', 'NB20' ,'NB23', 'NB26')


lis=list()
select=c('NB26','NB13','NB12','NB23','NB20')
for (i in select){
  tnname=nname[grepl(i,nname)]
  output(paste('bride8.',i,sep=''),bridge.ano[tnname],con$embedding[tnname,],spliced,unspliced)
  #lis[[i]]=con$plotGraph(groups=bridge.ano[tnname],title=i,alpha=0.3,size=0.5,plot.na=F,font.size = c(5, 6),palette=dannot.palf2)
}


a2=con$plotGraph(groups=bridge.ano,alpha=0.3,size=0.5,plot.na=F,font.size = c(5, 6),palette=dannot.palf2)
a2


lrow=2
lrow=3

b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
fout='F2C.png'
ggsave(fout,b,width = 2.3*lcol,height=2.3*lrow)





#tmp=bridge.ano[nname]
#nname1=names(tmp[tmp %in% c('Bridge','Neural')])

output('bride.all8',bridge.ano[nname],con$embedding[nname,],spliced,unspliced)






cname=names(annot)
cname=intersect(colnames(spliced),cname)


con$plotGraph(groups=annot[cname])

saveRDS(con$embedding[nname,],'bridge.flank.embedding.rds')
saveRDS(con$embedding[cname,],'all.embedding.rds')
saveRDS(annot,'all.annot.rds')







###

untreat=c('NB01', 'NB02', 'NB12', 'NB15', 'NB16', 'NB17', 'NB19', 'NB21' ,'NB26')

anoType=sample
index=anoType %in% untreat
anoType[index]='untreat'
anoType[!index]='treat'

table(anoType)




##

bridge=Toch(annot[annot %in% c('SCP-like')])
nname=names(bridge)

a2=con$plotGraph(groups=bridge,alpha=0.3,size=0.5,plot.na=F,font.size = c(5, 6),palette=dannot.palf2)
a2


load('/home/meisl/Workplace/neuroblastoma/Figures/CNV/CNV.0211.RData')

amplif=Toch(annot)
amplif[ names(amplif) %in% CNV.cells]='Amplified'
amplif[ !(names(amplif) %in% CNV.cells)]='NonAmplified'


con$plotGraph(colors=allScoreA[nname],plot.na=F,size=0.8,alpha=0.8)


#


nname1=nname[grepl('NB09',nname)]

test=amplif[nname1]
table(test)

DEheatmap(p2,test,'Bridge',removeGene=NULL,num=20)


con$plotGraph(gene='KRT19')

exp=t(p2$counts)


gene='KRT19'
print(gene)

# a=con$plotGraph(colors =exp[gene,cname],title=gene,alpha =alpha,plot.na=F,size=0.2)+ theme(plot.title = element_text(size=19,hjust = 0.5))
a=con$plotGraph(colors =exp[gene,nname],alpha =alpha,plot.na=F,size=0.5) #+annotate("text", x = -13.3+nchar(gene)*0.5, y=-20.9, label = gene,size=3.9)
a





f1b=con$plotGraph(groups=amplif[nname],plot.na=F,size=0.2,alpha=0.2,mark.groups=F,
                   show.legend=T,legend.pos=c(0.98, 0.98))+
  guides(colour = guide_legend(override.aes = list(size=6, alpha = 1)))
f1b

ggsave('Bridge.cell.amplification.pdf',f1b,width = 3,height=3)


sample.pal <- setNames(sample(rainbow(length(unique(sample[nname])),v=0.7,s=0.9)),unique(sample[nname]))
sample.palf <- function(n) return(sample.pal[1:n])

length(unique(sample[nname]))

p1 <- con$plotGraph(alpha=0.3,size=0.5,groups=sample[nname],plot.na=F,mark.groups=F,palette=sample.palf)+
  theme(legend.position="bottom", legend.box = "horizontal")+
  guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
p1
##

ggsave('Bridge.sample.pdf',p1,width = 4,height=5.1)


tab=table(sample[nname])







raw=readRDS('~/Workplace/neuroblastoma/raw.1210.rds')


exp=mergeDat2(raw)

#exp=mergeDat2(lapply(con$samples,function(x) t(x$counts)))



bridge=Toch(annot[annot %in% c('SCP-like','Noradrenergic','Mesenchymal','Pericytes','Myofibroblasts')])

bridge=Toch(annot[annot %in% c('Noradrenergic')])

nname=names(bridge)



amplif=Toch(annot)
amplif[ names(amplif) %in% CNV.cells]='Amplified'
amplif[ !(names(amplif) %in% CNV.cells)]='NonAmplified'

bridge.amp=amplif[nname]




nname=nname[grepl('NB09',nname)]



group=paste(bridge,bridge.amp)
names(group)=names(bridge)

group=group[nname]
unique(group)



nname=names(annot)
nname=nname[grepl('NB12',nname)]
group=Toch(sannot[nname])
table(group)


dat=as.matrix(exp[,nname])

ano=data.frame(nname,group)
table(ano[,2])
table(group)

dim(ano)
dim(dat)

table(colnames(dat)==ano[,1])

write.table(ano,'inferCNA.bridge.NB12.ano.txt',sep='\t',col.names=F,row.names=F,quote=F)
write.table(dat,'inferCNA.bridge.NB12.exp.txt',sep='\t',col.names=T,row.names=T,quote=F)

table(is.na(dat))


library("infercnv")



devtools::install_github("broadinstitute/infercnv")
devtools::install_github("broadinstitute/infercnv", ref="RELEASE_3_10")

fano='inferCNA.bridge.NB12.ano.txt'
fexp='inferCNA.bridge.NB12.exp.txt'
cname='bridge.NB12_new3'



infercnv_obj = CreateInfercnvObject(raw_counts_matrix=fexp,
                                    annotations_file=fano,
                                    delim="\t",
                                    gene_order_file="/home/meisl/tools/inferCNV/tests/hg19.sort.inferCNV.txt",
                                    ref_group_names=c('B cells',"Endothelial",'Myeloid','T cells')) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=cname, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             mask_nonDE_genes = T,
                             HMM=TRUE)




setwd('/home/meisl/Workplace/neuroblastoma/aCGH/NB09')








# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=fexp,
                                    annotations_file=fano,
                                    delim="\t",
                                    gene_order_file="/home/meisl/tools/inferCNV/tests/hg19.sort.inferCNV.txt",   #,'Epitheial_Hillock'
                                    ref_group_names=c('B cells',"Endothelial",'Myeloid','T cells'))






out_dir=cname
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             use_zscores = TRUE,
                             mask_nonDE_genes = T,
                             include.spike=F  # used for final scaling to fit range (0,2) centered at 1.
)








ncdcon <- readRDS("~pkharchenko/m/ninib/NB/ncdcon.rds")

sannot <- setNames(as.character(annot),names(annot));
sannot[grep("^T|ILC3|NK",sannot)] <- 'T cells'
sannot[grep("^B$|Plasma|pDC",sannot)] <- 'B cells'
sannot[grep("Macrop|Mono|mDC|Mast",sannot)] <- 'Myeloid'
sannot <- as.factor(sannot)

sannot.pal <- setNames(annot.pal[match(levels(sannot),names(annot.pal))],levels(sannot))
sannot.pal['T cells'] <- annot.pal['Tcyto']
sannot.pal['B cells'] <- annot.pal['B']
sannot.pal['Myeloid'] <- annot.pal['Macrophages']


ncdcon=Conos$new(ncdcon)
sannot.de <- ncdcon$getDifferentialGenes(groups=sannot,n.cores=30,append.auc=TRUE,z.threshold=0,upregulated.only=T)

source("/home/meisl/bin/conos/R/plot.R")

#options(error=dump.frames())
plotDEheatmap(ncdcon,sannot,sannot.de,n.genes.per.cluster = 20 ,show.gene.clusters=T,column.metadata=list(samples=ncdcon$getDatasetPerCell()), column.metadata.colors = list(clusters=sannot.pal,samples=sample.pal), order.clusters = T, use_raster=F,min.auc = 0.75, row.label.font.size = 3)


traceback()



tf <- droplevels(sannot[grep("Noradr|SCP|Mesench",sannot)]); tf <-factor(tf,levels=levels(tf)[c(2,3,1)])
pp <- plotDEheatmap(ncdcon,tf,sannot.de,n.genes.per.cluster = 10 ,show.gene.clusters=T,column.metadata=list(samples=ncdcon$getDatasetPerCell()),order.clusters = F, column.metadata.colors = list(clusters=sannot.pal,samples=sample.pal),additional.genes = c("SOX10"),min.auc=0.75,row.label.font.size = 10,split=T, split.gap=1,column_title = NULL)
pp

debugger()


browser()

#devtools::install_github('hms-dbmi/conos')



genes <- c("CD79B",'MS4A1','PLVAP','EPCAM','SPINK2','IL7R','LYZ','C1QB','CLEC10A','COL1A1','LUM','S100A9','IL1B','MYL9','GNLY','RTN1','MDK','NNAT','IRF4','MYH11','FCRL5','PLP1','GZMA','RORA','TIGIT')
genes <- c("CD3D",'CD7','CD79A','HLA-DRA','CD74','LYZ','STMN2','MDK','MYL9','MYH11','PECAM1','LUM','PLP1')
tf <- sannot[unlist(tapply(1:length(sannot),sannot,function(x) { if(length(x)>1e2) x <- sample(x,1e2); x}))]
pp <- plotDEheatmap(ncdcon,tf,sannot.de,n.genes.per.cluster = 20 ,show.gene.clusters=T,column.metadata=list(samples=ncdcon$getDatasetPerCell()), column.metadata.colors = list(clusters=sannot.pal,samples=sample.pal), order.clusters = T, use_raster=T, additional.genes = genes, labeled.gene.subset = genes, min.auc = 0.75)
pp







source("/home/pkharchenko/m/p2/conos/R/plot.R")
genes <- c("CD79B",'MS4A1','PLVAP','EPCAM','SPINK2','IL7R','LYZ','C1QB','CLEC10A','COL1A1','LUM','S100A9','IL1B','MYL9','GNLY','RTN1','MDK','NNAT','IRF4','MYH11','FCRL5','PLP1','GZMA','RORA','TIGIT')
genes <- c("CD3D",'CD7','CD79A','HLA-DRA','CD74','LYZ','STMN2','MDK','MYL9','MYH11','PECAM1','LUM','PLP1')
tf <- sannot[unlist(tapply(1:length(sannot),sannot,function(x) { if(length(x)>1e2) x <- sample(x,1e2); x}))]
pp <- plotDEheatmap(ncdcon,tf,sannot.de,n.genes.per.cluster = 20 ,show.gene.clusters=T,column.metadata=list(samples=ncdcon$getDatasetPerCell()), column.metadata.colors = list(clusters=sannot.pal,samples=sample.pal), order.clusters = T, use_raster=T, additional.genes = genes, labeled.gene.subset = genes, min.auc = 0.75)
pp



layout(mm, widths = mw, heights = mh)
