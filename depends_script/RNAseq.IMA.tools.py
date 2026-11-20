#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse
import click
import pandas as pd
from multiprocessing.pool import Pool


'''
@author: yueyao
@time: 2019/6/21
@file: RNAseq.IMA.tools.py
@mail: yueyaomail@gmail.com

'''

@click.group()
def main():
    pass

def run_r(cmd):
    os.system(cmd)

@click.command()
@click.option("-i","--input", help="all.genes.fpkm.xls file create by DEG step.")
@click.option("-g","--group", help="Group List create by DEG step.")
@click.option("-c","--compare", help="Compare List create by DEG step.")
@click.option("-o","--outdir", help="output directory.")
def cluster(input,group,compare,outdir):
    '''
    do cluster for deg result.
    '''
    if input is None or group is None or compare is None or outdir is None:
        logging.warning("\nRNAseq.IMA.tools.py cluster -i ../4.Comparison/all.genes.fpkm.xls -g ../4.Comparison/temp/GroupList.txt -c ../4.Comparison/temp/CompareList.txt -o .")
        sys.exit(1)

    groupdic={}
    agroup=[]
    asample=[]
    asorted=[]
    cmd=[]
    sampledic={}
    with open(group,'r') as f:
        for line in f.readlines():
            arr=line.strip().split()
            groupdic.setdefault(arr[1],arr[0])
            for ss in arr[1].split(","):
                sampledic.setdefault(ss,arr[0])
    with open(compare,'r') as f:
        for line in f.readlines():
            com=line.strip().split()
            control=groupdic.get(com[0])
            treat=groupdic.get(com[1])
            conl=len(com[0].split(","))*[control]
            treal=len(com[1].split(","))*[treat]
            gg=com[0].split(",")+com[1].split(",")
            ll=conl+treal
            asample=asample+gg
            asorted=asorted+[control,treat]
            gs=','.join('"' + item + '"' for item in gg)
            ns=','.join('"' + item + '"' for item in ll)
            sorted=','.join('"' + item + '"' for item in [control,treat])
            outfpkm=outdir+"/"+control+"_vs_"+treat+".fpkm.txt"
            RCODE='''
library(gmodels)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(ggrepel)

expr <- read.table("%(input)s",sep="\\t",header=T,check.names=F,quote="")

expr <- subset(expr, select = c("gene",%(samples)s))

keep <- apply(expr[,-1],1,sum) > 0
write.table(expr[keep,],file="%(outfpkm)s",sep="\\t",quote=F,row.names=F)
expr <- expr[,-1]

group <- c(%(group)s)
name <- c(%(samples)s)
c=data.frame(group,name,row.names = name)

data = t(as.matrix(log2(expr+1)))
data.pca = fast.prcomp(data,retx=T,scale=F,center=T)

a = summary(data.pca)
tmp = a[4]$importance

pro1 = as.numeric(sprintf("%(fh)s.3f",tmp[2,1]))*100
pro2 = as.numeric(sprintf("%(fh)s.3f",tmp[2,2]))*100

pc1 = as.data.frame(a$x)
pc=merge(pc1,c,by = "row.names",all=F)
names=pc[,1]
pc=pc[,-1]
row.names(pc)=names
sorted=c(%(sorted)s)

xlab = paste("PC1(",pro1,"%(fh)s)",sep="")
ylab = paste("PC2(",pro2,"%(fh)s)",sep="")
pca_nolab = ggplot(pc,aes(PC1,PC2)) +
geom_point(size=3,aes(color=group)) +
labs(x=xlab,y=ylab) +
geom_hline(yintercept=0,linetype="dashed",color="black") +
geom_vline(xintercept=0,linetype="dashed",color="black") +
theme_bw()+
  theme(axis.title.x=element_text(colour = "black",size=6),
        axis.title.y=element_text(colour = "black",size=6),
        panel.grid=element_blank(),
        axis.line=element_line(size=0.25,colour="black"),
        axis.text = element_text(colour = "black",size=6)
        ) +
scale_color_discrete(limits=sorted,guide=FALSE) +
labs(color="Sample") +stat_ellipse(geom="polygon",alpha=1/4,aes(fill=group)) +
scale_fill_discrete(limits=sorted,guide=FALSE) +geom_rug(aes(color=group))


#### change by yueyao at 2019.11.19
pca = pca_nolab+geom_text_repel(aes(label=names),size=2) 
#### change by yueyao at 2019.11.19


out1 = paste(sorted,collapse="_vs_")
out2 = paste("%(outdir)s","/PCA.",out1,sep="")
ggsave(paste(out2,".pdf",sep=""),pca,width=3.5,height=3.5)
ggsave(paste(out2,".png",sep=""),pca,width=3.5,height=3.5)

####change by yueyao at 2019.11.19
ggsave(paste(out2,".nolab.pdf",sep=""),pca_nolab,width=3.5,height=3.5)
ggsave(paste(out2,".nolab.png",sep=""),pca_nolab,width=3.5,height=3.5)
####change by yueyao at 2019.11.19

###################################################################################################
hc = hclust(dist(data),method="average")

color = rep(0,length(labels(hc)))

for(m in 1:length(c$name)){
    for(n in 1:length(sorted)){
      if(c$group[m] == sorted[n]){
         colnum = grep(c$name[m],labels(hc))
        color[colnum] = rainbow(length(sorted))[n]
      }
   }
}

dend = as.dendrogram(hc) %(ss)s set("labels_col",color) %(ss)s set("labels_cex",1) %(ss)s set("leaves_pch",19) %(ss)s set("leaves_cex",1) %(ss)s set("leaves_col",color)

out3 = paste("%(outdir)s","/hclust",".",out1,sep="")
rw = length(labels(hc))/2.5
pdf(file=paste(out3,".pdf",sep=""),width=rw,height=4)
par(mar=c(10,3,2,0))
plot(dend,ylab="Height")
dev.off()

rw = rw*80
png(filename=paste(out3,".png",sep=""),width=rw,height=400)
par(mar=c(10,3,2,0))
plot(dend,ylab="Height")
dev.off()

'''
            RCODE=RCODE % {'input':input,'samples':gs,'outfpkm':outfpkm,'group':ns,'sorted':sorted,'fh':"%",'ss':"%>%",'outdir':outdir}
            with open(outdir+'/'+control+"_vs_"+treat+".cluster.r",'w') as fw:
                fw.write(RCODE)
            rcmd="Rscript "+outdir+'/'+control+"_vs_"+treat+".cluster.r"
            cmd.append(rcmd)
        ALLRCODE='''
library(gmodels)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(ggrepel)

expr <- read.table("%(input)s",sep="\\t",header=T,check.names=F,quote="")
expr <- expr[,c(-1,-2)]

group <- c(%(group)s)
name <- c(%(samples)s)
c=data.frame(group,name,row.names = name)

data = t(as.matrix(log2(expr+1)))
data.pca = fast.prcomp(data,retx=T,scale=F,center=T)

a = summary(data.pca)
tmp = a[4]$importance

pro1 = as.numeric(sprintf("%(fh)s.3f",tmp[2,1]))*100
pro2 = as.numeric(sprintf("%(fh)s.3f",tmp[2,2]))*100

pc1 = as.data.frame(a$x)
pc=merge(pc1,c,by = "row.names",all=F)
names=pc[,1]
pc=pc[,-1]
row.names(pc)=names
sorted=c(%(sorted)s)

xlab = paste("PC1(",pro1,"%(fh)s)",sep="")
ylab = paste("PC2(",pro2,"%(fh)s)",sep="")
pca_nolab = ggplot(pc,aes(PC1,PC2)) +
geom_point(size=3,aes(color=group)) +
labs(x=xlab,y=ylab) +
geom_hline(yintercept=0,linetype="dashed",color="black") +
geom_vline(xintercept=0,linetype="dashed",color="black") +
theme_bw()+
  theme(axis.title.x=element_text(colour = "black",size=6),
        axis.title.y=element_text(colour = "black",size=6),
        panel.grid=element_blank(),
        axis.line=element_line(size=0.25,colour="black"),
        axis.text = element_text(colour = "black",size=6)
        ) +
scale_color_discrete(limits=sorted,guide=FALSE) +
labs(color="Sample") +stat_ellipse(geom="polygon",alpha=1/4,aes(fill=group)) +
scale_fill_discrete(limits=sorted,guide=FALSE) +geom_rug(aes(color=group))


#### change by yueyao at 2019.11.19
pca = pca_nolab+geom_text_repel(aes(label=names),size=2) 
#### change by yueyao at 2019.11.19

out1 = "allsamples"
out2 = paste("%(outdir)s","/PCA.",out1,sep="")
ggsave(paste(out2,".pdf",sep=""),pca,width=3.5,height=3.5)
ggsave(paste(out2,".png",sep=""),pca,width=3.5,height=3.5)


####change by yueyao at 2019.11.19
ggsave(paste(out2,".nolab.pdf",sep=""),pca_nolab,width=3.5,height=3.5)
ggsave(paste(out2,".nolab.png",sep=""),pca_nolab,width=3.5,height=3.5)
####change by yueyao at 2019.11.19


###################################################################################################
hc = hclust(dist(data),method="average")

color = rep(0,length(labels(hc)))

for(m in 1:length(c$name)){
    for(n in 1:length(sorted)){
      if(c$group[m] == sorted[n]){
         colnum = grep(c$name[m],labels(hc))
        color[colnum] = rainbow(length(sorted))[n]
      }
   }
}

dend = as.dendrogram(hc) %(ss)s set("labels_col",color) %(ss)s set("leaves_pch",19) %(ss)s set("leaves_cex",2) %(ss)s set("leaves_col",color)

out3 = paste("%(outdir)s","/hclust",".",out1,sep="")
rw = length(labels(hc))/2.5
pdf(file=paste(out3,".pdf",sep=""),width=rw,height=4)
par(mar=c(10,3,2,0))
plot(dend,ylab="Height")
dev.off()

rw = rw*80
png(filename=paste(out3,".png",sep=""),width=rw,height=400)
par(mar=c(10,3,2,0))
plot(dend,ylab="Height")
dev.off()

'''
        asample=set(asample)
        for js in asample:
            agroup.append(sampledic.get(js))
        asorted=set(asorted)
        ags=','.join('"' + item + '"' for item in asample)
        ans=','.join('"' + item + '"' for item in agroup)
        assorted=','.join('"' + item + '"' for item in asorted)
        ALLRCODE=ALLRCODE % {'input':input,'samples':ags,'group':ans,'sorted':assorted,'fh':"%",'ss':"%>%",'outdir':outdir}
        with open(outdir+"/Allsamples.r",'w') as fa:
            fa.write(ALLRCODE)
        rcmd = "Rscript " + outdir+"/Allsamples.r"
        cmd.append(rcmd)
        p = Pool(3)
        for s in cmd:
            p.apply_async(run_r, args=(s,))
        print("Waiting for all subprocess done...")
        p.close()
        p.join()
        print ("All subprocess done")

@click.command()
@click.option("-i","--input", help="A xx.psi file create by suppa.")
@click.option("-p","--prefix", help="output prefix.")
def statas(input,prefix):
    '''
    stat diffsplice result.
    '''
    ase = ['SE', 'A3', 'MX', 'A5', 'RI', 'AF', 'AL']
    for i in ase:
        out=prefix+"."+i+".xls"
        fw=open(out,'w')
        with open(input,'r') as f:
            for line in f.readlines()[1:]:
                line=line.strip().split()
                geneid=line[0]
                if geneid.split(':')[-1] == "-" and geneid.split(':')[0].endswith(i):
                    start1=int(geneid.split(':')[2].split('-')[0])
                    start2=int(geneid.split(':')[3].split('-')[0])
                    exonlen=start1-start2+1
                    fw.write(i+'\t'+str(exonlen)+"\n")
                elif geneid.split(':')[-1] == "+" and geneid.split(':')[0].endswith(i):
                    start1=int(geneid.split(':')[2].split('-')[1])
                    start2=int(geneid.split(':')[3].split('-')[1])
                    exonlen=start2-start1+1
                    fw.write(i+'\t'+str(exonlen)+"\n")
        fw.close()

@click.command()
@click.option("-i","--input", help="A xx.psi file create by suppa.")
@click.option("-p","--prefix", help="output prefix.")
def filterpsi(input,prefix):
    '''
     filter NA and 0 in all sample line
    '''
    ase=['SE','A3','MX','A5','RI','AF','AL']
    for i in ase:
        out=prefix+"."+i+".psi"
        fw=open(out,'w')
        with open(input,'r') as f:
            for line in f.readlines()[1:]:
                line=line.strip().split()
                geneid=line[0]
                if "NA" in line[1:] or "nan" in line[1:]:
                    continue
                data=map(eval, line[1:])
                if geneid.split(':')[0].endswith(i) and sum(data) >0.0:
                    fw.write("\t".join(line)+"\n")
        fw.close()


main.add_command(cluster)
main.add_command(statas)
main.add_command(filterpsi)

if __name__ == '__main__':
    main()
