#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import click
import yaml

'''
@author: yueyao
@time: 2019/7/15
@file: auto_RNA-Seq_analysis.v4.py
@mail: yueyaomail@gmail.com

'''

@click.group()
def main():
    pass


@click.command()
@click.option("-c","--conf", help="a configfile for snakemake pipeline.")
@click.option("-t","--threads", default=40,help="number of threads you want to use for this snakemake pipeline.")
def makeshell(conf,threads):
    if conf is None or threads is None:
        click.echo("This program is based on python 2.7 and it is used to create a main RNAseq snakemake file...")
        click.echo('Usage:\n\tpython auto.RNAseq.v1.py makeshell -c rnaseq.rna.yaml')
        sys.exit(1)
    config = conf
    with open(config,'r') as f:
        con=yaml.load(f,Loader=yaml.SafeLoader)
    outdir = con['outdir']
    peorse=con['PEorSE']
    steps = con['step']
    t=str(threads)
    if os.path.exists(outdir):
        pass
    else:
        os.mkdir(outdir)
    if os.path.exists(outdir+"/RNAseq.snakemake.py"):
        os.remove(outdir+"/RNAseq.snakemake.py")
    fs=open(outdir+"/RNAseq.snakemake.py",'a')
    # 如果step全为0可能会出现bug
    snakefile = '''#!/usr/bin/env python
# -*- coding:utf-8 -*-


SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
gnidx=config["database"]["genome"][0]
gtf=config["database"]["gtf"]
dependsdir=config['database']['depends']
spe=config['species']
ruledir=config["ruledir"]

peorse=config['PEorSE']
soapidx=config['database']['gene'][0]
repcompare=config['comparison']
group=config['group']
para2=config['DEseq2']
des=dependsdir+"/gene_description.txt"
readl=config['readlength']

FastQC_step=config['step']['FastQC']
ReadsRandom_step=config['step']['ReadsRandom']
Alignment_step=config['step']['Alignment']
Quantify_step=config['step']['Quantify']
Comparison_step=config['step']['Comparison']
Cluster_step=config['step']['Cluster']
GSEA_step=config['step']['GSEA']
GO_step=config['step']['GO']
KEGG_step=config['step']['KEGG']
REACTOME_step=config['step']['REACTOME']
PPI_step=config['step']['PPI']
GSVA_step=config['step']['GSVA']

rule all:
    input:
        '''
    fs.write(snakefile)
    if int(steps['FastQC']) == 1 and peorse == "PE":
        snakefile = '''expand('{outdir}/0.FastQC/{sample}/{sample}_1.finish',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/0.FastQC/{sample}/{sample}_2.finish',outdir=outdir,sample=SAMPLES),
        '''
        fs.write(snakefile)
    if int(steps['FastQC']) == 1 and peorse == "SE":
        snakefile = '''expand('{outdir}/0.FastQC/{sample}/{sample}.finish',outdir=outdir,sample=SAMPLES),
        '''
        fs.write(snakefile)
    if int(steps['ReadsRandom']) == 1:
        snakefile = '''expand('{outdir}/1.ReadsRandom/{sample}.soap.list',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/1.ReadsRandom/Reads_distribution/{sample}.ReadsRandom.png',outdir=outdir,sample=SAMPLES),
        '''
        fs.write(snakefile)
    if int(steps['Alignment']) == 1:
        snakefile = '''expand('{outdir}/2.Alignment/{sample}.bam',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/2.Alignment/{sample}.log',outdir=outdir,sample=SAMPLES),
        '''
        fs.write(snakefile)
    if int(steps['Quantify']) == 1:
        snakefile = '''expand('{outdir}/3.Quantify/DESeq2/{sample}/{sample}_abundance.txt',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/3.Quantify/DESeq2/{sample}/{sample}.gtf',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/3.Quantify/ensembl2symbol.txt',outdir=outdir),
        expand('{outdir}/3.Quantify/all.ensembl.fpkm.xls',outdir=outdir),
        expand('{outdir}/3.Quantify/gene_symbol_count_matrix.csv',outdir=outdir),
        expand('{outdir}/3.Quantify/all.genes.fpkm.xls',outdir=outdir),
        expand('{outdir}/3.Quantify/all.genes.tpm.xls',outdir=outdir),
        expand('{outdir}/3.Quantify/all.ensembl.tpm.xls',outdir=outdir),
        #expand('{outdir}/3.Quantify/all.genes.normalized.counts.xls',outdir=outdir),
        '''
        fs.write(snakefile)
    if int(steps['Comparison']) == 1:
        snakefile = '''expand('{outdir}/4.Comparison/comparisons.txt',outdir=outdir),
        expand('{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.add_des.xls',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.total_genes.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.DEG_stat.png', outdir=outdir, repcompare=repcompare),
        expand('{outdir}/4.Comparison/All_comparison.DEG_stat.png', outdir=outdir),
        expand('{outdir}/4.Comparison/DEG_stat.xls', outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.fpkm.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.no_dup_exp.fpkm.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.normalized.counts.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.normalized.counts.no_dup.xls',outdir=outdir),
        '''
        fs.write(snakefile)
    if int(steps['Cluster']) == 1:
        snakefile = '''expand('{outdir}/5.Cluster/PCA.{repcompare}.pdf', outdir=outdir, repcompare=repcompare),
        expand('{outdir}/5.Cluster/PCA.{repcompare}.png', outdir=outdir, repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.{repcompare}.pdf', outdir=outdir, repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.{repcompare}.png', outdir=outdir, repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.allsamples.pdf', outdir=outdir),
        expand('{outdir}/5.Cluster/PCA.allsamples.pdf', outdir=outdir),
        expand('{outdir}/5.Cluster/{repcompare}.fpkm.txt',outdir=outdir, repcompare=repcompare),
        '''
        fs.write(snakefile)
    if int(steps['GSEA']) == 1:
        snakefile = '''expand('{outdir}/6.GSEA/{repcompare}.class.cls',outdir=outdir, repcompare=repcompare),
        expand("{outdir}/6.GSEA/{repcompare}.Gsea.finish",outdir=outdir, repcompare=repcompare),
        '''
        fs.write(snakefile)
    if int(steps['GO']) == 1:
        snakefile = '''expand("{outdir}/8.GO/GO_enrichment_analysis.R",outdir=outdir, repcompare=repcompare),
        expand("{outdir}/8.GO/GO_heatmap.R",outdir=outdir, repcompare=repcompare),
        expand('{outdir}/8.GO/{repcompare}.molecular_function_enrichment.xls',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/8.GO/{repcompare}.biological_process_enrichment.xls',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/8.GO/{repcompare}.cellular_component_enrichment.xls',outdir=outdir, repcompare=repcompare),
        '''
        fs.write(snakefile)
    if int(steps['KEGG']) == 1:
        snakefile = '''expand('{outdir}/7.KEGG/{repcompare}.differentially_expressed_genes.list',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/7.KEGG/{repcompare}.kegg_enrichment.xls',outdir=outdir, repcompare=repcompare),
        expand("{outdir}/7.KEGG/kegg_heatmap.R",outdir=outdir, repcompare=repcompare),
        expand("{outdir}/7.KEGG/kegg_pathview.R",outdir=outdir, repcompare=repcompare),
        '''
        fs.write(snakefile)
    if int(steps['REACTOME']) == 1:
        snakefile = '''expand('{outdir}/9.Reactome/reactome_enrichment_analysis.R',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/9.Reactome/reactome_heatmap.R',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/9.Reactome/{repcompare}.reactome_enrichment.xls',outdir=outdir, repcompare=repcompare),
        '''
        fs.write(snakefile)
    if int(steps['PPI']) == 1:
        snakefile='''expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation.xls',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation_top_score_100.xls',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation.xls.100.fd',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.pdf',outdir=outdir, repcompare=repcompare),
        expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.png',outdir=outdir, repcompare=repcompare),
        '''
        fs.write(snakefile)

    if int(steps['GSVA']) == 1:
        snakefile='''expand("{outdir}/11.GSVA/all.gsva_score.xls",outdir=outdir),
        expand("{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.pdf",outdir=outdir,repcompare=repcompare),
        expand("{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.png",outdir=outdir,repcompare=repcompare),
        expand("{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.xls",outdir=outdir,repcompare=repcompare),
        expand("{outdir}/11.GSVA/{repcompare}.gsva_score.xls",outdir=outdir,repcompare=repcompare),
        '''
        fs.write(snakefile)

    snakefile = '''
#可以通过flag来控制是否需要运行某一个模块
if FastQC_step == 1 and peorse == "SE":
    include: ruledir+'/FastQC.IMA.single.snakemake.py'
if FastQC_step == 1 and peorse == "PE":
    include: ruledir+'/FastQC.IMA.paired.snakemake.py'
if ReadsRandom_step == 1 and peorse == "PE":
    include: ruledir+'/SOAP.IMA.paired.snakemake.py'
if ReadsRandom_step == 1 and peorse == "SE":
    include: ruledir+'/SOAP.IMA.single.snakemake.py'
if Alignment_step == 1 and peorse == "PE":
    include: ruledir+'/Align.IMA.paired.snakemake.py'
if Alignment_step == 1 and peorse == "SE":
    include: ruledir+'/Align.IMA.single.snakemake.py'
if Quantify_step == 1:
    include: ruledir+'/Quantify.IMA.snakemake.py'
if Comparison_step == 1:
    include: ruledir + '/Comparison.IMA.v2.snakemake.py'
if Cluster_step == 1 :
    include: ruledir + '/Cluster.IMA.v2.snakemake.py'
if GSEA_step == 1:
    include: ruledir + '/GSEA.IMA.snakemake.py'
if KEGG_step == 1:
    include: ruledir + '/KEGG.IMA.snakemake.py'
if REACTOME_step == 1:
    include: ruledir + '/REACTOME.IMA.snakemake.py'
if GO_step == 1:
    include: ruledir + '/GO.IMA.snakemake.py'
if PPI_step == 1:
    include: ruledir + '/PPI.IMA.snakemake.py'
if GSVA_step == 1:
    include: ruledir + '/GSVA.IMA.snakemake.py'
'''
    fs.write(snakefile)
    with open(outdir+"/cp.sh",'w') as cf:
        cpshell='''
#!/bin/bash

rm -rf Report && mkdir Report && cd Report
mkdir -p 1.QC/FastQC 1.QC/Reads_distribution 2.Cluster/PCA 2.Cluster/Hierarchical_clustering 3.DEG 4.Statistics 5.Expression 6.Function/GSEA 6.Function/KEGG 6.Function/GO 6.Function/Reactome 6.Function/PPI 6.Function/GSVA/
cp -R ../0.FastQC/* 1.QC/FastQC/
cp ../1.ReadsRandom/Reads_distribution/*png 1.QC/Reads_distribution/
cp ../4.Comparison/*.total_genes.xls ../4.Comparison/*.differentially_expressed_genes.xls 3.DEG/
cp ../4.Comparison/all.genes.fpkm.xls 5.Expression/
cp ../4.Comparison/all.genes.normalized.counts.xls 5.Expression/
cp ../3.Quantify/all.genes.tpm.xls 5.Expression/
cp ../5.Cluster/PCA.*.pdf ../5.Cluster/PCA.*.png 2.Cluster/PCA/
cp ../5.Cluster/hclust.*.pdf ../5.Cluster/hclust.*.png ../5.Cluster/hclust.*.svg 2.Cluster/Hierarchical_clustering/
grep "overall alignment rate" ../2.Alignment/*.log | perl -p -e 's%../2.Alignment/%%;s%.log:%\t%' > 4.Statistics/Alignment_rate.xls
cp ../4.Comparison/*.DEG_stat.pdf ../4.Comparison/*.DEG_stat.png ../4.Comparison/DEG_stat.xls 4.Statistics/
cp -R ../6.GSEA/* 6.Function/GSEA/
rm -rf 6.Function/GSEA/*.class.cls 6.Function/GSEA/*.txt 6.Function/GSEA/sep[0-9]*
cp -R ../7.KEGG/*.kegg_enrichment.xls ../7.KEGG/*.pdf ../7.KEGG/*.png ../7.KEGG/*.svg ../7.KEGG/*_heatmap/ ../7.KEGG/*_pathview/ 6.Function/KEGG/
rm 6.Function/KEGG/unfiltered.*
cp -R ../8.GO/*.pdf ../8.GO/*.png ../8.GO/*.svg ../8.GO/*.xls ../8.GO/*_heatmap 6.Function/GO/
cp -R ../9.Reactome/*.pdf ../9.Reactome/*.png ../9.Reactome/*.svg ../9.Reactome/*.xls ../9.Reactome/*_heatmap 6.Function/Reactome/
cp -R ../10.PPI/*.pdf ../10.PPI/*.png ../10.PPI/*.svg ../10.PPI/*.network.relation.xls ../10.PPI/*.network.relation_top_score_100.xls 6.Function/PPI/
cp -R ../11.GSVA/*.pdf ../11.GSVA/*.png ../11.GSVA/*.xls 6.Function/GSVA/
cd ..
#tar czvf Report.tar.gz Report
'''
        cf.write(cpshell)
    with open("run.sh",'w') as f:
        shell="/media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq/bin/snakemake -s %(s1)s --configfile %(s2)s -j %(s3)s --latency-wait 200 -k 1>logo.snakemake 2>loge.snakemake\n"
        shell=shell%{'s1':outdir+"/RNAseq.snakemake.py",'s2':config,'s3':t}
        f.write(shell)


main.add_command(makeshell)


if __name__ == '__main__':
    main()
