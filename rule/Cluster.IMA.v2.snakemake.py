#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     Cluster.IMA.v2.snakemake.py
   Description :
   Author :       yueyao
   date：          2019/12/11
-------------------------------------------------
   Change Activity:
                   2019/12/11:
-------------------------------------------------
"""
__author__ = 'yueyao'

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
gnidx=config["database"]["genome"][0]
gtf=config["database"]["gtf"]
dependsdir=config['database']['depends']
spe=config['species']

soapidx=config['database']['gene'][0]
repcompare=config['comparison']
group=config['group']
para2=config['DEseq2']
des=dependsdir+"/gene_description.txt"
readl=config['readlength']

rule all_cluster:
    input:
        expand('{outdir}/4.Comparison/all.genes.normalized.counts.no_dup.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/comparisons.txt',outdir=outdir),
    output:
        expand('{outdir}/5.Cluster/PCA.allsamples.pdf', outdir=outdir),
        expand('{outdir}/5.Cluster/hclust.allsamples.pdf', outdir=outdir),
    params:
        script = dependsdir,
        od = outdir + "/5.Cluster",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time = outdir + "/time/all_cluster.time",
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "group=`cut -f1,3 {input[1]} |sed 's/\\t/\\n/'|sort|uniq|perl -lane 'printf $_.\" \";'`\n" \
        "cat {params.script}/cluster_analysis.R | {params.R} --vanilla --slave --args {input} $group\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"

rule deg_cluster:
    input:
        expand('{outdir}/4.Comparison/all.genes.normalized.counts.no_dup.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/comparisons.txt',outdir=outdir),
    output:
        expand('{outdir}/5.Cluster/PCA.{repcompare}.pdf',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/PCA.{repcompare}.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.{repcompare}.pdf',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.{repcompare}.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/{repcompare}.fpkm.txt',outdir=outdir,repcompare=repcompare),
    threads:1
    params:
        script = dependsdir,
        od = outdir + "/5.Cluster",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time = outdir + "/time/Cluster.time",
        e = outdir + "/log/cluster.e",
        o = outdir + "/log/cluster.o",
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "cat {input[1]}|cut -f1,3|while read id;do cat {params.script}/cluster_analysis.R  |" \
        " {params.R} --vanilla --slave --args {input[0]} $id;done\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"