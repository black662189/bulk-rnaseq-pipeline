#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: KEGG.IMA.snakemake.py
@mail: yueyaomail@gmail.com

'''

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

dbs={
    'mmu':'org.Mm.eg.db',
    'hsa':'org.Hs.eg.db',
    'rno':'org.Rn.eg.db',
    'mcf':'org.Hs.eg.db',
	'mcc':'org.Hs.eg.db',
    'ssc':'org.Ss.eg.db',
    'ocu':'Orgdb',
    'dre':'org.Dr.eg.db',
	'cfa':'org.Cf.eg.db',
	'fca':'org.Cf.eg.db'
}
organism={
    'mmu':'mouse',
    'rno':'rat',
    'hsa':'human',
    'ssc':'pig',
    'mcf':'human',
	'mcc':'human',
    'ocu':'rabbit',
    'dre':'zebrafish',
	'cfa':'Canine',
	'fca':'cat'
}


rule prekeggenrich:
    input:
        expand('{dependsdir}/kegg_heatmap.R',dependsdir=dependsdir),
        expand('{dependsdir}/kegg_pathview.R',dependsdir=dependsdir),
    output:
        "{outdir}/7.KEGG/kegg_heatmap.R",
        "{outdir}/7.KEGG/kegg_pathview.R",
    params:
        db=dbs[spe],
        org=organism[spe],
        spe = spe,
        script = dependsdir,
        od = outdir+"/7.KEGG",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time = outdir + "/time/prekeggenrich.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "#cat {params.script}/kegg_pathway_download.R | {params.R} --vanilla --slave --args {params.spe}\n" \
        "cp {input[0]} {output[0]}\n" \
        "cp {input[1]} {output[1]}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"


rule kegg:
    input:
        '{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',
        '{outdir}/5.Cluster/{repcompare}.fpkm.txt',
        "{outdir}/7.KEGG/kegg_heatmap.R",
        "{outdir}/7.KEGG/kegg_pathview.R",
    output:
        '{outdir}/7.KEGG/{repcompare}.differentially_expressed_genes.list',
        '{outdir}/7.KEGG/{repcompare}.kegg_enrichment.xls'
    threads:1
    params:
        sp=spe,
        rscript=dependsdir+"/kegg_pathway_enrichment_analysis.R",
        od=outdir+"/7.KEGG",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time=outdir+"/time/{repcompare}_kegg.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "cut -f 1 {input[0]} >{output[0]}\n" \
        "{params.Rscript} {params.rscript} {output[0]} {params.sp} {wildcards.repcompare}\n" \
        "mkdir -p {params.od}/{wildcards.repcompare}_heatmap && cd {params.od}/{wildcards.repcompare}_heatmap\n" \
        "control=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[0];print $a;'`\n" \
        "treat=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[1];print $a;'`\n" \
        "cat {params.od}/kegg_heatmap.R | {params.R} --vanilla --slave --args {input[1]} {output[1]} $control $treat\n" \
        "mkdir -p {params.od}/{wildcards.repcompare}_pathview && cd {params.od}/{wildcards.repcompare}_pathview\n" \
        "cat {params.od}/kegg_pathview.R | {params.R} --vanilla --slave --args {input[0]} {output[1]} {params.sp}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
