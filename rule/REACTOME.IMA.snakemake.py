#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: REACTOME.IMA.snakemake.py
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
    'ocu':'org.Mm.eg.db',
    'dre':'org.Dr.eg.db',
	'cfa':'org.Hs.eg.db',
    'fca':'org.Hs.eg.db'
}
organism={
    'mmu':'mouse',
    'rno':'rat',
    'hsa':'human',
    'ssc':'human',
    'mcf':'human',
	'mcc':'human',
    'ocu':'mouse',
    'dre':'zebrafish',
	'cfa':'human',
    'fca':'human'
}

rule prereactomeenrich:
    input:
        expand('{dependsdir}/reactome_enrichment_analysis.R',dependsdir=dependsdir),
        expand('{dependsdir}/reactome_heatmap.R',dependsdir=dependsdir),
    output:
        '{outdir}/9.Reactome/reactome_enrichment_analysis.R',
        '{outdir}/9.Reactome/reactome_heatmap.R'
    params:
        db=dbs[spe],
        org=organism[spe]
    shell:
        "cp {input[0]} {output[0]}\n" \
        "cp {input[1]} {output[1]}\n" \
        "perl -i -p -e 's/library\(Orgdb\)/library\({params.db}\)/;s/AnnotationDbi::" \
        "select\(Orgdb/AnnotationDbi::select\({params.db}/;s/organism=organism/organism=\"{params.org}\"/' {output[0]}\n" \


rule reactome:
    input:
        '{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',
        '{outdir}/5.Cluster/{repcompare}.fpkm.txt',
        '{outdir}/9.Reactome/reactome_enrichment_analysis.R',
        '{outdir}/9.Reactome/reactome_heatmap.R'
    output:
        '{outdir}/9.Reactome/{repcompare}.reactome_enrichment.xls'
    params:
        od=outdir+'/9.Reactome',
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    threads:1
    log:
        time=outdir+"/time/{repcompare}_reactome.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "cat reactome_enrichment_analysis.R | {params.R}  --vanilla --slave --args {input[0]} {wildcards.repcompare}\n" \
        "mkdir -p {wildcards.repcompare}_heatmap && cd {wildcards.repcompare}_heatmap\n" \
        "control=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[0];print $a;'`\n" \
        "treat=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[1];print $a;'`\n" \
        "cat ../reactome_heatmap.R | {params.R}  --vanilla --slave --args {input[1]} {output} $control $treat\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
