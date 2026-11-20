#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: GO.IMA.snakemake.py
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
    'cfa':'org.Cf.eg.db',
    'fca':'org.Cf.eg.db'
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

rule pregoenrich:
    input:
        expand('{dependsdir}/GO_enrichment_analysis.R',dependsdir=dependsdir),
        expand('{dependsdir}/GO_heatmap.R',dependsdir=dependsdir),
    output:
        "{outdir}/8.GO/GO_enrichment_analysis.R",
        "{outdir}/8.GO/GO_heatmap.R",
    params:
        db=dbs[spe],
        org=organism[spe]
    shell:
        "cp {input[0]} {output[0]}\n" \
        "cp {input[1]} {output[1]}\n" \
        "perl -i -p -e 's/library\(Orgdb\)/library\({params.db}\)/;s/OrgDb=Orgdb/OrgDb={params.db}/' {output[0]}\n" \

rule go:
    input:
        '{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',
        '{outdir}/5.Cluster/{repcompare}.fpkm.txt',
        "{outdir}/8.GO/GO_enrichment_analysis.R",
        "{outdir}/8.GO/GO_heatmap.R",
    output:
        '{outdir}/8.GO/{repcompare}.molecular_function_enrichment.xls',
        '{outdir}/8.GO/{repcompare}.biological_process_enrichment.xls',
        '{outdir}/8.GO/{repcompare}.cellular_component_enrichment.xls',
    threads:1
    params:
        od=outdir+"/8.GO",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time=outdir+"/time/{repcompare}_go.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "control=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[0];print $a;'`\n" \
        "treat=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[1];print $a;'`\n" \
        "cat GO_enrichment_analysis.R | {params.R} --vanilla --slave --args {input[0]} {wildcards.repcompare}\n" \
        "mkdir -p {wildcards.repcompare}_heatmap && cd {wildcards.repcompare}_heatmap\n" \
        "mkdir -p molecular_function biological_process cellular_component\n" \
        "cd molecular_function\n" \
        "cat ../../GO_heatmap.R | {params.R} --vanilla --slave --args {input[1]} {output[0]} $control $treat\n" \
        "cd ../biological_process\n" \
        "cat ../../GO_heatmap.R | {params.R} --vanilla --slave --args {input[1]} {output[1]} $control $treat\n" \
        "cd ../cellular_component\n" \
        "cat ../../GO_heatmap.R | {params.R} --vanilla --slave --args {input[1]} {output[2]} $control $treat\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
