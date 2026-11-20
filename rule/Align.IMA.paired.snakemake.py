#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: Align.IMA.paired.snakemake.py
@mail: yueyaomail@gmail.com

'''

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


# rule all:
#     input:
#         expand('{outdir}/2.Alignment/{sample}.sam',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/2.Alignment/{sample}.bam',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/2.Alignment/{sample}.log',outdir=outdir,sample=SAMPLES),

rule alignment:
    input:
        fq1=lambda wildcards: SAMPLES[wildcards.sample][0],
        fq2=lambda wildcards: SAMPLES[wildcards.sample][1],
    output:
        temp('{outdir}/2.Alignment/{sample}.sam'),
        '{outdir}/2.Alignment/{sample}.bam',
        '{outdir}/2.Alignment/{sample}.log'
    params:
        hisat2=soft["hisat2"][0],
        samtools=soft["samtools"][0],
        para=soft["hisat2"][1],
        od=outdir+"/2.Alignment",
        genome=gnidx,
    log:
        e=outdir+"/log/alignment_{sample}.e",
        o=outdir+"/log/alignment_{sample}.o",
        time=outdir+"/time/alignment_{sample}.time"
    threads: 12
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "{params.hisat2} {params.para} -p {threads} -1 {input.fq1} -2 {input.fq2} -x {params.genome} 2>" \
        "{output[2]} -S {output[0]}\n" \
        "{params.samtools} sort -@ {threads} -O BAM -o {output[1]} {output[0]}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
