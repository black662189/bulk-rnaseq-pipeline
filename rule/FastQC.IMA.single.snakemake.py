#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: FastQC.IMA.single.snakemake.py
@mail: yueyaomail@gmail.com

'''

# configfile: "rnaseq.yaml"
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

#
# rule all:
#     input:
#         expand('{outdir}/0.FastQC/{sample}/{sample}_1_fastqc.html',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/0.FastQC/{sample}/{sample}_1_fastqc.zip',outdir=outdir,sample=SAMPLES)

rule fastqc:
    input:
        fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
    output:
        '{outdir}/0.FastQC/{sample}/{sample}.finish',
    params:
        od=outdir+"/0.FastQC/{sample}",
        fastqc=soft["fastqc"][0],
    log:
        time = outdir + "/time/fastqc_{sample}.time"
    threads:12
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "mkdir -p {params.od}\n" \
        "{params.fastqc} -t {threads} -o {params.od} {input.fq1}\n" \
        "touch {output}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
