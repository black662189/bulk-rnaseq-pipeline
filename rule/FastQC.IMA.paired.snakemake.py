#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: FastQC.IMA.paired.snakemake.py
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


rule fastqc:
    input:
        fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
        fq2=lambda wildcards: SAMPLES[wildcards.sample][1],
    output:
        '{outdir}/0.FastQC/{sample}/{sample}_1.finish',
        '{outdir}/0.FastQC/{sample}/{sample}_2.finish',
    params:
        od=outdir+"/0.FastQC/{sample}",
        fastqc=soft["fastqc"][0],
    log:
        time = outdir + "/time/fastqc_{sample}.time"
    threads:12
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "mkdir -p {params.od}\n" \
        "{params.fastqc} -t {threads} -o {params.od} {input.fq1} {input.fq2}\n" \
        "touch {output[0]}\n" \
        "touch {output[1]}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
