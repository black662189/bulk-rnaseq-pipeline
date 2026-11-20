#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: SOAP.IMA.single.snakemake.py
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

# rule all:
#     input:
#         expand('{outdir}/1.ReadsRandom/{sample}.soap.list',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/1.ReadsRandom/Reads_distribution/{sample}.ReadsRandom.png',outdir=outdir,sample=SAMPLES),


rule soap:
    input:
        fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
    output:
        temp('{outdir}/1.ReadsRandom/{sample}.soap.gz'),
        temp('{outdir}/1.ReadsRandom/{sample}.soap.list'),
    threads:24
    params:
        soap=soft['soap'][0],
        soap_idx=soapidx,
        soap_para=soft['soap'][1],
        od=outdir+'/1.ReadsRandom/',
    log:
        time = outdir + "/time/soap_{sample}.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "{params.soap} -a {input.fq1} -D {params.soap_idx} {params.soap_para} -p {threads} -o {output[0]} -2 unpaired_{wildcards.sample}.soap.gz\n" \
        "echo \"{output[0]}\" > {output[1]}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"


rule soappic:
    input:
        '{outdir}/1.ReadsRandom/{sample}.soap.gz',
        '{outdir}/1.ReadsRandom/{sample}.soap.list',
    output:
        '{outdir}/1.ReadsRandom/Reads_distribution/{sample}.ReadsRandom.png'
    threads:1
    log:
        time = outdir + "/time/reads_dis_{sample}.time"
    params:
        od = outdir + '/1.ReadsRandom/Reads_distribution',
        script = dependsdir + "/ReadsRandomInGene.pl",
        mm10 = dependsdir + "/mRNA_length.txt"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}/\n" \
        "perl {params.script} {params.mm10} {input[1]} {wildcards.sample}\n" \
        "convert {params.od}/{wildcards.sample}.ReadsRandom.svg {output}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"