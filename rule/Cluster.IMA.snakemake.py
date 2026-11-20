#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: Cluster.IMA.snakemake.py
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


rule cluster:
    input:
        expand('{outdir}/4.Comparison/all.genes.fpkm.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/temp/GroupList.txt',outdir=outdir),
        expand('{outdir}/4.Comparison/temp/CompareList.txt',outdir=outdir),
    output:
        expand('{outdir}/5.Cluster/PCA.{repcompare}.pdf',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/PCA.{repcompare}.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.{repcompare}.pdf',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/hclust.{repcompare}.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/5.Cluster/PCA.allsamples.pdf',outdir=outdir),
        expand('{outdir}/5.Cluster/hclust.allsamples.pdf',outdir=outdir),
        expand('{outdir}/5.Cluster/{repcompare}.fpkm.txt',outdir=outdir,repcompare=repcompare),
    threads:3
    params:
        script = '/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/depends_script',
        od = outdir + "/5.Cluster"
    log:
        time = outdir + "/time/Cluster.time",
        e = outdir + "/log/cluster.e",
        o = outdir + "/log/cluster.o",
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "python {params.script}/RNAseq.IMA.tools.py cluster -i {input[0]} -g {input[1]}" \
        " -c {input[2]} -o {params.od}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
