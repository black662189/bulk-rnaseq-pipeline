#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/10
@file: PPI.IMA.snakemake.py
@mail: yueyaomail@gmail.com

'''

SAMPLES = config["samples"]
soft=config["software"]
python=soft['python']
outdir=config["outdir"]
gnidx=config["database"]["genome"][0]
gtf=config["database"]["gtf"]
dependsdir=config['database']['depends']
ppipath=config['database']['ppi']
spe=config['species']

soapidx=config['database']['gene'][0]
repcompare=config['comparison']
group=config['group']
para2=config['DEseq2']
des=dependsdir+"/gene_description.txt"
readl=config['readlength']



def getgroup(repcompare,outdir):
    with open(outdir+'/group.txt','w') as f:
        for line in repcompare:
            f.write(line+"\n")
getgroup(repcompare,outdir)


# rule all:
#     input:
#         expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation.xls',outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation_top_score_100.xls',outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation.xls.100.fd',outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.pdf',outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.png',outdir=outdir, repcompare=repcompare),

rule ppi:
    input:
        '{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',
    output:
        '{outdir}/10.PPI/{repcompare}.differentially_expressed_genes',
        '{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation.xls',
        '{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation_top_score_100.xls',
        '{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.relation.xls.100.fd',
        '{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.pdf',
        '{outdir}/10.PPI/{repcompare}.differentially_expressed_genes.network.png'
    threads:1
    log:
        time=outdir+'/time/{repcompare}_ppi.time',
    params:
        od=outdir+'/10.PPI',
        python3=python[1],
        python2=python[0],
        mmppi=ppipath,
        ppiscript=dependsdir+'/PPI',
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "ln -sf {input} {output[0]}\n" \
        "{params.python3} {params.ppiscript}/get_deg_ppi.py {output[0]} {params.mmppi}\n" \
        "set +e pipefail;awk 'NR>1 {{print $0}}' {output[1]} |sort -nrk 5|head -100 >{output[2]}\n" \
        "{params.python2} {params.ppiscript}/fd_relaion.py {output[2]} {output[0]} >{output[3]}\n" \
        "{params.Rscript} {params.ppiscript}/PPI.R {output[2]} {output[3]} {output[4]}\n" \
#        "convert -density 300 -resize 30% {output[4]} {output[5]}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"

