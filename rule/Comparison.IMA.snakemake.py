#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: Comparison.IMA.snakemake.py
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


# def create_pre(outdir,comparelist,group):
#     group_str = ""
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)
#     with open(outdir + '/GroupList.txt', 'w') as gl:
#         for j in group.keys():
#             gl.write(j + "\t" + group[j] + "\n")
#             group_str += j + ":" + group[j] + " "
#     group_str = group_str.rstrip()
#     repcompare_str = ','.join(comparelist).replace('_vs_', '&')
#     with open(outdir + '/conf.lst', 'w') as conf:
#         conf.write('GeneDiff_Method = DEseq2' + '\n')
#         conf.write('DEseq2_Group = ' + group_str + '\n')
#         conf.write('DEseq2_VS = ' + repcompare_str + '\n')
#         conf.write('DEseq2_Filter = DEseq2' + '\n\n')
#     with open(outdir + '/CompareList.txt', 'w') as cl:
#         for clt in comparelist:
#             arr = clt.split('_vs_')
#             cl.write(group[arr[0]] + "\t" + group[arr[1]] + "\n")
#
# compath=outdir+"/4.Comparison/temp"
# create_pre(compath,repcompare,group)

rule pre:
    input:
        expand('{outdir}/3.Quantify/gene_symbol_count_matrix.csv',outdir=outdir)
    output:
        grouplist='{outdir}/4.Comparison/temp/GroupList.txt',
        comparelist='{outdir}/4.Comparison/temp/CompareList.txt',
        conf='{outdir}/4.Comparison/temp/conf.lst'
    params:
        sample=SAMPLES
    run:
        redir = os.path.dirname(output.grouplist)
        os.makedirs(redir,exist_ok=True)
        group_str=""
        with open(redir+'/GroupList.txt','w') as gl:
            for j in group.keys():
                gl.write(j+"\t"+group[j]+"\n")
                group_str+=j+":"+group[j]+" "
        group_str=group_str.rstrip()
        repcompare_str=','.join(repcompare).replace('_vs_','&')
        with open(redir+'/conf.lst','w') as conf:
            conf.write('GeneDiff_Method = DEseq2' + '\n')
            conf.write('DEseq2_Group = ' + group_str + '\n')
            conf.write('DEseq2_VS = ' + repcompare_str + '\n')
            conf.write('DEseq2_Filter = DEseq2' + '\n\n')
        with open (redir+'/CompareList.txt','w') as cl:
            for clt in repcompare:
                arr = clt.split('_vs_')
                cl.write(group[arr[0]]+"\t"+group[arr[1]]+"\n")


rule deg:
    input:
        exp = expand('{outdir}/3.Quantify/gene_symbol_count_matrix.csv',outdir=outdir),
        grouplist = expand('{outdir}/4.Comparison/temp/GroupList.txt',outdir=outdir),
        comparelist = expand('{outdir}/4.Comparison/temp/CompareList.txt',outdir=outdir),
        fpkm=expand('{outdir}/3.Quantify/all.genes.fpkm.xls',outdir=outdir)
    output:
        expand('{outdir}/4.Comparison/{repcompare}.total_genes.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.DEG_stat.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/All_comparison.DEG_stat.png',outdir=outdir),
        expand('{outdir}/4.Comparison/DEG_stat.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.fpkm.xls',outdir=outdir)
    threads:1
    params:
        script = '/media/sdc/yueyao/pipeline/RNAseq/depends_dir/RNAseq',
        para = para2,
        desc = des,
        od = outdir + "/4.Comparison"
    log:
        time = outdir + "/time/DEseq2.time",
        e = outdir + "/log/deg.e",
        o = outdir + "/log/deg.o",
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n"
        "perl {params.script}/DEseq2.IMA.RNAseq.pl -count {input.exp} -diff {input.comparelist} -group {input.grouplist} " \
        "{params.para} -desc {params.desc} -outdir {params.od} 1>{log.o} 2>{log.e}\n" \
        "cd {params.od}\n" \
        "ln -sf {input.fpkm} .\n" \
        "python {params.script}/add_description.py {params.desc} all.genes.fpkm.xls temp0\n" \
        "mv temp0 all.genes.fpkm.xls\n" \
        "cat {params.script}/DEG_histogram.R | R --vanilla --slave\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"