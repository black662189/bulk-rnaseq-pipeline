#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     Comparison.IMA.v2.snakemake.py
   Description :
   Author :       yueyao
   date：          2019/12/11
-------------------------------------------------
   Change Activity:
                   2019/12/11:
-------------------------------------------------
"""
__author__ = 'yueyao'

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
#         expand('{outdir}/4.Comparison/comparisons.txt',outdir=outdir),
#         expand('{outdir}/4.Comparison/{repcompare}.total_genes.xls', outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.add_des.xls',outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls', outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/4.Comparison/{repcompare}.DEG_stat.png', outdir=outdir, repcompare=repcompare),
#         expand('{outdir}/4.Comparison/All_comparison.DEG_stat.png', outdir=outdir),
#         expand('{outdir}/4.Comparison/DEG_stat.xls', outdir=outdir),
#         expand('{outdir}/4.Comparison/all.genes.fpkm.xls', outdir=outdir)


import os

rule pre_comparison:
    input:
        expand('{outdir}/3.Quantify/gene_symbol_count_matrix.csv',outdir=outdir),
    output:
        comparison_list = '{outdir}/4.Comparison/comparisons.txt',
    threads:1
    run:
        redir = os.path.dirname(output.comparison_list)
        os.makedirs(redir, exist_ok=True)
        group_dic={}
        for v,k in group.items():
            k_list = k.split(",")
            k_num_list = [n.split("-")[1] for n in k_list]
            group_dic[v] = k_num_list
        with open(redir+"/comparisons.txt",'w') as f:
            for compare_line in repcompare:
                compare_list = compare_line.split("_vs_")
                f.write(compare_list[0]+"\t"+",".join(group_dic.get(compare_list[0]))+"\t"+compare_list[1]+"\t"+",".join(group_dic.get(compare_list[1]))+"\n")

rule deg:
    input:
        exp = expand('{outdir}/3.Quantify/gene_symbol_count_matrix.csv',outdir=outdir),
        comparison = expand('{outdir}/4.Comparison/comparisons.txt',outdir=outdir),
        fpkm=expand('{outdir}/3.Quantify/all.genes.fpkm.xls',outdir=outdir)
    output:
        expand('{outdir}/4.Comparison/{repcompare}.total_genes.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/{repcompare}.DEG_stat.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/4.Comparison/All_comparison.DEG_stat.png',outdir=outdir),
        expand('{outdir}/4.Comparison/DEG_stat.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.fpkm.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.no_dup_exp.fpkm.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.normalized.counts.xls',outdir=outdir),
        expand('{outdir}/4.Comparison/all.genes.normalized.counts.no_dup.xls',outdir=outdir),
    params:
        script = dependsdir,
        foldchange = para2['foldchange'],
        qvalue = para2['padj'],
        desc = des,
        od = outdir + "/4.Comparison",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time = outdir + "/time/deseq2.time",
        e = outdir + "/log/deseq2.e",
        o = outdir + "/log/deseq2.o",
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "cat {params.script}/DESeq2.R| {params.R} --vanilla --slave --args {input.exp} {input.comparison} {params.foldchange} {params.qvalue} 1>{log.o} 2>{log.e}\n" \
        "ln -sf {input.fpkm} .\n" \
        "cat {params.script}/filter_dup_low_exp.R | {params.R} --vanilla --slave --args {input.fpkm} all.genes.no_dup_exp.fpkm.xls\n" \
        "cat {params.script}/filter_dup.R | {params.R} --vanilla --slave --args all.genes.normalized.counts.xls all.genes.normalized.counts.no_dup.xls\n" \
        "python {params.script}/add_description.py {params.desc} all.genes.no_dup_exp.fpkm.xls temp0\n" \
        "mv temp0 all.genes.no_dup_exp.fpkm.xls\n" \
        "python {params.script}/add_description.py {params.desc} all.genes.normalized.counts.no_dup.xls temp0\n" \
        "mv temp0 all.genes.normalized.counts.no_dup.xls\n" \
        "cat {params.script}/DEG_histogram.R | {params.R} --vanilla --slave\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"

rule add_des:
    input:
        '{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.xls'
    output:
        '{outdir}/4.Comparison/{repcompare}.differentially_expressed_genes.add_des.xls'
    params:
        script = '/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/depends_script',
        desc = des,
        od = outdir + "/4.Comparison"
    log:
        time = outdir + "/time/add_des.time",
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "python {params.script}/add_description.py {params.desc} {input} {output}\n" \
        "perl -i -p -e 's/gene_id\\tNA/gene_id\\tdescription/g' {output}\n" \
        "cp {output} {input}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
