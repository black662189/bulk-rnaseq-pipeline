#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: yueyao
@time: 2019/7/15
@file: Quantify.IMA.snakemake.py
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



rule stringtie:
    input:
        '{outdir}/2.Alignment/{sample}.bam'
    output:
        '{outdir}/3.Quantify/DESeq2/{sample}/{sample}_abundance.txt',
        '{outdir}/3.Quantify/DESeq2/{sample}/{sample}.gtf'
    threads:12
    params:
        R=soft["R"][0],
        Rscript=soft["Rscript"][0],
        stringtie=soft["stringtie"][0]
    log:
        e = outdir + "/log/stringtie_{sample}.e",
        o = outdir + "/log/stringtie_{sample}.o",
        time=outdir+"/time/stringtie_{sample}.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "{params.stringtie} -p {threads} -G {gtf} -B -e -o {output[1]} -A {output[0]} -l {wildcards.sample} {input} 1>{log.o} 2>{log.e}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"

rule mergefpkm:
    input:
        expand('{outdir}/3.Quantify/DESeq2/{sample}/{sample}.gtf',outdir=outdir,sample=SAMPLES)
    output:
        '{outdir}/3.Quantify/ensembl2symbol.txt',
        '{outdir}/3.Quantify/all.ensembl.fpkm.xls',
        '{outdir}/3.Quantify/gene_symbol_count_matrix.csv',
        '{outdir}/3.Quantify/all.genes.fpkm.xls',
        '{outdir}/3.Quantify/gene_count_matrix.csv',
        '{outdir}/3.Quantify/transcript_count_matrix.csv',
        '{outdir}/3.Quantify/all.ensembl.tpm.xls',
        '{outdir}/3.Quantify/all.genes.tpm.xls',
    threads:1
    log:
        time=outdir+"/time/mergefpkm.time"
    params:
        l=readl,
        script="/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/depends_script",
        od=outdir+"/3.Quantify/DESeq2",
        od2=outdir+"/3.Quantify",
        sd="DESeq2",
        dependpath=dependsdir,
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od2}\n" \
        "/media/nfs/nfs02/yuzhij/biosoft/prepDE.py -l {params.l} -i DESeq2/\n" \
        "grep -P \"\\t(pseudogene|gene)\\t\" {gtf} | perl -p -e \'if(/gene_id \"([^\"]+)\";.+gene_name \"([^\"]+)\"/){{$_=\"$1\\t$2\\n\"}} elsif(/ID=(gene-[^;]+);.+\\bgene=([^;]+)/){{$_=\"$1\\t$2\\n\"}} else{{$_=\"\"}}\'|awk '{{if(NF==2) print $1\"\\t\"$2}}' > {output[0]}\n" \
        "perl {params.script}/merge_stringtie_gene_FPKM_single_table.pl DESeq2/*/*_abundance.txt" \
        " > {output[1]}\n" \
        "perl -i -p -e 's%{params.sd}/%%g' {output[1]}\n" \
        "perl {params.script}/merge_stringtie_gene_TPM_single_table.pl DESeq2/*/*_abundance.txt > {output[6]}\n" \
        "perl -i -p -e 's%{params.sd}/%%g' {output[6]}\n" \
        "cat {params.dependpath}/ensembl2symbol.R | {params.R} --vanilla --slave\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"

rule all_deseq_normal:
    input:
        '{outdir}/3.Quantify/gene_symbol_count_matrix.csv'
    output:
        '{outdir}/3.Quantify/all.genes.normalized.counts.xls'
    params:
        script="/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/depends_script",
        od = outdir+"/3.Quantify",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time=outdir+"/time/mergefpkm.time"
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "cat {params.script}/DESeq2.all.normalized.R | {params.R} --vanilla --slave --args {input}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
