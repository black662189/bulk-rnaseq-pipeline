#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
-------------------------------------------------
   File Name：     GSVA.IMA.snakemake.py
   Description :
   Author :       yueyao
   date：          2019/11/29
-------------------------------------------------
   Change Activity:
                   2019/11/29:
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

dbgmt={
    'mmu':'mmu.kegg_optimal.symbol.gmt',
    'hsa': 'hsa.kegg_optimal.symbol.gmt',
    'rno': 'rno.kegg_optimal.symbol.gmt',
    'mcf': 'mcf.kegg_optimal.symbol.gmt',
	'mcc': 'mcc.kegg_optimal.symbol.gmt',
    'ssc': 'ssc.kegg_optimal.symbol.gmt',
    'ocu': 'ocu.kegg_optimal.symbol.gmt',
    'dre':'kegg.gmt',
	'cfa':'cfa.kegg_optimal.symbol.gmt',
	'fca':'fca.kegg_optimal.symbol.gmt'
}



# rule all:
#     input:
#         expand("{outdir}/11.GSVA/all.gsva_score.xls",outdir=outdir),
#         expand("{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.pdf",outdir=outdir,repcompare=repcompare),
#         expand("{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.png",outdir=outdir,repcompare=repcompare),
#         expand("{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.xls",outdir=outdir,repcompare=repcompare),
#         expand("{outdir}/11.GSVA/{repcompare}.gsva_score.xls",outdir=outdir,repcompare=repcompare),


rule gsva_all:
    input:
        '{outdir}/4.Comparison/all.genes.normalized.counts.no_dup.xls',
    output:
        "{outdir}/11.GSVA/all.gsva_score.xls",
    params:
        gsva=soft['gsva'][0],
        od=outdir+"/11.GSVA",
        gmt=dependsdir+"/"+dbgmt[spe],
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time = outdir + "/time/gsva_all.time",
        e = outdir + "/log/gsva_all.e",
        o = outdir + "/log/gsva_all.o",
    threads:1
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "{params.Rscript} {params.gsva} {params.gmt} {input}\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"

rule gsva_deg:
    input:
        '{outdir}/4.Comparison/all.genes.normalized.counts.no_dup.xls',
        "{outdir}/11.GSVA/all.gsva_score.xls",
    output:
        "{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.pdf",
        "{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.png",
        "{outdir}/11.GSVA/DEs_{repcompare}.limma_output.filter.xls",
        "{outdir}/11.GSVA/{repcompare}.gsva_score.xls",
    params:
        gsva = soft['gsva'][1],
        gsva_all = soft['gsva'][0],
        od = outdir + "/11.GSVA",
        gmt = dependsdir + "/" + dbgmt[spe],
        all_gsva = outdir + "/11.GSVA/all.gsva_score.xls",
        R=soft["R"][0],
        Rscript=soft["Rscript"][0]
    log:
        time = outdir + "/time/gsva_{repcompare}.time",
        e = outdir + "/log/gsva_{repcompare}.e",
        o = outdir + "/log/gsva_{repcompare}.o",
    threads:1
    shell:
        "echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "cd {params.od}\n" \
        "control=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[0];print $a;'`\n" \
        "treat=`perl -e '$a=(split/_vs_/,\"{wildcards.repcompare}\")[1];print $a;'`\n" \
        "{params.Rscript} {params.gsva} {params.gmt} {input[0]} $control $treat\n" \
        "echo finish [ `date +%F` `date +%T` ] >>{log.time}"
