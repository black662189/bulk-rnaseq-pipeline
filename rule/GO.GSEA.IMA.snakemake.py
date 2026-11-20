#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

'''
@author: dianyu liu
@time: 2023/9/21
@file: GO.GSEA.IMA.snakemake.py
@mail: ldyyyy21@foxmail.com

'''

rule pregsea:
    input:
        '{outdir}/5.Cluster/{repcompare}.fpkm.txt'
    output:
        '{outdir}/7.GO_GSEA/{repcompare}.class.cls'
    params:
        od=outdir+'/7.GO_GSEA'
    run:
        for ss in repcompare:
            with open(params.od+'/'+ss+".class.cls", 'w') as f:
                control=ss.split('_vs_')[0]
                treat=ss.split('_vs_')[1]
                controlnum=len(group.get(control).split(","))
                treatnum=len(group.get(treat).split(","))
                totalsample=controlnum+treatnum
                f.write(str(totalsample)+" 2 1\n")
                f.write("#"+control+" "+treat+"\n")
                f.write(" ".join([control]*controlnum+[treat]*treatnum)+"\n")

dbgmt={
    'mmu':'c2.cp.kegg.v7.5.1.symbols.MMU.gmt',
    'hsa': 'c2.cp.kegg.v7.5.1.symbols.gmt',
    'rno': 'c2.cp.kegg.v7.5.1.symbols.RNO.gmt',
    'mcf':'c2.cp.kegg.v7.5.1.symbols.MFA.gmt',
	'mcc':'mcc.kegg_optimal.symbol.gmt',
    'ssc':'c2.cp.kegg.v7.5.1.symbols.SSC.gmt',
    'ocu':'ocu.kegg_optimal.symbol.gmt',
    'dre':'kegg.gmt',
	'cfa':'c2.cp.kegg.v7.5.1.symbols.CFA.gmt',
	'fca':'c2.cp.kegg.v7.5.1.symbols.FCA.gmt'
}

rule gsea:
    input:
        '{outdir}/5.Cluster/{repcompare}.fpkm.txt',
        '{outdir}/7.GO_GSEA/{repcompare}.class.cls'
    output:
        "{outdir}/7.GO_GSEA/{repcompare}.Gsea.finish",
    params:
        gsea=soft['gsea'][0],
        od=outdir+"/7.GO_GSEA/{repcompare}",
        gmt=dependsdir+"/"+dbgmt[spe]
    log:
        time = outdir + "/time/gsea.time",
        e = outdir + "/log/go_gsea_{repcompare}.e",
        o = outdir + "/log/go_gsea_{repcompare}.o",
    threads:1
    shell:
        #"echo start [ `date +%F` `date +%T` ] >{log.time}\n" \
        "java -Xmx2048m -cp {params.gsea} xtools.gsea.Gsea -res {input[0]} -cls {input[1]} -gmx {params.gmt} " \
        "-out {params.od} -rpt_label {wildcards.repcompare} -collapse false -mode Max_probe -norm meandiv -nperm 1000 " \
        "-permute gene_set-rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Signal2Noise " \
        "-sort real-order descending -include_only_symbols true -make_sets true -median false -num 100 " \
        "-plot_top_x 5000 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 10 -zip_report false -gui false\n" \
        "touch {output}\n" \
        #"echo finish [ `date +%F` `date +%T` ] >>{log.time}"
