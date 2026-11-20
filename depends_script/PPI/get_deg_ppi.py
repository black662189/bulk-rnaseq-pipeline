#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import logging
import argparse

""" 
@author:yueyao 
@file: get_deg_ppi.py 
@time: 2019/03/26 
"""

def main(diff,ppidata):
    deglist=[]
    degfilename=os.path.basename(diff)
    degppi=open(degfilename+".network.relation.xls",'w')
    with open(diff,'r') as f:
        for line in f.readlines()[1:]:
            id = line.strip().split("\t")[0]
            deglist.append(id)
    degppi.write("gene1\tgene2\tprotein_cluster1\tprotein_cluster2\tscore\n")
    with open(ppidata,'r') as f:
        for line in f.readlines():
            line=line.strip()
            arr=line.split("\t")
            if arr[0] in deglist and arr[1] in deglist:
                degppi.write(arr[0]+"\t"+arr[1]+"\t"+line+"\n")

if __name__ == '__main__':

    if len(sys.argv) == 1:
        print("python "+ sys.argv[0]+ " <deg.lst> <rename.uniq.10090.protein.links.v10.5.txt>")
        sys.exit()
    deg=sys.argv[1]
    ppidata=sys.argv[2]
    main(deg,ppidata)


