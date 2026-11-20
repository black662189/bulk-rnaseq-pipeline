#!/usr/bin/python
# -*- coding:utf-8 -*-

from sys import argv
script, descfile, input, output = argv

dict = {}
description = open(descfile,"r")
for line1 in description:
	line1 = line1.strip("\n")
	data1 = line1.split("\t")
	
	dict[data1[0]] = data1[1]

inf = open(input,"r")
outf = open(output,"w")
for line2 in inf:
	line2 = line2.strip("\n")
	data2 = line2.split("\t")
	gene = data2.pop(0)
	data2 = "\t".join(data2)
	
	if dict.get(gene):
		out = gene + "\t" + dict[gene] + "\t" + data2 + "\n"
		outf.write(out)
	else:
		out = gene + "\tNA\t" + data2 + "\n"
		outf.write(out)
