#!/bin/bash
#准备文件
cp -r /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/image/ .

# 生成新图片
/media/nfs/nfs02/yuzhij/.conda/envs/RNA-seq/bin/python /media/nfs/nfs02/yuzhij/Project/HTP/png_make.py

# 替换图片到image目录
if [ -f "analysis_output/sample_info.png" ] && [ -f "analysis_output/analysis_strategy.png" ]; then
    cp analysis_output/sample_info.png image/
    cp analysis_output/analysis_strategy.png image/
    echo "图片替换成功"
else
    echo "警告: 图片生成失败，使用原图片"
fi
rm -rf image/analysis_output
#生成image文件
cd image
for i in `cat ../group.txt` 
do
rm -rf $i 
mkdir $i
cp ../Report/2.Cluster/PCA/PCA\.$i\.png $i/pca\.png
cp ../Report/2.Cluster/PCA/PCA\.$i\.pdf $i/pca\.pdf
cp ../Report/2.Cluster/Hierarchical_clustering/hclust\.$i\.png $i/cluster\.png
cp ../Report/2.Cluster/Hierarchical_clustering/hclust\.$i\.svg $i/cluster\.svg
cp ../Report/2.Cluster/Hierarchical_clustering/hclust\.$i\.pdf $i/cluster\.pdf
cp ../Report/6.Function/KEGG/$i\.kegg_enrichment_category\.png $i/kegg.png
cp ../Report/6.Function/KEGG/$i\.kegg_enrichment_category\.svg $i/kegg.svg
cp ../Report/6.Function/KEGG/$i\.kegg_enrichment_category\.pdf $i/kegg.pdf
cp ../Report/6.Function/Reactome/$i\.reactome_enrichment.svg $i/reactome.svg
cp ../Report/6.Function/Reactome/$i\.reactome_enrichment.png $i/reactome.png
cp ../Report/6.Function/Reactome/$i\.reactome_enrichment.pdf $i/reactome.pdf
cp ../Report/6.Function/GO/$i\.biological_process_enrichment.png $i/biological_process.png
cp ../Report/6.Function/GO/$i\.biological_process_enrichment.svg $i/biological_process.svg
cp ../Report/6.Function/GO/$i\.biological_process_enrichment.pdf $i/biological_process.pdf
cp ../Report/6.Function/GO/$i\.cellular_component_enrichment.svg $i/cellular_component.svg
cp ../Report/6.Function/GO/$i\.cellular_component_enrichment.png $i/cellular_component.png
cp ../Report/6.Function/GO/$i\.cellular_component_enrichment.pdf $i/cellular_component.pdf
cp ../Report/6.Function/GO/$i\.molecular_function_enrichment.svg $i/molecular_function.svg
cp ../Report/6.Function/GO/$i\.molecular_function_enrichment.png $i/molecular_function.png
cp ../Report/6.Function/GO/$i\.molecular_function_enrichment.pdf $i/molecular_function.pdf
cp ../Report/4.Statistics/$i\.DEG_stat\.png $i/kegg_histogram.png
cp ../Report/4.Statistics/$i\.DEG_stat\.pdf $i/kegg_histogram.pdf
### add by yueyao at 2019/12/5
cp ../Report/6.Function/PPI/$i\.differentially_expressed_genes.network.pdf $i/ppi.differentially_expressed_genes.network.pdf
cp ../Report/6.Function/PPI/$i\.differentially_expressed_genes.network.png $i/ppi.differentially_expressed_genes.network.png
cp ../Report/6.Function/PPI/$i\.differentially_expressed_genes.network.svg $i/ppi.differentially_expressed_genes.network.svg
cp ../Report/6.Function/GSVA/DEs_$i\.limma_output.filter.pdf $i/gsva.limma_output.filter.pdf
cp ../Report/6.Function/GSVA/DEs_$i\.limma_output.filter.png $i/gsva.limma_output.filter.png
### add by yueyao at 2019/12/5

done
cd ../

#生成html脚本
rm group_tmp && rm RNA_seq_report.html
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/part1 >RNA_seq_report.html
var2=3
num1=3.3
num2=1
class_seed=0
id_seed=0
for i in `cat group.txt`
do
#生成目录可变部分
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/part2 | 
sed "s/3.3/$num1/g" | 
sed "s/GROUP/$i/g" | 
sed "s/#group1/#group$num2/g" | 
sed "s/#pca1/#pca$num2/g" | 
sed "s/#exp1/#exp$num2/g" | 
sed "s/#gsea1/#gsea$num2/g" | 
sed "s/#kegg1/#kegg$num2/g" |
sed "s/#gsva1/#gsva$num2/g" |
sed "s/#go1/#go$num2/g" |
sed "s/#ppi1/#ppi$num2/g" |
sed "s/#reactome1/#reactome$num2/g"  >>RNA_seq_report.html
#生成结果部分
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part1 | 
sed "s/image\/group1/image\/$i/g" | 
sed "s/3.3/$num1/g" | 
sed "s/GROUP/$i/g" | 
sed "s/group1/group$num2/g" | 
sed "s/pca1/pca$num2/g" | 
sed "s/exp1/exp$num2/g" >>group_tmp
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/pipline/DEG_table.R | /media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R --vanilla --slave --args $i
sed '1d' up_DEG_table.txt | 
sed 's/^/\@/g' | 
sed 's/$/<\/td>\n<\/tr>/g' | 
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp

cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part2 | sed "s/3.3/$num1/g" >>group_tmp
sed '1d' down_DEG_table.txt | 
sed 's/^/\@/g' | 
sed 's/$/<\/td>\n<\/tr>/g' | 
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp

#kegg part
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part3 | 
sed "s/image\/group1/image\/$i/g" | 
sed "s/3.3/$num1/g" | 
sed "s/group1/group$num2/g" | 
sed "s/kegg1/kegg$num2/g" >>group_tmp
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/pipline/kegg_table.R | /media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R --vanilla --slave --args $i
sed '1d' kegg_table.txt | 
sed 's/^/\@/g' | 
sed 's/$/<\/td>\n<\/tr>/g' | 
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp
rm kegg_table.txt


#part4 GSVA add by yueyao at 2019/12/5
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part11 |
sed "s/image\/group1/image\/$i/g" |
sed "s/3.3/$num1/g" |
sed "s/group1/group$num2/g" |
sed "s/gsva1/gsva$num2/g" >> group_tmp
# add by yueyao at 2019/12/5


#part5 GSEA
sample1=`echo $i | sed 's/_vs_/\t/g' | cut -f 1`
sample2=`echo $i | sed 's/_vs_/\t/g' | cut -f 2`
#path1=`ls Report/6.Function/GSEA/$i/*/gsea_report_for_$sample1\_*.html`
#path2=`ls Report/6.Function/GSEA/$i/*/gsea_report_for_$sample2\_*.html`
path1=`ls Report/6.Function/GSEA/$i/*/gsea_report_for*html | grep "gsea_report_for_${sample1}_[0-9]*.html"`
path2=`ls Report/6.Function/GSEA/$i/*/gsea_report_for*html | grep "gsea_report_for_${sample2}_[0-9]*.html"`

cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part5 |
sed "s/3.3/$num1/g" |
sed "s/gsea1/gsea$num2/g" |
sed "s@path1@$path1@g" |
sed "s@path2@$path2@g" |
sed "s@Sample1@$sample1@" |
sed "s@Sample2@$sample2@" >>group_tmp

#part6 GO
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part6 |
sed "s/image\/group1/image\/$i/g" |
sed "s/3.3/$num1/g" |
sed "s/go1/go$num2/g" >>group_tmp

cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/pipline/GO_table.R | /media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R --vanilla --slave --args $i
sed '1d' GO_table_biological_process.txt |
sed 's/^/\@/g' |
sed 's/$/<\/td>\n<\/tr>/g' |
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part7 >>group_tmp
sed '1d' GO_table_cellular_component.txt |
sed 's/^/\@/g' |
sed 's/$/<\/td>\n<\/tr>/g' |
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part8 >>group_tmp
sed '1d' GO_table_molecular_function.txt |
sed 's/^/\@/g' |
sed 's/$/<\/td>\n<\/tr>/g' |
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part9 >>group_tmp

rm GO_table*
#生成button相关class id编号
#ClassType
for ((aa=1;aa<=5;aa++))
do
	sed -i "s/ClassType_classtype${aa}/ClassType_$((aa+class_seed))/g" group_tmp
done
#IdType
for ((bb=1;bb<=12;bb++))
do
	sed -i "s/IdType_idtype${bb}/IdType_$((bb+class_seed))/g" group_tmp
done


#part7 Reactome通路表格
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/reactome_part |
sed "s/image\/group1/image\/$i/g" |
sed "s/3.3/$num1/g" |
sed "s/group1/group$num2/g" |
sed "s/reactome1/reactome$num2/g" >>group_tmp
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/pipline/reactom_table.R | /media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R --vanilla --slave --args $i
sed '1d' reactom_table.txt |
sed 's/^/\@/g' |
sed 's/$/<\/td>\n<\/tr>/g' |
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp
rm reactom_table.txt
#reactome的尾部
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part4 >>group_tmp


#part8 PPI add by yueyao at 2019/12/5
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part10 |
sed "s/image\/group1/image\/$i/g" |
sed "s/3.3/$num1/g" |
sed "s/ppi1/ppi$num2/g" >>group_tmp

cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/pipline/ppi_table.R | /media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R --vanilla --slave --args $i
cat ppi_table.txt|
sed 's/^/\@/g' |
sed 's/$/<\/td>\n<\/tr>/g' |
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>group_tmp
rm ppi_table.txt

cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/group_part12 |
sed "s/image\/group1/image\/$i/g" |
sed "s/3.3/$num1/g" |
sed "s/group1/group$num2/g" |
sed "s/ppi1/ppi$num2/g" >> group_tmp
### add by yueyao at 2019/12/5

var1=3.
var2=`echo $var2 + 1 |bc`
num1=${var1}${var2}
let num2=num2+1
let class_seed=class_seed+5
let id_seed=id_seed+12

done

cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/part3 >> RNA_seq_report.html
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/pipline/fastq_table.R | /media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R --vanilla --slave
sed 's/^/\@/g' QC_table.txt | 
sed 's/$/<\/td>\n<\/tr>/g' | 
awk '{gsub(/\t/,"</td>\n<td>")}{gsub(/\@/,"<tr>\n<td>")}1' >>RNA_seq_report.html
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/part5 >> RNA_seq_report.html
cat group_tmp >> RNA_seq_report.html
cat /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/part7 >> RNA_seq_report.html

rm group_tmp up_DEG_table.txt down_DEG_table.txt kegg_table.txt QC_table.txt GO_table_* reactom_table.txt
tar zcvf html_report.tar.gz image RNA_seq_report.html Report
