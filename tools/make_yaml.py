import click
import os
import glob

@click.group()
def main():
    pass

@click.command()
@click.option('-y','--yaml_file_path',default='null',show_default=True,
              help='input raw yaml file path')
@click.option('-s','--species',required=True,
              help='input species')
@click.option('-f','--foldchange',default='2',show_default=True,
              help='input fold change')
@click.option('-e','--end_type',default='SE',show_default=True,
              help='Single-end input SE; paired-end input PE')
@click.option('-r','--rawdata',required=True,
              help='input raw data path')
@click.option('-sa','--sample_info',default='sample.txt',show_default=True,
              help='input sample info file path')
def yaml(yaml_file_path, species, foldchange, end_type, rawdata, sample_info):
    #get yaml file
    if yaml_file_path == 'null':
        yaml_file = '/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/depends_' + species + '/' + species + '.yaml'
    else:
        yaml_file = yaml_file_path

    # rawdata end with '/'
    if rawdata.endswith('/'):
        rawdata = rawdata
    else:
        rawdata = rawdata + '/'
    # get output path
    outpath = os.getcwd()

    out = open(species + '_complemented.yaml', 'w')
    #part1 add raw yaml
    with open(yaml_file) as file:
        for line in file:
            keep = line
            lines = line.strip('\n').split(':')
            if lines[0] ==  'outdir':
                keep = 'outdir: ' + outpath + '\n'
            elif lines[0] == 'PEorSE':
                keep = 'PEorSE: ' + end_type + '\n'
            else:
                keep = keep
            out.writelines(keep)

    #read sample info file
    comparison = []
    group = {}
    DESeq2 = 'DEseq2:\n foldchange: ' + foldchange + '\n padj: 0.05\n'
    samples = {}

    with open(sample_info) as file:
        for line in file:
            lines=line.strip('\n').split('\t')
            if lines[0] == 'comparison':
                comparison.append(lines[1] + '_vs_' + lines[2])
            else:
                if end_type == 'SE':
                    gua = glob.glob(rawdata + lines[0] + '/*.fq*')
                    samples[lines[2]] = gua
                else:
                    gua1 = glob.glob(rawdata + lines[0] + '/*_1.*fq*')
                    gua2 = glob.glob(rawdata + lines[0] + '/*_2.*fq*')
                    samples[lines[2]] = [gua1, gua2]

                if lines[1] in group:
                    group[lines[1]].append(lines[2])
                else:
                    group[lines[1]] = [lines[2]]

    #part2 add comparison
    out.writelines('comparison:\n')
    out.writelines(['  - ' + x + '\n' for x in comparison])
    #part3 add group
    out.writelines('group:\n')
    out.writelines(['  ' + k + ': ' + (',').join(v) + '\n' for k,v in group.items()])
    #part4 add DESeq2
    out.writelines(DESeq2 + '\n')
    #part5 add samples
    out.writelines('samples:\n')
    if end_type == 'SE':
        out.writelines('  ' + k + ':\n' + '    - ' + v[0] + '\n' for k,v in samples.items())
    else:
        out.writelines('  ' + k + ':\n' + '    - ' + v[0][0] + '\n' +  '    - ' + v[1][0] + '\n' for k, v in samples.items())
    out.close()

main.add_command(yaml)

if __name__ == '__main__':
    main()
