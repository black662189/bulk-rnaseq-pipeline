#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

my $usage=<<INFO;
Description:
	Draw SVG picture for randomness of reads distrubution in genes.
	Need genes' length infomation and genes' soap result.
	Support 'gzip' formatted soap result.
Date:
	2010-10-26
	2010-11-05
	2011-07-13
Author:
	He Zengquan, hezengquan\@genomics.org.cn
	nixiaoming,  nixiaoming\@genomics.org.cn
Modifier:
	Huang Kai, huangkai\@genomics.org.cn
Usage: 
	perl $0 [options] <genelength> <soaplist> <prefix> 
	
	<genelength> *  gene length information, formatted as 'geneID \\t length'
	<soaplist> *    a list of soap result files, each line has one file with its path
	<prefix> *      result prefix, the result file will be named as 'prefix.ReadsRandom.svg'
	
	options:

	-n <int>      how much columns in the figure [100]
	-h|help       output help information to screen
	-width <int>  svg width [800]
	-height <int> svg height [600]
	-color <str>  color of the line[black]
	e.g.:
	perl $0 genelen.list XX.genesoap.list xxx.ReadsRandom.svg -n 20 

INFO

my ($lengthFile,$soaplist,$result,$cNum,$help,$W,$H,$lineColor);
my $debug;
GetOptions(
	"n:i"=>\$cNum,
	"h|help"=>\$help,
	"width=i"=>\$W,
	"height=i"=>\$H,
	"debug" => \$debug,
	"color" => \$lineColor,
);

$lengthFile = shift;
$soaplist = shift;
$result = shift;
$W ||= 800;
$H ||= 600;
$lineColor ||= "black";

die "$usage" if($help||!$lengthFile||!$soaplist||!$result);

$result .=".ReadsRandom.svg";

$cNum ||= 100;   ##columns
my %gene;

my @colValue;      ## $colValue[i] is a counter of the value between i*1/cNum and (i+1)*1/cNum
for (my $i=0; $i<=$cNum ; $i++)
{
	$colValue[$i]=0;
}

## record geneID and the length information, $gene{ID}=length
open(my $GENE,$lengthFile) || die "can't open $lengthFile";
showLog("reading $lengthFile");
while (<$GENE>)
{
	next if(/^\s*$/ || /^\#/);
	chomp;
	my @temp=split(/\t/,$_);
	$gene{$temp[0]}=$temp[1];
}
close($GENE);


## read soap file

##*********************************************************************************
## We get each area's total reads number, draw columns for each area in order.
## If the randomness is good, the reads number will be contiguous, even horizontal.
##*********************************************************************************

open(my $SOAPLIST,$soaplist) || die "can't open $soaplist";
while (<$SOAPLIST>) {
	next if(/^\s*$/ || /^\#/);
	chomp;

	my $soapfile=$1  if(/^(\S+)/);
	my $SOAP;
	my $fileformat=`file $soapfile `;
	chomp $fileformat;
	if ($fileformat =~ /gzip compressed data/) {
		open $SOAP, " gzip -dc $soapfile |" or die "can't open the soap file  $soapfile";
	}elsif ($fileformat =~ /bzip2 compressed data/){
		open $SOAP, " bzip2 -cd $soapfile |" or die "can't open the soap file  $soapfile";
	}else{
		open $SOAP, $soapfile  or die "can't open the soap file  $soapfile";
	}

	showLog("reading $soapfile");
	while (<$SOAP>)
	{
		next if(/^\s*$/ || /^\#/);
		chomp;
		my @temp = split(/\t/,$_);
		if(!exists $gene{$temp[7]})
		{
			print STDERR "can't find the length of gene: $temp[7]\n";
			next;
		}
		my $ratio = $temp[8]/$gene{$temp[7]};
		if($ratio>1)
		{
			print STDERR "$temp[8] is larger than the length of $temp[7]\n";
			next;
		}
		CountReadsNumber(\@colValue,$temp[8],$temp[5],$gene{$temp[7]},$temp[3]);

	}
	close($SOAP);
}
close($SOAPLIST);

# draw SVG picture
$colValue[$cNum-1]+=pop @colValue;### 将 比对到 序列 最后一个 点的加到 最后一个区间

for (my $i=0; $i<$cNum ; $i++)
{
	$colValue[$i] = int($colValue[$i]+0.5);
	print "$i\t$colValue[$i]\n" if (defined $debug);
}

drawSVG($result,\@colValue,$W,$H);

showLog("Finish! $result");

############
sub showLog
{
	my ($info) = @_;
	my @times = localtime; # sec, min, hour, day, month, year
	printf ("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}
############
# count reads number for each columns
sub CountReadsNumber
{
	my ($value,$start,$map_len,$gene_len,$hit)=@_;
	my $end=$start+$map_len-1;
	if ($end > $gene_len)
	{
		$end=$gene_len;
		$map_len=$gene_len-$start+1;
	}
	my $add=1/($map_len*$hit);
	for (my $i=$start;$i <= $end;$i++)
	{
		my $index=int($cNum*$i/$gene_len);
		$$value[$index]+=$add;
	}
}

sub drawSVG
{
	my ($file,$value,$width,$height) = @_;
	my @values = @$value;

	my $rulerLen = 5;

	my $maxValue = Max(@values);
	
	# max power, highest decimal place, max value in figure
	my $maxPower = int( log($maxValue)/log(10) );
	my $maxTimes = int( $maxValue/(10**$maxPower) )+1;
	my $maxH = $maxTimes*(10**$maxPower);

	# margin edge
	my ($edgeTop,$edgeBottom,$edgeLeft,$edgeRight) = ($height/10,$height*9/10,$width/10,$width*9/10);
	# main area
	my $mainLeft = $width*2.5/10;
	my $mainBottom = $height*7/10;

	my $columnWidth = ($edgeRight-$mainLeft)/$cNum;

	my $XaxisTextx = ($mainLeft+$edgeRight)/2;
	my $XaxisTexty = $height*9/10;
	my $YaxisTestx = $edgeLeft;
	my $YaxisTesty = ($edgeTop+$mainBottom)/2;

	my $windowsy = $height*245/300;
	my $windowstext = $cNum."windows";

	my $svg = <<SVGCODE;  # code
<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg width="$width" height="$height" version="1.1" xmlns="http://www.w3.org/2000/svg">
<style type="text/css">
text {font-family: Arial; fill:black;}
line {stroke-width:1; stroke:black;}
</style>

<polyline fill="none" stroke="black" stroke-width="2" 
	points = "$mainLeft,$edgeTop $mainLeft,$mainBottom $edgeRight,$mainBottom" />

<g font-weight="bold" text-anchor="middle">
	<text font-size="10" x="$XaxisTextx" y="$windowsy">$windowstext</text>
	<text font-size="16" x="$XaxisTextx" y="$XaxisTexty">Relative Position in Gene(5\'->3\')</text>
	<text font-size="16" x="$YaxisTestx" y="$YaxisTesty" transform="rotate(-90,$YaxisTestx,$YaxisTesty)">Number of Reads</text>
</g>

SVGCODE

	# draw main area, columns
	$svg.= "<polygon fill=\"none\" stroke=\"$lineColor\" stroke-width=\"1\" \n\tpoints = \"$mainLeft,$mainBottom \n";
	foreach my $i (0..($cNum-1))
	{
		my $thisH = (1-$values[$i]/$maxH)*($mainBottom-$edgeTop)+$edgeTop;
		my $x1 = $mainLeft+$i*$columnWidth;
		my $x2 = $mainLeft+($i+1)*$columnWidth;
		$svg.="\t\t$x1,$thisH $x2,$thisH\n";
	}
	$svg.="\t\t$edgeRight,$mainBottom\"/>\n\n";


	# draw x y axis, rulers
	$svg.="<g stroke=\"black\" stroke-width=\"0.7\">\n";
	my $marks="<g font-size=\"10\" >\n";
	foreach my $i (0,2,4,6,8,10)
	{
		my $x = $mainLeft+($edgeRight-$mainLeft)*$i/10;
		my $value = sprintf("%.1f",$i/10);
		my $y2 = $mainBottom+$rulerLen;
		my $texty = $y2+2*$rulerLen;
		$svg.="\t<line x1=\"$x\" y1=\"$mainBottom\" x2=\"$x\" y2=\"$y2\" />\n";
		$marks.="\t<text x=\"$x\" y=\"$texty\" text-anchor=\"middle\">$value</text>\n";
	}
	foreach my $i (0..$maxTimes)
	{
		my $y = (1-$i/$maxTimes)*($mainBottom-$edgeTop)+$edgeTop;
		my $value = $i*(10**$maxPower);
		$value = ($maxPower<6) ? $value:sprintf("%.1e",$value);
		my $x1 = $mainLeft-$rulerLen;
		my $textx = $x1-2*$rulerLen;
		$svg.="\t<line x1=\"$x1\" y1=\"$y\" x2=\"$mainLeft\" y2=\"$y\" />\n";
		$marks.="\t<text x=\"$textx\" y=\"$y\" text-anchor=\"end\">$value</text>\n";
	}
	$svg.="</g>\n\n";
	$marks.="</g>\n";
	$svg.=$marks;
	$svg.="</svg>\n";



	# print SVG code to result file
	open (my $SVG, ">$file") || die $!;
	print $SVG $svg;
	close($SVG);
}

sub Max{
	my $max=shift;
	foreach  (@_) {
		$max=$_ if($_>$max);
	}
	return $max;
}
