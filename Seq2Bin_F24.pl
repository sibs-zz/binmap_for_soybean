######### EXAMPLE ########
##perl Seq2Bin.pl Lane1.AAT.PE.fastq.rlt jap_v4_length_list
##########################
#use strict;
#use warnings;
my @hang;
my @array;
my @chrom;
my @bin;
my @len_detail;
my @temp;

my $parent1=0;
my $parent2=0;
my $win_size=0;
my $head=0;
my $head_edge=0;
my $max_n=0;
my $hetero_12="";
my $hetero_key="";
my $hetero_start="";
my $hetero_end="";
my $line="";
my $line1="";
my $chromosome1="";
my $cstart="";
my $cend="";
my $number1="";
my $chrom_len="";
my $lane_n=0;
my $bstart="";
my $bend="";
my $origin="";
my $number="";
my $c=0;
my $h=0;
my $len=0;
my $round=0;
my $name="";
# Optional ARGV[5]: minimum bin size (bp); default 300000
my $min_bin = (defined $ARGV[5] && $ARGV[5] =~ /^\d+$/) ? int($ARGV[5]) : 300000;
#####################################################
#Judge the edges of every bin.#
#####################################################
open INPUT,"<$ARGV[0]" or die "$!";

my $n=0;
my $m=1;
$chrom[0]=0;
my $chromosome="chromosome01";
while (<INPUT>) {
	$line = $_;
	chomp($line);
	@hang = split(/\t/,$line);
	$lane_n=$#hang;
	for (0..($lane_n)){
		$array[$n][$_]=$hang[$_];
	}
	$n++;
}
if ($n<5000){
	$win_size=15;
}elsif($n>=5000 && $n<10000){
	$win_size=25;
}elsif($n>=10000 && $n<20000){
	$win_size=39;
}elsif($n>=20000 && $n<100000){
	$win_size=59;
}elsif($n>=100000){
	$win_size=99;
}
my $half_win=int($win_size/2);
my $win_size_dom=int($win_size*$ARGV[2]);
#my $head_tail_dom=int($ARGV[2]*$half_win);
#print $win_size_dom."\n";
#print $head_tail_dom."\n";
open INPUT2,"<$ARGV[1]" or die "$!";
while (<INPUT2>) {
	$line1=$_; 
	chomp($line1);
	($chromosome1,$chrom_len)=split(/\s+/,$line1);
	$number1=$chromosome1;
	$number1=~s/chromosome//;
	$number1=~s/^0//;
	$len_detail[$number1]=$chrom_len;
}
my $all_chromosome=$number1;
close INPUT2;

for (0..$n) {
	if ($array[$_][2] eq "$chromosome"){
		$chrom[$m]++;
	}else{
		$chromosome=$array[$_][2];
		$chrom[$m+1]++;
		$m++;
	}
}
for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1];
}

#for (1..$all_chromosome){
#	print "\$chrom[$_]=$chrom[$_]\n";
#}

#calculate ratio and dominant#
for (1..$all_chromosome){
	$c=$_;
	for (($chrom[$c-1]+$half_win)..($chrom[$c]-$half_win-1)) {
		my $win_start=$_;
		for ((($win_start)-$half_win)..($win_start+$half_win)){
			if ($array[$_][0] eq "P1"){
				$parent1++;
			}elsif($array[$_][0] eq "P2"){
				$parent2++;
			}
		}
		
		if ($parent1>$parent2){
			$array[$win_start][5]="parent1";
		}elsif ($parent1<$parent2){
			$array[$win_start][5]="parent2";
		}
		if ($parent1>$win_size_dom){
			$array[$win_start][6]="parent1";
			$array[$win_start][7]=$parent1.":".$parent2;
		}elsif ($parent1<($win_size-$win_size_dom)){
			$array[$win_start][6]="parent2";
			$array[$win_start][7]=$parent1.":".$parent2;
		}else {
			$array[$win_start][6]="hetero";
			$array[$win_start][7]=$parent1.":".$parent2;
		}
		$parent1=0;
		$parent2=0;
	}
	
	
	
#head half bin of every chromosome#

	for (0..$half_win-1){
		$range=$_;
		$head_tail_dom=int($ARGV[2]*(2*$range+1));
#		print $head_tail_dom."\n";
		for (($chrom[$c-1])..($chrom[$c-1]+$range*2)) {
			$hhh=$_;
			if ($array[$hhh][0] eq "P1"){
				$parent1++;
			}elsif($array[$hhh][0] eq "P2"){
				$parent2++;
			}
		}

		if ($parent1>$parent2){
			$array[$chrom[$c-1]+$range][5]="parent1";
		}elsif ($parent1<$parent2){
			$array[$chrom[$c-1]+$range][5]="parent2";
		}
		if ($parent1>$head_tail_dom){
#			$array[$chrom[$c-1]+$range][6]="parent1";
			$array[$chrom[$c-1]+$range][7]=$parent1.":".$parent2;
		}elsif ($parent1<(2*$range+1-$head_tail_dom)){
#			$array[$chrom[$c-1]+$range][6]="parent2";
			$array[$chrom[$c-1]+$range][7]=$parent1.":".$parent2;
		}else {
#			$array[$chrom[$c-1]+$range][6]="hetero";
			$array[$chrom[$c-1]+$range][7]=$parent1.":".$parent2;
		}
		$parent1=0;
		$parent2=0;
		$head_tail_dom=0;
	}
#head[6]	
	$head_tail_dom=int($ARGV[2]*$half_win);
	for (($chrom[$c-1])..($chrom[$c-1]+$half_win-1)) {
		$hhh=$_;
		if ($array[$hhh][0] eq "P1"){
			$parent1++;
		}elsif($array[$hhh][0] eq "P2"){
			$parent2++;
		}
	}
		
	if ($parent1>$head_tail_dom){
		for (($chrom[$c-1])..($chrom[$c-1]+$half_win-1)) {
			$array[$_][6]="parent1";
		}
	}elsif ($parent1<($half_win-$head_tail_dom)){
		for (($chrom[$c-1])..($chrom[$c-1]+$half_win-1)) {
			$array[$_][6]="parent2";
		}
	}else {
		for (($chrom[$c-1])..($chrom[$c-1]+$half_win-1)) {
			$array[$_][6]="hetero";
		}
	}
	for (($chrom[$c-1])..($chrom[$c-1]+4)){
		$array[$_][5]=$array[$chrom[$c-1]][6];
	}
	
	$array[$chrom[$c-1]][9]=$array[$chrom[$c-1]][6];
	$parent1=0;
	$parent2=0;

#tail half bin of every chromosome#
#for ($ii=$chrom[$c]-1;$ii>=$chrom[$c-1];$ii--) {
	for ($ii=$half_win-1;$ii>=0;$ii--){
		$range=$ii;
		$head_tail_dom=int($ARGV[2]*(2*$range+1));
#		print $head_tail_dom."\n";
		for (($chrom[$c]-1-$range*2)..($chrom[$c]-1)) {
			$hhh=$_;
			if ($array[$hhh][0] eq "P1"){
				$parent1++;
			}elsif($array[$hhh][0] eq "P2"){
				$parent2++;
			}
		}

		if ($parent1>$parent2){
			$array[$chrom[$c]-1-$range][5]="parent1";
		}elsif ($parent1<$parent2){
			$array[$chrom[$c]-1-$range][5]="parent2";
		}
		if ($parent1>$head_tail_dom){
			$array[$chrom[$c]-1-$range][6]="parent1";
			$array[$chrom[$c]-1-$range][7]=$parent1.":".$parent2;
		}elsif ($parent1<(2*$range+1-$head_tail_dom)){
			$array[$chrom[$c]-1-$range][6]="parent2";
			$array[$chrom[$c]-1-$range][7]=$parent1.":".$parent2;
		}else {
			$array[$chrom[$c]-1-$range][6]="hetero";
			$array[$chrom[$c]-1-$range][7]=$parent1.":".$parent2;
		}
		$parent1=0;
		$parent2=0;
		$head_tail_dom=0;
	}
	
#tail[6]
	$head_tail_dom=int($ARGV[2]*$half_win);
	for (($chrom[$c]-$half_win)..($chrom[$c]-1)) {
		$hhh=$_;
		if ($array[$hhh][0] eq "P1"){
			$parent1++;
		}elsif($array[$hhh][0] eq "P2"){
			$parent2++;
		}
	}
		
	if ($parent1>$parent2){
		for (($chrom[$c]-$half_win)..($chrom[$c]-1)) {
			$array[$_][5]="parent1";
		}
	}elsif ($parent1<$parent2){
		for (($chrom[$c]-$half_win)..($chrom[$c]-1)) {
			$array[$_][5]="parent2";
		}
	}
	if ($parent1>$head_tail_dom){
		for (($chrom[$c]-$half_win)..($chrom[$c]-1)) {
			$array[$_][6]="parent1";
		}
	}elsif ($parent1<($half_win-$head_tail_dom)){
		for (($chrom[$c]-$half_win)..($chrom[$c]-1)) {
			$array[$_][6]="parent2";
		}
	}else {
		for (($chrom[$c]-$half_win)..($chrom[$c]-1)) {
			$array[$_][6]="hetero";
		}
	}
	for (($chrom[$c]-5)..($chrom[$c]-1)){
		$array[$_][5]=$array[$chrom[$c]-1][6];
	}
	$array[$chrom[$c]-1][9]=$array[$chrom[$c]-1][6];
	$parent1=0;
	$parent2=0;
}
	
###### get raw edge ######
for (0..$chrom[$all_chromosome]) {
	my $edge=$_;
	if (($array[$edge][5] ne $array[$edge+1][5]) && ($array[$edge][2] eq $array[$edge+1][2])){
		$array[$edge][9]=$array[$edge][5];
		$array[$edge+1][9]=$array[$edge+1][5];
	}
}
### Adjust the edges for parent1/parent2
my $key1=1;
my $key2=1;
my $max=0;
RRR666: for (0..$chrom[$all_chromosome]) {
	my $edge=$_;
	if (($array[$edge][9] ne "") && ($array[$edge+1][9] ne "") && ($array[$edge][9] ne $array[$edge+1][9]) && ($array[$edge][2] eq $array[$edge+1][2])){
		$adjust_size=int($half_win/2);
		for (1..$adjust_size){
			if (($array[$edge-$_][9] ne "") or ($array[$edge+$_+1][9] ne "")){
				next RRR666;
			}
		}
		$edge_1=$array[$edge][9];
		$edge_2=$array[$edge+1][9];
		$array[$edge][10]=$edge_1;
		$array[$edge+1][10]=$edge_2;
		$array[$edge][9]="";
		$array[$edge+1][9]="";
		
RRR6:		for(-$adjust_size..$adjust_size){
			if ($array[$edge][2] ne $array[$edge+$_][2]){
				next RRR6;
			}
			if ($array[$edge+$_][9] ne ""){
				next RRR6;
			}
			$temp=$_;
#			print $array[$edge][2]."\t".$array[$edge][1]."\t".$array[$edge][3]."\n";
			for(1..$adjust_size-1){
				if ($array[$edge-$_+$temp][0] eq $array[$edge+$temp][0]){
					$key1++;
				}elsif($array[$edge-$_+$temp][0] ne $array[$edge+$temp][0]){
#					last;
				}
			}
			for(1..$adjust_size-1){
				if ($array[$edge+$_+$temp][0] ne $array[$edge+$temp][0]){
					$key2++;
				}elsif($array[$edge+$_+$temp][0] eq $array[$edge+$temp][0]){
#					last;
				}
			}
			$array[$edge+$temp][8]=$key1+$key2;
			if ($array[$edge+$temp][8]>=$max){
				$max=$array[$edge+$temp][8];
				$max_n=$edge+$temp;
			}
			$key1=1;
			$key2=1;
		}
		if (($array[$max_n][9] eq "") && ($array[$max_n+1][9] eq "")){
			$array[$max_n][9]=$edge_1;
			$array[$max_n+1][9]=$edge_2;
			$max=0;
		}else{
			$array[$edge][9]=$array[$edge][10];
			$array[$edge+1][9]=$array[$edge+1][10];
			$max=0;
		}
#		elsif($array[$edge][8]<5){
#			$array[$edge][9]=$array[$edge][5];
#			$array[$edge+1][9]=$array[$edge+1][5];
#			$max=0;
#		}
	}
}


### Judge the heterous region.
RRR3: for (0..$chrom[$all_chromosome]) {
	if (($array[$_][6] eq "hetero") && ($array[$_-1][6] ne "hetero")){
		$hetero_start=$_;
		$hetero_key=0;
		$hetero_12=$array[$hetero_start][5];
		next RRR3;
	}
	if (($array[$_][6] eq "hetero") && ($array[$_-1][6] eq "hetero")){
		if ($hetero_12 eq $array[$_][5]){
			next RRR3;
		}elsif ($hetero_12 ne $array[$_][5]){
			$temp[$hetero_key]=$_;
			$hetero_key++;
			
			$hetero_12=$array[$_][5];
			next RRR3;
		}
	}
	if (($array[$_][6] ne "hetero") && ($array[$_-1][6] eq "hetero")){
		$hetero_end=$_-1;
		if ($hetero_key==2 && ($temp[1]-$temp[0]<int($half_win/2) )){
			$array[$hetero_start-1][9]=$array[$hetero_start-1][6];
			$array[$hetero_end+1][9]=$array[$hetero_end+1][6];
			$array[$hetero_start][9]="hetero";
			$array[$hetero_end][9]="hetero";
			for ($hetero_start+1..$hetero_end-1){
				$array[$_][9]="";
			}
		}
		if ($hetero_key>2){
			$array[$hetero_start-1][9]=$array[$hetero_start-1][6];
			$array[$hetero_end+1][9]=$array[$hetero_end+1][6];
			$array[$hetero_start][9]="hetero";
			$array[$hetero_end][9]="hetero";
			for ($hetero_start+1..$hetero_end-1){
				$array[$_][9]="";
			}
		}
		$hetero_key=0;
		$hetero_start="";
		$hetero_end="";
		next RRR3;
	}
}

#hetero_head
my $hetero_head;
LLL1: for (1..$all_chromosome){
	$c=$_;
	$hetero_head_start="";
	$hetero_head_end="";
RRR31:	for ($chrom[$c-1]..$chrom[$c]-1) {
		if (($array[$_][6] ne "hetero") && ($_==$chrom[$c-1])){
			next LLL1;
		}
		if (($array[$_][6] eq "hetero") && ($_==$chrom[$c-1])){
#			print $array[$_][2]."\t".$array[$_][1]."\n";
			$hetero_head_start=$_;
			$array[$_][9]="hetero";
			next RRR31;
		}
		if (($array[$_][6] ne "hetero") && ($array[$_-1][6] eq "hetero")){
#			print "end\t".$array[$_][2]."\t".$array[$_][1]."\n";
			$hetero_head_end=$_-1;
			$array[$_-1][9]="hetero";
			$array[$_][9]=$array[$_][6];
			for ($hetero_head_start+1..$hetero_head_end-1){
				$array[$_][9]="";
			}
			$hetero_head_start="";
			$hetero_head_end="";
			next LLL1;
		}
		
	}
}
		
#hetero_tail
LLL2: for (1..$all_chromosome){
#LLL2: for (1..1){
	$c=$_;
	$hetero_tail_start="";
	$hetero_tail_end="";
#	print "$chrom[$c]-1"."\t".$chrom[$c-1]."\n";
RRR32:	for ($ii=$chrom[$c]-1;$ii>=$chrom[$c-1];$ii--) {
#		print "\$ii=".$ii."\n";
		if (($array[$ii][6] ne "hetero") && ($ii==$chrom[$c]-1)){
#			print $array[$ii][2]."\t".$array[$ii][1]."\t".$array[$ii][6]."\n";
			next LLL2;
		}
		if (($array[$ii][6] eq "hetero") && ($ii==$chrom[$c]-1)){
#			print $array[$ii][2]."\t".$array[$ii][1]."\n";
			$hetero_tail_end=$ii;
			$array[$ii][9]="hetero";
			next RRR32;
		}
		if (($array[$ii][6] eq "hetero") && ($array[$ii-1][6] ne "hetero")){
#			print "start\t".$array[$ii][2]."\t".$array[$ii][1]."\n";
			$hetero_tail_start=$ii;
#			print $hetero_tail_start."\n";
			$array[$ii][9]="hetero";
			$array[$ii-1][9]=$array[$ii-1][6];
			for ($hetero_tail_start+1..$hetero_tail_end-1){
				$array[$_][9]="";
			}
			$hetero_tail_start="";
			$hetero_tail_end="";
			next LLL2;
		}
		
	}
}
	
### Adjust the edges for hetero/parent1 and hetero/parent2
my $key1=1;
my $key2=1;
my $max=0;
RRR666: for (0..$chrom[$all_chromosome]) {
	my $edge=$_;
	if (($array[$edge][9] ne "") && ($array[$edge+1][9] ne "") && (($array[$edge][9] eq "hetero") or ($array[$edge+1][9] eq "hetero")) && ($array[$edge][9] ne $array[$edge+1][9]) && ($array[$edge][2] eq $array[$edge+1][2])){
		$adjust_size=int($half_win/2);
		for (1..$adjust_size){
			if (($array[$edge-$_][9] ne "") or ($array[$edge+$_+1][9] ne "")){
				next RRR666;
			}
		}
		$edge_1=$array[$edge][9];
		$edge_2=$array[$edge+1][9];
		$array[$edge][10]=$edge_1;
		$array[$edge+1][10]=$edge_2;
		$array[$edge][9]="";
		$array[$edge+1][9]="";
		
RRR6:		for(-$adjust_size..$adjust_size){
			if ($array[$edge][2] ne $array[$edge+$_][2]){
				next RRR6;
			}
			if ($array[$edge+$_][9] ne ""){
				next RRR6;
			}
			$temp=$_;
#			print $array[$edge][2]."\t".$array[$edge][1]."\t".$array[$edge][3]."\n";
			for(1..$adjust_size-1){
				if ($array[$edge-$_+$temp][0] eq $array[$edge+$temp][0]){
					$key1++;
				}elsif($array[$edge-$_+$temp][0] ne $array[$edge+$temp][0]){
					$key1--;
					last;
				}
			}
			for(1..$adjust_size-1){
				if ($array[$edge+$_+$temp][0] ne $array[$edge+$temp][0]){
					$key2++;
				}elsif($array[$edge+$_+$temp][0] eq $array[$edge+$temp][0]){
					$key2--;
					last;
				}
			}
			$array[$edge+$temp][8]=$key1+$key2;
			if ($array[$edge+$temp][8]>=$max){
				$max=$array[$edge+$temp][8];
				$max_n=$edge+$temp;
			}
			$key1=1;
			$key2=1;
		}
		if (($array[$max_n][9] eq "") && ($array[$max_n+1][9] eq "")){
			$array[$max_n][9]=$edge_1;
			$array[$max_n+1][9]=$edge_2;
			$max=0;
		}else{
			$array[$edge][9]=$array[$edge][10];
			$array[$edge+1][9]=$array[$edge+1][10];
			$max=0;
		}
#		elsif($array[$edge][8]<5){
#			$array[$edge][9]=$array[$edge][5];
#			$array[$edge+1][9]=$array[$edge+1][5];
#			$max=0;
#		}
	}
}


#for (1..$all_chromosome){
#	$c=$_;
#	$array[$chrom[$c-1]][9]=$array[$chrom[$c-1]][6];
#	$array[$chrom[$c]-1][9]=$array[$chrom[$c]-1][6];
#}
for (0..$chrom[$all_chromosome]) {
	if (($array[$_][9] ne "") && ($array[$_-1][9] ne "") && ($array[$_+1][9] ne "")){
		print "Warning\t".$array[$_][2]."\t".$array[$_][1]."\t".$array[$_][3]."\t".$array[$_][9]."\n";
	}
}
open OUT1, ">$ARGV[0].edge";

for (0..$chrom[$all_chromosome]){
	$h=$_;
	for (0..10){
		if ($_==0){
			print OUT1 $array[$h][$_];
		}else{
			print OUT1 "\t".$array[$h][$_];
		}
	}
	print OUT1 "\n";
}			
			
close INPUT;
close OUT1;

#####################################################
#Get every bin based on the edges.#
#####################################################

open IN2,"<$ARGV[0].edge" or die "$!";
my $filename=$ARGV[0].".win$win_size.edge";
my ($prefix,$cenfix,$suffix)=split(/\./,$filename);
open OUT2, ">$prefix.$cenfix.bin";
$n=0;
$m=1;
@array=();
@chrom=();
$chrom[0]=0;
@bin=();
my $ij="";
my $sn=0;
$chromosome="chromosome01";
while (<IN2>) {
	$line = $_;
	chomp($line);
	@hang = split(/\t/,$line);
	$lane_n=9;
	for (0..($lane_n)){
		$array[$n][$_]=$hang[$_];
	}
	$n++;
}
close IN2;
for (0..$n) {
	if ($array[$_][2] eq "$chromosome"){
		$chrom[$m]++;
	}else{
		$chromosome=$array[$_][2];
		$chrom[$m+1]++;
		$m++;
	}
}
for (1..$all_chromosome){
	$chrom[$_]=$chrom[$_]+$chrom[$_-1];
}

#for (1..$all_chromosome){
#	print "\$chrom[$_]=$chrom[$_]\n";
#}

for (1..$all_chromosome){
	$c=$_;
	
	for ($chrom[$c-1]..$chrom[$c]-1) {
		if ($_==$chrom[$c-1]){
			$bin[$sn][0]=$array[$_][2];
			$bin[$sn][1]=1;
			$ij=$array[$_][9];
		}
#		if (($array[$_][9] ne "") && ($array[$_+1][9] ne "") && ($array[$_][9] eq $array[$_+1][9])){
#			print $_."\t".$array[$_][9]."\t".$array[$_+1][9]."\n";
#		}
#		if (($array[$_][9] ne "") && ($array[$_+1][9] eq "") && ($array[$_-1][9] eq "")){
#			$ij=$array[$_][9];
#		}
		if (($array[$_][9] ne "") && ($array[$_+1][9] ne "") && ($array[$_][9] ne $array[$_+1][9])&& ($array[$_][2] eq $array[$_+1][2])){
			$bin[$sn][2]=int(($array[$_][3]+$array[$_+1][3])/2);
			$bin[$sn][3]=$ij;
			if ($ij ne $array[$_][9]){
#				print $array[$_][2]."\t".$array[$_][3]."\n";
			}
			$bin[$sn][4]=$array[$_][1];
			$bin[$sn][5]=$bin[$sn][2]-$bin[$sn][1]+1;
			$sn++;
			$bin[$sn][0]=$array[$_][2];
			$bin[$sn][1]=$bin[$sn-1][2]+1;
			$ij=$array[$_+1][9];
		}
		if ($_==$chrom[$c]-1){
			$bin[$sn][2]=$len_detail[$c];
			$bin[$sn][3]=$ij;
			$bin[$sn][4]="chr_end";
			$bin[$sn][5]=$bin[$sn][2]-$bin[$sn][1]+1;
			$sn++;
		}
	}
}

for (0..($sn-1)){
	$h=$_;
	for (0..5){
		if ($_==0){
			print OUT2 $bin[$h][$_];
		}else{
			print OUT2 "\t".$bin[$h][$_];
		}
	}
	print OUT2 "\n";
}
close OUT2;

#####################################################
# Filtering bins which are smaller than 300K.#
#####################################################
$n=0;
for(1..100){
	$round=$_;
	open IN4,"<$prefix.$cenfix.bin" or die "$!";
	$n=0;
	while (<IN4>) {
		$line = $_;
		@hang = split(/\s+/,$line);
		$lane_n=4;
		for (0..($lane_n)){
			$bin[$n][$_]=$hang[$_];
		}
		$n++;
	}
	close IN4;
	open OUT3, ">$prefix.$cenfix.bin";
	my $skip=0;
	for (0..($n-1)){
		$h=$_+$skip;
		if ($h<$n-1) {
			if (($bin[$h+1][2]-$bin[$h+1][1]+1)>=$min_bin){
				for (0..4){
					if ($_==0){
						print OUT3 $bin[$h][$_];
					}else{
						print OUT3 "\t".$bin[$h][$_];
					}
				}
				$len=$bin[$h][2]-$bin[$h][1]+1;
				print OUT3 "\t".$len;
				print OUT3 "\n";
			}elsif((($bin[$h+1][2]-$bin[$h+1][1]+1)<$min_bin) && ($bin[$h+2][3] eq $bin[$h][3])){
				print OUT3 $bin[$h][0]."\t";
				print OUT3 $bin[$h][1]."\t";
				if (($bin[$h+1][4] ne "chr_end") && ($bin[$h+1][1] ne "1")){
					print OUT3 $bin[$h+2][2]."\t";
					print OUT3 $bin[$h+2][3]."\t";
					print OUT3 $bin[$h+2][4]."\t";
					$len=$bin[$h+2][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					$skip=$skip+2;
				}elsif ($bin[$h+1][1] eq "1"){
					print OUT3 $bin[$h][2]."\t";
					print OUT3 $bin[$h][3]."\t";
					print OUT3 $bin[$h][4]."\t";
					$len=$bin[$h][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
				}elsif ($bin[$h+1][4] eq "chr_end"){
					print OUT3 $bin[$h][2]."\t";
					print OUT3 $bin[$h][3]."\t";
					print OUT3 $bin[$h][4]."\t";
					$len=$bin[$h][2]-$bin[$h][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					
					print OUT3 $bin[$h+1][0]."\t";
					print OUT3 $bin[$h+1][1]."\t";
					print OUT3 $bin[$h+1][2]."\t";
					print OUT3 $bin[$h+1][3]."\t";
					print OUT3 $bin[$h+1][4]."\t";
					$len=$bin[$h+1][2]-$bin[$h+1][1]+1;
					print OUT3 $len;
					print OUT3 "\n";
					$skip=$skip+1;
				}
			}elsif((($bin[$h+1][2]-$bin[$h+1][1]+1)<$min_bin) && ($bin[$h+2][3] ne $bin[$h][3])){
				for (0..4){
					if ($_==0){
						print OUT3 $bin[$h][$_];
					}else{
						print OUT3 "\t".$bin[$h][$_];
					}
				}
				$len=$bin[$h][2]-$bin[$h][1]+1;
				print OUT3 "\t".$len;
				print OUT3 "\n";
			}
		}
		elsif ($h==$n-1) {
			for (0..4){
				if ($_==0){
					print OUT3 $bin[$h][$_];
				}else{
					print OUT3 "\t".$bin[$h][$_];
				}
			}
			$len=$bin[$h][2]-$bin[$h][1]+1;
			print OUT3 "\t".$len;
			print OUT3 "\n";
		}
	}
	close OUT3;
}

#####################################################
#Draw a figure in PNG format indicating bins and SNPs of the RIL.
#####################################################
# GD is optional: binmap_pipeline.py draws the figure itself.
eval { require GD; GD->import(); 1 } or do {
	warn "GD.pm not available; bin calling finished, skip Perl PNG plot.\n";
	exit 0;
};
open OUT4, ">$prefix.$cenfix.combine.png";
my $left = 400;
my $im = new GD::Image(5000+$left,2000); 
my $black = $im->colorAllocate(0,0,0); 
my $white = $im->colorAllocate(255,255,255); 
my $red = $im->colorAllocate(255,0,0);
my $blue = $im->colorAllocate(0,0,255);
my $green = $im->colorAllocate(0,255,0);
my $hetero = $im->colorAllocate(255,215,0);
my $background=$im->colorAllocate(200,200,200);
my $parent1_name = $ARGV[3];
my $parent2_name = $ARGV[4];
$im->fill(10,10,$white); 
$im->rectangle(0+$left,20,5000+$left,30,$black);
$im->fill(2500,25,$black);
$im->rectangle(3600+$left,670,3700+$left,700,$red);
$im->fill(3650+$left,685,$red);
$im->stringFT($red,'/data6/home/zqiang/db/Fonts/helr45w.ttf',20,0,3710+$left,700,"$parent1_name");
#$im->string(gdGiantFont,3550,685,"P1",$red);

$im->rectangle(3600+$left,720,3700+$left,750,$blue);
$im->fill(3650+$left,735,$blue);
$im->stringFT($blue,'/data6/home/zqiang/db/Fonts/helr45w.ttf',20,0,3710+$left,750,"$parent2_name");
#$im->string(gdGiantFont,3550,735,"P2",$blue);

$im->rectangle(3600+$left,770,3700+$left,800,$hetero);
$im->fill(3650+$left,785,$hetero);
$im->stringFT($hetero,'/data6/home/zqiang/db/Fonts/helr45w.ttf',20,0,3710+$left,800,"hetero");
#$im->string(gdGiantFont,3550,785,"P3",$orange);
for (1..50){
	$im->line($_*100+$left,5,$_*100+$left,20,$black);
}
for (1..9){
	my $scale_temp=$_*5;
	my $scale=$scale_temp." Mb";
	$im->stringFT($black,'/data6/home/zqiang/db/Fonts/helr45w.ttf',20,0,$_*500-25+$left,75,"$scale");
	#$im->string(gdGiantFont,$_*500-25,50,"$scale",$black);
}

open INPUT2,"<$ARGV[1]" or die "$!";
while (<INPUT2>) {
	$line1=$_; 
	chomp($line1);
	($chromosome1,$chrom_len)=split(/\s+/,$line1);
	$info=$chromosome1;
	$number1=$chromosome1;
	$number1=~s/chromosome//;
	$number1=~s/^0//;
	$cend=int(($chrom_len/10000)+0.5);
	$im->rectangle(0+$left,$number1*120,$cend+$left,($number1*120+100),$background);
	$im->fill($cend*0.5,($number1*120+50),$background);
	$im->stringFT($black,'/data6/home/zqiang/db/Fonts/helr45w.ttf',30,0,100,($number1*120+50),"$info"); 
}
close INPUT2;

open IN5,"<$prefix.$cenfix.bin" or die "$!";
while (<IN5>) {
	$line=$_; 
	chomp($line);
	($chromosome,$bstart,$bend,$origin,$name)=split(/\s+/,$line);
	$bstart=int($bstart/10000)+$left;
	$bend=int($bend/10000)+$left;
	$number=$chromosome;
	$number=~s/chromosome//;
	$number=~s/^0//;
	if ($origin eq "parent1"){
		$im->rectangle($bstart,$number*120+21,$bend,($number*120+40),$red);
		$im->fill(($bstart+$bend)*0.5,($number*120+30),$red); 
	}elsif ($origin eq "parent2"){
		$im->rectangle($bstart,$number*120+21,$bend,($number*120+40),$blue);
		$im->fill(($bstart+$bend)*0.5,($number*120+30),$blue); 
	}
	elsif ($origin eq "hetero"){
		$im->rectangle($bstart,$number*120+21,$bend,($number*120+40),$hetero);
		$im->fill(($bstart+$bend)*0.5,($number*120+30),$hetero); 
	}
}
close IN5;

open INPUT4,"<$ARGV[0]" or die "$!";
while (<INPUT4>) {
	$line=$_; 
	chomp($line);
	($origin,$name,$chromosome,$bstart,$bend)=split(/\s+/,$line);
	$bstart=int($bstart/10000)+$left;
	$bend=int($bend/10000)+$left;
	$number=$chromosome;
	$number=~s/chromosome//;
	$number=~s/^0//;
	if ($origin eq "P1"){
		$im->rectangle($bstart,$number*120+51,$bend,($number*120+70),$red);
		$im->fill(($bstart+$bend)*0.5,($number*120+60),$red); 
	}elsif ($origin eq "P2"){
		$im->rectangle($bstart,$number*120+71,$bend,($number*120+90),$blue);
		$im->fill(($bstart+$bend)*0.5,($number*120+90),$blue); 
	}
}
close INPUT4;
	binmode OUT4; 
	print OUT4 $im->png; 
close OUT4;
exit; 
