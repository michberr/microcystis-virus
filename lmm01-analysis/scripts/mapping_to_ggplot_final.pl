#!/usr/bin/perl
use strict;
use autodie;
use Getopt::Long;
use Cwd;
my $h='';
my $dir='';
my $genome_asked='';
my $out_file='';
my $cmd='';
my $out='';
GetOptions ('help' => \$h, 'h' => \$h, 'b=s'=>\$dir, 'g=s'=>\$genome_asked, 'o=s'=>\$out_file);
if ($h==1 || $dir eq "" || $out_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to cluster new contigs
# Arguments :
-b : bam file directory
-o : out file
-g : (facultative) Genome Id to generate the RP for (leave blank for all genomes available)\n";
	die "\n";
}

my $default_th_TOV=0; # Th on the number of reads mapped for a sample to get out of the "Other TOV" category
my %th_reads_TOV=("PSAHP1_revcom_reperm_ORF00031"=>50); #  We used 50 by default, 500 for HP1
# my $palette_file="Color_palette.tsv";
my $code_other="Other_TOV90";


my $k=0;
my %counts;
my %check_seq;
my %store_col;
my @list=<$dir*_sorted.bam>;
open my $s1,">",$out_file;
print $s1 "Genome,Start,Stop,Id,Sample\n";
foreach my $bam_file (@list){
	#if ($bam_file=~/All_TOV_90/){next;} # We don't want a recruitment from the merged bam files of all samples
	$bam_file=~/(\w+)-bam\/(\w+)-(\w+)-align.bam_sorted.bam/;
	my $sample_name=$2;
	print "$bam_file\t$sample_name\n";
	my %store_hit=();
	my %count=();
	open my $bam,"samtools view $bam_file |";
	while(<$bam>){
		chomp($_);
		my @tab=split("\t",$_);
		if (($tab[2] eq $genome_asked) || ($genome_asked eq "")){
			my $genome=$tab[2];
			my $start=$tab[3];
			my $length=length($tab[9]);
			my $stop=$start+$length;
			my $match=0;
			my $mismatch=0;
			my $md_index=11;
			for (my $j=11;$j<=$#tab;$j++){
				if ($tab[$j]=~/^MD/){$md_index=$j;}
			}
			$tab[$md_index]=~/MD:[^:]+:(.*)/;
			$tab[$md_index]=$1;
			my $store_tmp_md=$1;
			while (length($tab[$md_index])>0){
				if ($tab[$md_index]=~/^\D/){
					while($tab[$md_index]=~/^\D(.*)/){
						$tab[$md_index]=$1;
						$mismatch++;
					}
				}
				elsif ($tab[$md_index]=~/^(\d+)(.*)/){
					$match+=$1;
					$tab[$md_index]=$2;
				}
				else{
					print "I can't deal with $tab[$md_index] !!! - $_\n";
					<STDIN>;
				}
			}
			my $ali_length=$match+$mismatch;
# 			print "$tab[$md_index]\t$match\t$mismatch\n";
			if (length($tab[9])>$ali_length){$ali_length=length($tab[9]);}
			my $id=int ($match/($match+$mismatch)*10000)/100;
			$count{$genome}++;
			$store_hit{$genome}{"Hit_".$count{$genome}}="\"$genome\",$start,$stop,$id";
			}
		}
		close $bam;
		my $tag=0;
		foreach my $genome (keys %count){
			my $th=$default_th_TOV;
			my $sample_name_g=$sample_name;
			if (defined($th_reads_TOV{$genome})){$th=$th_reads_TOV{$genome}}
			if ($count{$genome}<$th){
				$sample_name_g=$code_other;
			}
			else{
				$tag=1;
			}
			if ($count{$genome}>0){
				foreach(keys %{$store_hit{$genome}}){
					print $s1 "$store_hit{$genome}{$_},$sample_name_g\n";
				}
			}
			print "$sample_name_g - $genome - $count{$genome}\n";
		}
		
}
close $s1;

