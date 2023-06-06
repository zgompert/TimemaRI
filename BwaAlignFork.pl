#!/usr/bin/perl
#
# run bwa mem 
#



use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

# knulli
#my $genome = "/uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_knulli/mod_hic_output.fasta";
## green striped 
#my $genome = "/uufs/chpc.utah.edu/common/home/u6000989/data/timema/hic_genomes/t_crist/timema_cristinae_12Jun2019_lu3Hs.fasta";
# chumash
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group1/data/timema/hic_genomes/t_chumash/hic/timema_chumash_29Feb2020_N4ago.fasta";

$in = shift(@ARGV);
open(IN, $in) or die;
while(<IN>){
	chomp;
	m/^(\S+)/;
	$base = $1;
	$fq = "$base".".fq.gz";
	$pm->start and next; ## fork
	system "bwa mem -t 2 -k 19 -r 1.5 -R \'\@RG\\tID:timema-"."$ind\\tPL:ILLUMINA\\tLB:timema-"."$ind\\tSM:timema-"."$ind"."\' $genome $fq > aln_south_$base.sam\n";
	$pm->finish;
}

$pm->wait_all_children;



