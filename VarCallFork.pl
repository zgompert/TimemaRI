#!/usr/bin/perl
#
# variant calling 
#



use Parallel::ForkManager;
my $max = 10;
my $pm = Parallel::ForkManager->new($max);

open(IN, "blist") or die "failed to read\n";
while(<IN>){
	chomp;
	@line = split(/\s+/,$_);
	$bams = $line[0];
	$genome = $line[1];
	$bams =~ m/bams_([a-z]+)/;
	$base = $1;
	$pm->start and next; ## fork
	system "bcftools mpileup -b $bams -d 500 --ignore-RG -f $genome -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o timema_$base"."_variants.vcf\n";

	$pm->finish;
}

$pm->wait_all_children;



