#!/usr/bin/perl -w

use strict;
my $inc    = 5;     #5 bp increment for each loop
my $min    = 14400; # scan start 14400
my $max    = 14500; # scan end 15400
my $width  = 290;   # target width
my $margin = 100;   # region for primers
my $analysis_name  = 'MEGAPRIMER_v1';
my $results_folder = 'MEGAPRIMER_v1_results';
my $p3         = 'primer3_run_params.p3';
my $msa_aa     = 'MEGAPRIMER_v1_PolB_Refs_plus_Tara_plus_outgroup_OTV1_non_truncated_prot.aln';
my $msa_nuc    = 'MEGAPRIMER_v1_PolB_Refs_plus_Tara_plus_outgroup_OTV1_non_truncated_nuc.aln';
my $tree       = 'MEGAPRIMER_v1_PolB.nw';
my $outgroup   = 'PHYCO_Ostreococcus_tauri_virus_1'; #optional
my $extra_nucs = 'MEGAPRIMER_v1_PolB_Refs_plus_Tara_plus_outgroup_OTV1_all.fna'; # sequences to test for e-PCR
my $weights    = 'MEGAPRIMER_v1_PolB_Tara_contigs_abundances.csv';  #optional
my $distances  = '-O MEGAPRIMER_v1_PolB.nw.ingroup_phylogram.distances'; #optional
my $cpu             = 1; 
my $degen           = 256;
my $unsuccess_seqs  = 20;     # before closing primer pair, eg 200
my $max_primer_eval = 1000;    #1000
my $cache           = "-c i"; #"-c ifpncbmwzy"
my $t_coffee        = 't_coffee '; # path to t_coffee

my $jobs = 0;

for (my $pos = $min; $pos <= $max ; $pos += $inc){
	my $ps = `ps aux| grep " ./UjiLity.pl" | grep -v grep | wc -l`;
	my $wait = 0;
	while ($ps >= $cpu ){
		print "Waiting...." unless $wait;
		sleep 2;
		$wait += 2;
		$ps = `ps aux| grep " ./UjiLity.pl" | grep -v grep  | wc -l`;
	}
	print " $wait secs\n" if $wait;
	$jobs++;
	print "====Job $jobs window start : $pos\n";
	system("mkdir targets") unless -d "./targets";
	system("mkdir outputs") unless -d "./outputs";
	my $endpos = $pos + $width + 2 * $margin;
	my $ren = $endpos;$ren++;
	my $rst = int($pos / 3);
	$rst = 1 unless $rst;
	$ren = int($ren / 3);
	system("$t_coffee -other_pg seq_reformat -in $msa_aa  -action +extract_block cons $rst $ren     > targets/$analysis_name.$jobs.$pos.$endpos.faa");
	system("$t_coffee -other_pg seq_reformat -in $msa_nuc -action +extract_block cons $pos $endpos  > targets/$analysis_name.$jobs.$pos.$endpos.fna");
	system("./UjiLity.pl -k bottomup -d $results_folder/$analysis_name.$jobs -t $tree $distances -p targets/$analysis_name.$jobs.$pos.$endpos.faa -n targets/$analysis_name.$jobs.$pos.$endpos.fna" .
	       " -m $p3 -z $extra_nucs ".($weights?"-a $weights":'')." -s target -y $degen -x $max_primer_eval -i 0 -j 100 -l 0 -q prealigned -h $margin-".($margin+$width).
	       " -w 0 -e 45 -E 51 -u $unsuccess_seqs -D \"+\" $cache ".($outgroup?"-o $outgroup":'')." > outputs/$analysis_name.$jobs \&");
}

