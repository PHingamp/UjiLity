#!/usr/bin/perl -w

use strict;
use Math::Round;

# ./UjiLity_scan_report.pl 'Uji_v2' Uji_results

my $motif = $ARGV[0];
my $resdir = $ARGV[1];
my $totnb = $ARGV[2] / 100;
my $results = "$resdir/$motif*"; 
my %jobs;
my %max_tm;
my %max;
my %ba_pos; my %ba_range; my %ba_pp;  my %ba_eva; my %ba_dg;
my $positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^E\\s+:"`;
my @positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
#MEGAPRIMER_v1_results//MEGAPRIMER_v1.1/UjiLity_run_output.txt:14:E	: 51
	(my $run, my $job, my $tm) = ($line =~ /_v(\d+\.)(\d+)\/.+:E\s+:\s+(\d+)\s*$/);
	if ( $job ){
		$max_tm{"$run$job"} = "$tm";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}
my %min_tm;
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^e\\s+:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $tm) = ($line =~ /_v(\d+\.)(\d+)\/.+:e\s+:\s+(\d+)\s*$/);
	if ( $job ){
		$min_tm{"$run$job"} = "$tm";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}
my %max_degen;
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^y\\s+:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $tm) = ($line =~ /_v(\d+\.)(\d+)\/.+:y\s+:\s+(\d+)\s*$/);
	if ( $job ){
		$max_degen{"$run$job"} = "$tm";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}
my %max_p3;
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^x\\s+:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $tm) = ($line =~ /_v(\d+\.)(\d+)\/.+:x\s+:\s+(\d+)\s*$/);
	if ( $job ){
		$max_p3{"$run$job"} = "$tm";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}
my %max_neg;
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^u\\s+:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $tm) = ($line =~ /_v(\d+\.)(\d+)\/.+:u\s+:\s+(\d+)\s*$/);
	if ( $job ){
		$max_neg{"$run$job"} = "$tm";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}
my %p3;
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^m\\s+:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $tm) = ($line =~ /_v(\d+\.)(\d+)\/.+:m\s+:\s+(\S+)\s*$/);
	if ( $job ){
		$p3{"$run$job"} = "$tm";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}

$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H "Total nb of seqs in ampl. clades"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){
	(my $run, my $job, my $pos) = ($line =~ /_v(\d+\.)(\d+).+\s(\d+)\s*$/);
	$ba_pos{"$run$job"} = $pos;
	$max{$run} = $job unless exists $max{$run};
	$max{$run} = $job if $job > $max{$run};
}

$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H -P "^n\\s+:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $st, my $en) = ($line =~ /_v(\d+\.)(\d+)\/.+\d+\.(\d+)\.(\d+)\.fna\s*$/);
	if ( $job ){
		$ba_range{"$run$job"} = "$st\t$en";
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
		#warn "ok";
	}
}
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H "Total nb of amplifiable clusters"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $pos) = ($line =~ /_v(\d+\.)(\d+).+\s(\d+)\s*$/);
	$ba_pp{"$run$job"} = $pos;
	$max{$run} = $job unless exists $max{$run};
	$max{$run} = $job if $job > $max{$run};
}
$positive = `find $results -name 'UjiLity_run_output.txt' | xargs grep -n -H "nb double ctg hits:"`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	(my $run, my $job, my $pos) = ($line =~ /_v(\d+\.)(\d+).+nb double ctg hits:\s+(\d+)/);
	$ba_eva{"$run$job"} = $pos;
	$max{$run} = $job unless exists $max{$run};
	$max{$run} = $job if $job > $max{$run};
}
$positive = `find $results/ -name "*.primers.FRrc.fasta" | xargs cat`;
@positives = split /\n/, $positive;
foreach my $line (@positives){#print "$line\n";
	#>BUt53.0_PP1_L DEGEN=1 Tm=53.70
	#GGCCACACCCACATCTAGAATTTAT
	if ((my $run, my $job, my $pos) = ($line =~ /^>\S+_v(\d+\.)(\d+).+DEGEN=(\d+)/)){
		$ba_dg{"$run$job"} = 0 unless exists $ba_dg{"$run$job"};
		$ba_dg{"$run$job"} += $pos;
		$max{$run} = $job unless exists $max{$run};
		$max{$run} = $job if $job > $max{$run};
	}
}
print "JOB\tSTART\tEND\tINCL\tPERC\tDETECTED\tNB_PRIMER_PAIRS\tMIN_TM\tMAX_TM\tMAX_DEGEN\tP3_CONF\tMAX_P3\tMAX_NEG\tSUM_DEGEN\n";
foreach my $run (sort keys %max){
	for (my $j=1;$j<=$max{$run};$j++){
		my $i = "$run$j";
		$ba_range{$i} = '' unless exists $ba_range{$i};
		$ba_pos{$i} = '' unless exists $ba_pos{$i};
		$ba_eva{$i} = '' unless exists $ba_eva{$i};
		$ba_pp{$i} = '' unless exists $ba_pp{$i};
		$ba_dg{$i} = '' unless exists $ba_dg{$i};
		$min_tm{$i} = '' unless exists $min_tm{$i};
		$max_tm{$i} = '' unless exists $max_tm{$i};
		$max_degen{$i} = '' unless exists $max_degen{$i};
		$p3{$i} = '' unless exists $p3{$i};
		$max_p3{$i} = '' unless exists $max_p3{$i};
		$max_neg{$i} = '' unless exists $max_neg{$i};
		if($i !~ /^1.(3|4)/ || $ba_range{$i}){
			print "$i\t$ba_range{$i}\t$ba_pos{$i}".($ba_pos{$i} && $ba_pos{$i}>2500?'*':'')."\t".($ba_pos{$i}?nearest(.01,$ba_pos{$i}/$totnb):'')."\t$ba_eva{$i}\t$ba_pp{$i}\t$min_tm{$i}\t$max_tm{$i}\t$max_degen{$i}\t$p3{$i}\t$max_p3{$i}\t$max_neg{$i}\t$ba_dg{$i}\n";
		}
	}
}
