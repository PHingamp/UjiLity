#!/usr/bin/perl -w

=pod
    		UjiLity - UJI muLtiple prImer uTilitY
    A script to design sets of oligonucleotide primers with controlled
    degeneracy aimed at amplifying marker genes with very high sequence
    diversity (e.g. virus diversity studies).
    
    The script was written in the spring of 2015 in Uji, Japan, by 
    Pascal Hingamp as a Kyoto University visiting scholar in the lab
    of Hiro Ogata whose team members provided instrumental insights during
    development: Hiro Ogata, Takashi Yoshida, Susumu Goto, Tomoko Mihara
    & Yosuke Nishimura. The resulting primer design strategy was validated 
    by an experimental test of a UjiLity designed set of primers (called
    MEGAPRIMER) that demonstrated detection of giant viruses diversity
    in marine water: https://www.mdpi.com/1999-4915/10/9/496
    
    		Version 2.17 (2022-05-28) pascal.hingamp@univ-amu.fr
    		
    Warning: as often, this script was never meant to grow so large.
    It therefore rather monolithic, severely lacking in structure,
    as well as begging for richer comments... It does however do the job.
    The STDOUT output is very terse, basically showing progress, whilst
    the STDERR output is over verbose and should be redirected to a file.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
=cut



use strict 'vars';
use vars qw($opt_t $opt_o $opt_O $opt_p $opt_n $opt_d $opt_m $opt_b $opt_c $opt_s $opt_g $opt_e $opt_E $opt_a $opt_z $opt_y $opt_u $opt_k $opt_x $opt_i $opt_j $opt_f $opt_r $opt_l $opt_w $opt_q $opt_h $opt_S $opt_T $opt_D);
use File::Basename;
use FileHandle;
use Getopt::Std;
use File::Copy qw(copy);
use IO::Handle;
use POSIX qw(strftime);
use Parallel::ForkManager;
use Math::Round;
use constant VERSION => '2.17';
STDERR->autoflush(1);
my $debug = 0;
my $seqkit = 'seqkit --threads 1 ';
my $t_coffee = 't_coffee ';
my $mafft = 'mafft --quiet --thread 1 '; 
my $blastx = 'blastx -num_threads 1 ';
my $primer3_exec = 'primer3_core -strict_tags ';
my $dnaMate = './bin/dnaMATE ';
my $dreg = 'fuzznuc -supper1 -auto -rformat excel ';
my $seqret = 'seqret ';
my $goalign = './bin/goalign_amd64_linux ';
my $nw_utils = './bin';
my ($sm,$pm,$col,$max_ambig,$max_degen,$leaf_primers,$ti,$ppnb,$runID) = (0,0,0,0,0,0,0,0,'');
my %stats;my $manager;my ($css,$orn,$compiled_primers) = ('','',''); my %oligo_tm_cache;
my (@colors, %pp2color, %amb_tbl, %positives,%envHits,%envCounts,%abundances,%dist,%todo,%primer_listing,@ingroupNames);

init();
if ($opt_k eq 'topdown'){
	get_primers_topdown("$opt_d/$opt_t","ingroup");	
	stats();
	evaluate_primers();
	wrap_up();
}elsif($opt_k eq 'bottomup'){
	get_primers_bottomup();
	stats();
	evaluate_primers();
	wrap_up();
}elsif($opt_k eq 'evaluate'){
	evaluate_primers();
	wrap_up();
}elsif($opt_k eq 'degentm'){
	calculate_tm();
}
close(STDERR);
############################################################################################################################
sub divide_tree{
	my $tree = shift;
	my @leaves = get_leaves($tree);
	# if tree to divide has one leaf, then reached end, stop dividing.
	if (scalar(@leaves) < 2) {
		print STDERR "========================================================\n".
			     "Reached the leaf, can't divide further!\n";
		return;
	}
	print STDERR "========================================================\n".
		     "Dividing: $tree\n";# . `$nw_utils/nw_display -e r $tree`."\n";
	my %dists;
	# get distances from the root of each tree node
	foreach my $line (split /\n/, `$nw_utils/nw_distance -s l -n $tree`){
		my($node,$dist) = split /\t/, $line;
		$dists{$node}=$dist;
		$dists{$node}+=1 unless $node =~ /^N\d+$/;#dirty hack to lengthen the leaves versus internal nodes
						    #so that internal nodes are chosen first in case of equal
						    #distance from root
	}
	# sort distances from smallest to greatest; the first is the subclade node itself
	# so ignore it, the second node is the first subclade "A" after division
	my ($null, $min, $min2, $min3) = (sort { $dists{$a} <=> $dists{$b} } keys %dists);
	$debug && print STDERR "The closest nodes to the root are: $min & $min2! Then $min3\n";
	do_cmd("$nw_utils/nw_clade $tree $min > $tree.$min");
	# the -s switch in the next command returns the sister clade of clade "A" above
	# ie the second clade "B" after division
	do_cmd("$nw_utils/nw_clade -s $tree $min > $tree.$min.sis");
	#do_cmd("$nw_utils/nw_clade $tree $min2 > $tree.$min2");
	# to get the name of the sister node "B", it is the closest to root within its own subclade "B"
	%dists=();
	foreach my $line (split /\n/, `$nw_utils/nw_distance -s l -n $tree.$min.sis`){
		my($node,$dist) = split /\t/, $line;
		$dists{$node}=$dist;
	}
	my ($min_sis) = (sort { $dists{$a} <=> $dists{$b} } keys %dists);
	$debug && print STDERR "Sister node to $min is: $min_sis!\n";
	rename("$tree.$min.sis","$tree.$min_sis") || die "can't rename $tree.$min.sis $tree.$min_sis $!";
	# search primers for the two descendant subclades; if they don't succeed, they will recursively 
	# call "divide_tree" on the subtrees... 
	my ($L1,$L2);
	if ($opt_S eq 'cladeSize'){
		$L1 = scalar(get_leaves("$tree.$min"));
		$L2 = scalar(get_leaves("$tree.$min_sis"));
	}elsif ($opt_S eq 'abundance'){
		$L1=$L2=0;
		foreach my $id (get_leaves("$tree.$min")){
			$L1 += $abundances{$id};
		}
		foreach my $id (get_leaves("$tree.$min_sis")){
			$L2 += $abundances{$id};
		}
	}else{
		#order by node numbering
		($L1) = ($min_sis =~ /(\d+)$/);
		($L2) = ($min  =~ /(\d+)$/);
		
	}
	if ($L1 > $L2){
		get_primers_topdown($tree,$min);
		get_primers_topdown($tree,$min_sis);
	}else{
		get_primers_topdown($tree,$min_sis);
		get_primers_topdown($tree,$min);
	}
}
############################################################################################################################
sub get_primers_topdown{
	my $parent = shift;
	my $current_node = shift;
	my $tree = "$parent.$current_node";
	#my @sequences = @ingroupNames;
	#my($nodename, $dirs, $suffix) = fileparse($tree);
	#my $code = $nodename;
	#$code =~ s/^.+\.ingroup\.N/N/;
	my $original_current_node = $current_node;
	unless ($current_node =~ /^N\d+/) { 
		$current_node = "L".($leaf_primers++);
	}
	my @leaves_orig = get_leaves($tree);
	my @leaves;
	$stats{'clade_total'}++;
	#my $ppnb = 0;
	print STDERR "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n".
		     "Getting primers for subtree $original_current_node (".scalar(@leaves_orig)." leaves):\n" . `$nw_utils/nw_display $tree`."\n";
	my $pruned = 0;
	foreach my $id (@leaves_orig){
			push (@leaves, $id) unless exists($positives{$id});
			$pruned++ if exists($positives{$id});
	}
	print STDERR " Pruned $pruned leaves which were already detected by previous primers.\n";
	if (scalar(@leaves) > 1){
		my $nucs = sub_fasta($tree,$opt_n,"$tree.fna",\@leaves);
		if ($opt_q eq 'align'){	
			#align protein sequences
			my $prots = sub_fasta($tree,$opt_p,"$tree.faa",\@leaves);
			#system("mv $tree.aln.is.done $prots.aln.p.cached") if -e "$tree.aln.is.done";
			do_cmd_cached("$mafft --maxiterate 1000 --$opt_g $prots > $prots.aln","$prots.aln",'p');
			#thread nuc seqs in prot aln
			pal2nal("$prots.aln", $nucs, '-output','clustal', '-outfile', "$nucs.aln");
		}else{
			#extract subset from prealigned nuc file
			#do_cmd("t_coffee -other_pg seq_reformat -in $nucs -action +rm_gap 100 -output clustal > $nucs.aln");
			#do_cmd("$t_coffee -other_pg seq_reformat -in $nucs -output clustal > $nucs.aln");
			#do_cmd("$seqret -osformat2 clustal -stdout -auto -sequence $nucs > $nucs.aln");
			do_cmd("$goalign reformat clustal -i $nucs > $nucs.aln");
		}
		#make consensus
		if(cached("$tree.cons",'c') && cached("$tree.cons.mask",'c')){
			$debug && print STDERR "\nUsing cached: $tree.cons & $tree.cons.mask\n";
		}else{
			consensus("$nucs.aln",$opt_j,"$tree.cons","$tree.cons.mask");
			do_cmd("touch $tree.cons.c.cached");do_cmd("touch $tree.cons.mask.c.cached");
		}
		#do_cmd_cached("$consensus -a $nucs.aln -t 95 -o $tree.cons","$tree.cons",'c');
	}elsif(scalar(@leaves) == 1){
		my $nucs = sub_fasta($tree,$opt_n,"$tree.cons",\@leaves);
		do_cmd_cached("touch $tree.cons.mask","$tree.cons.mask",'c');
	}else{
		print STDERR "This tree is empty after pruning, finished:)\n";
		#$stats{'clade_positives'}++;
		return;
	}
	my %domst;my %domen;
	if (defined $opt_b){
		#locate target regions with blast
		do_cmd_cached("$blastx -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -query $tree.cons -out $tree.blastx -db $opt_b","$tree.blastx",'b');
		open BLAST, "$tree.blastx" || die "$!";

		while (<BLAST>){
			my @cols = split /\t/;
			$cols[1] =~ s/_\d+$//;
			if (!exists($domst{$cols[1]})){
				$domst{$cols[1]}= $cols[6];# start pos
				$domen{$cols[1]}= $cols[7];# end pos
				$debug && print STDERR "Found Domain $cols[1]: $domst{$cols[1]}-$domen{$cols[1]}\n";
			}
		}
		close BLAST;
	}elsif(defined $opt_h && $opt_h){
		if(my ($st,$en) = ($opt_h =~ /^(\d+)\-(\d+)/)){
			$domst{'Dom1'}= $st;# start pos
			$domen{'Dom1'}= $en;# end pos
			$debug && print STDERR "Using Domain Dom1: $st-$en\n";
		}else{
			die "Expecting -h in the form '123-456' !";
		}
	}else{
		die "Neither -b nor -h were provided !";
	}
	my ($primers_pairs_found_ref,$primers_found_ref);
	#don't even bother if less than 1 domain have been detected by BLAST?
	if (scalar(keys %domst) < 1){
		print STDERR "**** Not even 1 domain detected by BLASTx, not even bothering to predict primers! ****\n"
	}else{
		($primers_pairs_found_ref,$primers_found_ref) = compute_primers_Primer3Plus($tree,$current_node,\%domst,\%domen);
	}
	if (scalar(keys %$primers_pairs_found_ref) > 0){
	#if primer pairs found, stop
		my %matef;my %mater;my %mates;
		$ppnb++;
		my $ppcode=$runID.'_'.$original_current_node."_PP$ppnb";
		print STDERR "\n______________________________________________________________\n".
			     "Found ".scalar(keys %$primers_pairs_found_ref).
			     " domains with a possible PCR, finished this clade.\nPrimer pair $ppcode\n".
			     "______________________________________________________________\n";  
		foreach my $id (@leaves){
			$positives{$id} = $ppcode;
		}
		$primer_listing{$ppcode}++;
		my $primer_r=$primers_found_ref->{'PRIMER_RIGHT.SEQ'};
		my $primer_f=$primers_found_ref->{'PRIMER_LEFT.SEQ'};
		$compiled_primers.=$primer_r.$primer_f;
		my $fastar .= ">$ppcode"."_R DEGEN=" .$primers_found_ref->{'PRIMER_RIGHT.DEGEN'}.' Tm='.$primers_found_ref->{'PRIMER_RIGHT.TM'}."\n$primer_r\n";
		my $fastaf .= ">$ppcode"."_L DEGEN=" .$primers_found_ref->{'PRIMER_LEFT.DEGEN'} .' Tm='.$primers_found_ref->{'PRIMER_LEFT.TM'} ."\n$primer_f\n";
		primers($fastar,'R');
		primers($fastaf,'F');
		print STDERR "$fastaf$fastar";
		my $prc = `echo \"$primer_r\" | revseq -filter`;
		($prc) = ($prc =~ /:\n(\S+)/);
		$primer_f = pattern($primer_f);
		$prc = pattern($prc);
		do_cmd_cached("$dreg -sequence $opt_d/$opt_t.ingroup.fna -pattern \"$primer_f\" -outfile $tree.dreg.f","$tree.dreg.f",'w');
		do_cmd_cached("$dreg -sequence $opt_d/$opt_t.ingroup.fna -pattern \"$prc\" -outfile $tree.dreg.r","$tree.dreg.r",'w');
		open IN,"$tree.dreg.f" or die "$!";
		while (my $line = <IN>){
			chomp($line);
			my @cols = split (/\t/, $line);
			next if $cols[0] eq 'SeqName';
			$matef{$cols[0]}=$cols[1];
			$mates{$cols[0]}=1;
		}
		close IN;
		open IN,"$tree.dreg.r" or die "$!";
		while (my $line = <IN>){
			chomp($line);
			my @cols = split (/\t/, $line);
			next if $cols[0] eq 'SeqName';
			$mater{$cols[0]}=$cols[1];
			$mates{$cols[0]}=1;
		}
		close IN;
		my @are_hit;
		foreach my $id (keys %mates){
			if (exists($mater{$id}) && exists($matef{$id}) && $mater{$id} > 0 && $matef{$id} > 0 && $mater{$id} - $matef{$id} < 700){ #TODO make this as script parameter!
				push @are_hit,$id;
			 }
		}
		my @extra;
		if (scalar(@are_hit) < scalar(@leaves)){
			print STDERR "OOPS, DREG only found ".scalar(@are_hit)." of the ".scalar(@leaves)." targets...\n";
		}elsif(scalar(@are_hit) > scalar(@leaves)){
			print STDERR "WOW, DREG found ".scalar(@are_hit).", more than the ".scalar(@leaves)." targets...\n";
			foreach my $id (@are_hit){
				next if (grep {$_ eq $id} @leaves );
				if (exists($positives{$id})){
					$debug && print STDERR "\t $id was already successfully targeted...\n";
					push @extra, $id;
				}else{
					$debug && print STDERR "\t $id had not yet been targeted, adding to this Primer Pair...\n";
					$positives{$id} = $ppcode;
					push @leaves, $id;
				}
			}
		}
		$col+=1;$col=0 if $col >= scalar(@colors);
		$pp2color{$ppcode}=$col;
		add_css('Clade',$original_current_node);
		foreach my $id (@leaves){
			print STDERR "OUCH, Intended target $id was not found by DREG !?!\n" unless grep {$_ eq $id} @are_hit;
		}
		$stats{'clade_positives'}++;
		$stats{'clade_positives_leaves_total'} += scalar(@leaves);
		$stats{'clade_positives_list'} .= ($stats{'clade_positives_list'}?'|':'').join('|',@leaves);
		$envHits{$ppcode}=0 unless exists($envHits{$ppcode});
		$envCounts{$ppcode}=0 unless exists($envCounts{$ppcode});
		summary($stats{'clade_positives'}."\t$ppcode\t".scalar(@leaves)."\t".
			$envHits{$ppcode}."\t".$envCounts{$ppcode}."\t".
			scalar(keys %$primers_pairs_found_ref)."\t".join(',',(sort keys %$primers_pairs_found_ref))."\t".
			$primers_found_ref->{'amplicon_length'}."\t".$primers_found_ref->{'PRIMER_LEFT.TM'}."\t".
			$primers_found_ref->{'PRIMER_LEFT.DEGEN'}."\t".$primers_found_ref->{'PRIMER_RIGHT.TM'}."\t".
			$primers_found_ref->{'PRIMER_RIGHT.DEGEN'}.
			"\t$ppcode".'_L:'.$primers_found_ref->{'PRIMER_LEFT.SEQ'}."\t$ppcode".'_R:'.
			$primers_found_ref->{'PRIMER_RIGHT.SEQ'}."\t".
			$primers_found_ref->{'left'}."\t".
			$primers_found_ref->{'right'}."\t".
			$primers_found_ref->{'pairs'}."\t".
			join(',',@leaves)."\t".
			join(',',@extra)."\n"
			);
			print STDOUT scalar(@leaves)." sequences detected, ".(scalar(@ingroupNames)-$stats{'clade_positives_leaves_total'})." sequences to go...\n";
	}else{
		#if not, divide clade further into subclades
		print STDERR "______________________________________________________________\n".
			     "No suitable primers found, let's try subclades...\n".
			     "______________________________________________________________\n";
		$stats{'clade_negatives'}++;
		$stats{'clade_negatives_leaves_total'} += scalar(@leaves);
		$stats{'clade_negatives_list'}.= ($stats{'clade_negatives_list'}?'|':'').join('|',@leaves) if scalar(@leaves)==1;
		#summary("$stats{'clade_negatives'}\t$current_node\t".(scalar(@leaves)>1?'no_internal_node':'no_terminal_leaf').
		#	"\t".scalar(@leaves)."\t".join(',',@leaves)."\n",'.reject');
		divide_tree($tree);
	}
}

############################################################################################################################

sub get_primers_bottomup{
	my $tree = "$opt_d/$opt_t.ingroup";
	#my @sequences = @ingroupNames;
	foreach my $ss (@ingroupNames){
		$todo{$ss} = 1;
	}	
	#align protein sequences
	#my $allprots = sub_fasta($tree,$opt_p,"$tree.faa",\@ingroupNames);
	#do_cmd_cached("$mafft --maxiterate 1000 --$opt_g $allprots > $allprots.aln","$allprots.aln",'p');
	#thread nuc seqs in prot aln
	my $allnucs = sub_fasta($tree,$opt_n,"$tree.fna",\@ingroupNames);
	#pal2nal("$allprots.aln", $allnucs, '-output','clustal', '-outfile', "$allnucs.aln");
	#my $ppnb = 0;
	foreach my $seq (@ingroupNames){$abundances{$seq} = 1 unless exists($abundances{$seq})}
	#sort sequences by abundance
	my @sorted_by_abundance = sort {$abundances{$b} <=> $abundances{$a} or $a cmp $b} @ingroupNames;
	print STDOUT "=======================================================\n";
	print STDOUT scalar(@ingroupNames)." sequences to target for amplification.\n";
	my $singleSeqFail = 0;
	SEQ: foreach my $seqid (@sorted_by_abundance) {
		print STDERR "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		if ($positives{$seqid}){
			print STDERR "Sequence $seqid is already detected by Primer Pair $positives{$seqid}. Go to next most abundant sequence.\n".
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
			next SEQ;
		}
		if ($ppnb > 0){
			my $ppcodep=$runID."_PP$ppnb";
			if (! -f "$opt_d/$opt_t.".$ppcodep."_arch.tar.gz"){
				do_cmd("tar cfz $opt_d/$opt_t.".$ppcodep."_arch.tar.gz $opt_d/$opt_t.$ppcodep.*");
				do_cmd("rm $opt_d/$opt_t.$ppcodep.*");
			}
		}
		$ppnb++;
		my $ppcode=$runID."_PP$ppnb";

		print STDOUT "-------------------------------------------------------\n";
		print STDOUT "Primer pair $ppcode\n";
		$stats{'clade_total'}++;
		print STDERR "Cluster $ppcode to be seeded with $seqid (abundance: $abundances{$seqid}).\n".
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		my @still_undetected = grep {!exists($positives{$_})} (@ingroupNames);
		#sort remaining sequences by order of proximity in guide tree
		my @totry = sort {$dist{$seqid}{$a} <=> $dist{$seqid}{$b} or $a cmp $b} (@still_undetected);
		my (@leaves,@prec_leaves);
		my ($primers_pairs_found_ref,$primers_found_ref,$prev_primers_pairs_found_ref,$prev_primers_found_ref);
		my $iter = 0;my $never_give_up=0;
		my $ign=0;
		AGGREGATE: while (1){
			my $add = shift @totry;
			if (exists($positives{$add})){#shouldn't happen because @still_undetected = grep {!exists($positives{$_})} (@ingroupNames);
				$ign++;
				next AGGREGATE;
			}
			print STDERR "======================================================================================================\n";
			print STDERR "Ignoring $ign closest homologs which already have been detected \n";
			if ($iter > 0 && ! -f $tree."_arch.tar"){
				do_cmd("tar cf ".$tree."_arch.tar $tree.*");
				do_cmd("rm $tree.*");
			}
			$iter++;
			$tree = "$opt_d/$opt_t.$ppcode.$iter";
			print STDERR "Adding to cluster $ppcode the next closest homolog N°$iter: $add (abundance: $abundances{$add}).\n";
			push (@leaves, $add);
			if (scalar(@leaves) > 1){
				my $nucs = sub_fasta($tree,$opt_n,"$tree.fna",\@leaves);
				if ($opt_q eq 'align'){	
					#align protein sequences
					my $prots = sub_fasta($tree,$opt_p,"$tree.faa",\@leaves);
					#system("mv $tree.aln.is.done $prots.aln.p.cached") if -e "$tree.aln.is.done";
					do_cmd_cached("$mafft --maxiterate 1000 --$opt_g $prots > $prots.aln","$prots.aln",'p');
					#thread nuc seqs in prot aln
					pal2nal("$prots.aln", $nucs, '-output','clustal', '-outfile', "$nucs.aln");
				}else{
					#extract subset from prealigned nuc file
					#do_cmd("t_coffee -other_pg seq_reformat -in $nucs -action +rm_gap 100 -output clustal > $nucs.aln");
					#do_cmd("$t_coffee -other_pg seq_reformat -in $nucs -output clustal > $nucs.aln");
					#do_cmd("$seqret -osformat2 clustal -stdout -auto -sequence $nucs > $nucs.aln");
					do_cmd("$goalign reformat clustal -i $nucs > $nucs.aln");
				}				
				#make consensus
				if(cached("$tree.cons",'c') && cached("$tree.cons.mask",'c')){
					$debug && print STDERR "\nUsing cached: $tree.cons & $tree.cons.mask\n";
				}else{
					consensus("$nucs.aln",$opt_j,"$tree.cons","$tree.cons.mask");
					do_cmd("touch $tree.cons.c.cached");do_cmd("touch $tree.cons.mask.c.cached");
				}
				#do_cmd_cached("$consensus -a $nucs.aln -t 95 -o $tree.cons","$tree.cons",'c');
			}elsif(scalar(@leaves) == 1){
				my $nucs = sub_fasta($tree,$opt_n,"$tree.cons",\@leaves);
				do_cmd_cached("touch $tree.cons.mask","$tree.cons.mask",'c');
			}else{
				print STDERR "This tree is empty after pruning, finished:)\n";
				#$stats{'clade_positives'}++;
				next SEQ;
			}
			my %domst;my %domen;
			if (defined $opt_b){
				#locate target regions with blast
				do_cmd_cached("$blastx -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -query $tree.cons -out $tree.blastx -db $opt_b","$tree.blastx",'b');
				open BLAST, "$tree.blastx" || die "$!";

				while (<BLAST>){
					my @cols = split /\t/;
					$cols[1] =~ s/_\d+$//;
					if (!exists($domst{$cols[1]})){
						$domst{$cols[1]}= $cols[6];# start pos
						$domen{$cols[1]}= $cols[7];# end pos
						$debug && print STDERR "Found Domain $cols[1]: $domst{$cols[1]}-$domen{$cols[1]}\n";
					}
				}
				close BLAST;
			}elsif(defined $opt_h && $opt_h){
				if(my ($st,$en) = ($opt_h =~ /^(\d+)\-(\d+)/)){
					$domst{'Dom1'}= $st;# start pos
					$domen{'Dom1'}= $en;# end pos
					$debug && print STDERR "Using Domain Dom1: $st-$en\n";
				}else{
					die "Expecting -h in the form '123-456' !";
				}
			}else{
				die "Neither -b nor -h were provided !";
			}
			#don't even bother if less than 1 domain have been detected by BLAST?
			my ($test_primers_pairs_found_ref,$test_primers_found_ref);
			if (scalar(keys %domst) < 1){
				print STDERR "**** Not even 1 domain detected by BLASTx, not even bothering to predict primers! ****\n"
			}else{
				($test_primers_pairs_found_ref,$test_primers_found_ref) = compute_primers_Primer3Plus($tree,$ppcode,\%domst,\%domen);
			}
			print STDERR "\n======================================================================================================\n";
			if (scalar(keys %$test_primers_pairs_found_ref) > 0 && scalar (@totry) > 0){
				print STDERR "=====> primers found, so try adding more sequences to current cluster...\n";
				($prev_primers_pairs_found_ref,$prev_primers_found_ref,@prec_leaves) = ($test_primers_pairs_found_ref,$test_primers_found_ref,@leaves);
				$never_give_up=0;
				next AGGREGATE;
			}elsif(scalar(keys %$test_primers_pairs_found_ref) > 0 && scalar (@totry) == 0){
				print STDERR "=====> primers found, no more sequences to add to current cluster, so wrap up this cluster!\n";
				($primers_pairs_found_ref,$primers_found_ref) = ($test_primers_pairs_found_ref,$test_primers_found_ref);
				last AGGREGATE;
			}elsif(scalar(keys %$test_primers_pairs_found_ref) == 0 && $iter == 1){
				print STDERR "=====> primers not found even with only the single seed sequence, abandon this seed...\n";
				@leaves = grep {$_ ne $add} @leaves;
				$ppnb--;
				$singleSeqFail++;
				last AGGREGATE;
			}elsif(scalar(keys %$test_primers_pairs_found_ref) == 0 && scalar (@totry) > 0 && $never_give_up < $opt_u){
				print STDERR "=====> primers not found, but try adding next ".$never_give_up."th similar sequence...\n";
				@leaves = grep {$_ ne $add} @leaves;
				$never_give_up++;
				next AGGREGATE;
			}else{
				print STDERR "OK, tried $opt_u further sequences with no success, giving up hope that any more will work!\n" if $never_give_up == $opt_u;
				print STDERR "=====> primers not found, no more sequences to add, so revert to previous successful cluster cluster & end!\n";
				($primers_pairs_found_ref,$primers_found_ref,@leaves) = ($prev_primers_pairs_found_ref,$prev_primers_found_ref,@prec_leaves);
				last AGGREGATE;
			}
		}
		#############################################################################################
		#if primer pairs found, stop
		if (scalar(keys %$primers_pairs_found_ref) > 0){
			my %matef;my %mater;my %mates;
			foreach my $id (@leaves){
				$positives{$id} = $ppcode;
			}
			print STDERR "\n________________________________________________________________________\n".
				     "Found ".scalar(keys %$primers_pairs_found_ref).
				     " domains with a possible PCR, new Primer Pair: $ppcode.\n".
				     "________________________________________________________________________\n\n";
			$primer_listing{$ppcode}++;
			my $primer_r=$primers_found_ref->{'PRIMER_RIGHT.SEQ'};
			my $primer_f=$primers_found_ref->{'PRIMER_LEFT.SEQ'};
			$compiled_primers.=$primer_r.$primer_f;
			my $fastar .= ">$ppcode"."_R DEGEN=" .$primers_found_ref->{'PRIMER_RIGHT.DEGEN'}.' Tm='.$primers_found_ref->{'PRIMER_RIGHT.TM'}."\n$primer_r\n";
			my $fastaf .= ">$ppcode"."_L DEGEN=" .$primers_found_ref->{'PRIMER_LEFT.DEGEN'} .' Tm='.$primers_found_ref->{'PRIMER_LEFT.TM'} ."\n$primer_f\n";
			primers($fastar,'R');
			primers($fastaf,'F');
			print STDERR "$fastaf$fastar";
			my $prc = `echo \"$primer_r\" | revseq -filter`;
			($prc) = ($prc =~ /:\n(\S+)/);
			$primer_f = pattern($primer_f);
			$prc = pattern($prc);
			do_cmd_cached("$dreg -sequence $opt_d/$opt_t.ingroup.fna -pattern \"$primer_f\" -outfile $tree.dreg.f","$tree.dreg.f",'w');
			do_cmd_cached("$dreg -sequence $opt_d/$opt_t.ingroup.fna -pattern \"$prc\" -outfile $tree.dreg.r","$tree.dreg.r",'w');
			open IN,"$tree.dreg.f" or die "$!";
			while (my $line = <IN>){
				chomp($line);
				my @cols = split (/\t/, $line);
				next if $cols[0] eq 'SeqName';
				$matef{$cols[0]}=$cols[1];
				$mates{$cols[0]}=1;
			}
			close IN;
			open IN,"$tree.dreg.r" or die "$!";
			while (my $line = <IN>){
				chomp($line);
				my @cols = split (/\t/, $line);
				next if $cols[0] eq 'SeqName';
				$mater{$cols[0]}=$cols[1];
				$mates{$cols[0]}=1;
			}
			close IN;
			my @are_hit;
			foreach my $id (keys %mates){
				if (exists($mater{$id}) && exists($matef{$id}) && $mater{$id} > 0 && $matef{$id} > 0 && $mater{$id} - $matef{$id} < 1500){ #TODO make this as script parameter!
					push @are_hit,$id;
				 }
			}
			my @extra;
			if (scalar(@are_hit) < scalar(@leaves)){
				print STDERR "OOPS, DREG only found ".scalar(@are_hit)." of the ".scalar(@leaves)." targets...\n";
			}elsif(scalar(@are_hit) > scalar(@leaves)){
				print STDERR "WOW, DREG found ".scalar(@are_hit).", more than the ".scalar(@leaves)." targets...\n";
				foreach my $id (@are_hit){
					next if (grep {$_ eq $id} @leaves );
					if (exists($positives{$id})){
						$debug && print STDERR "\t $id was already successfully targeted...\n";
						push @extra, $id;
					}else{
						$debug && print STDERR "\t $id had not yet been targeted, adding to this Primer Pair...\n";
						$positives{$id} = $ppcode;
						push @leaves, $id;
					}
				}
			}
			$col+=1;$col=0 if $col >= scalar(@colors);
			$pp2color{$ppcode}=$col;
			foreach my $id (@leaves){
				print STDERR "OUCH, Intended target $id was not found by DREG !?!\n" unless grep {$_ eq $id} @are_hit;
				add_css('Individual',$id);
			}
			$stats{'clade_positives'}++;
			$stats{'clade_positives_leaves_total'} += scalar(@leaves);
			$stats{'clade_positives_list'} .= ($stats{'clade_positives_list'}?'|':'').join('|',@leaves);
			$envHits{$ppcode}=0 unless exists($envHits{$ppcode});
			$envCounts{$ppcode}=0 unless exists($envCounts{$ppcode});
			summary($stats{'clade_positives'}."\t$ppcode\t".scalar(@leaves)."\t".
				$envHits{$ppcode}."\t".$envCounts{$ppcode}."\t".
				scalar(keys %$primers_pairs_found_ref)."\t".join(',',(sort keys %$primers_pairs_found_ref))."\t".
				$primers_found_ref->{'amplicon_length'}."\t".$primers_found_ref->{'PRIMER_LEFT.TM'}."\t".
				$primers_found_ref->{'PRIMER_LEFT.DEGEN'}."\t".$primers_found_ref->{'PRIMER_RIGHT.TM'}."\t".
				$primers_found_ref->{'PRIMER_RIGHT.DEGEN'}.
				"\t$ppcode".'_L:'.$primers_found_ref->{'PRIMER_LEFT.SEQ'}."\t$ppcode".'_R:'.
				$primers_found_ref->{'PRIMER_RIGHT.SEQ'}."\t".
				$primers_found_ref->{'left'}."\t".
				$primers_found_ref->{'right'}."\t".
				$primers_found_ref->{'pairs'}."\t".
				join(',',@leaves)."\t".
				join(',',@extra)."\n"
				);
			print STDOUT scalar(@leaves)." sequences detected, ".(scalar(@ingroupNames)-$stats{'clade_positives_leaves_total'})." sequences to go...\n";
		}else{
			#if not, abandon ship :)
			print STDOUT "Unexpectidly, this route led nowhere :( Move on to next attempt...\n";
			print STDERR "______________________________________________________________\n".
				     "No suitable primers found, what happened ????...\n".
				     "______________________________________________________________\n";
			$stats{'clade_negatives'}++;
			$stats{'clade_negatives_leaves_total'} += scalar(@leaves);
			$stats{'clade_negatives_list'}.= ($stats{'clade_negatives_list'}?'|':'').join('|',@leaves) if scalar(@leaves)==1;
			if ($singleSeqFail >= scalar(@ingroupNames)/4){
				print STDOUT "/!\\ Over 25% target single sequences not amplifyable, abandon this run :(\n";
				print STDERR "/!\\ Over 25% target single sequences not amplifyable, abandon this run :(\n";
				last SEQ;
			}
			#summary("$stats{'clade_negatives'}\t$ppcode\t".(scalar(@leaves)>1?'no_internal_node':'no_terminal_leaf').
			#	"\t".scalar(@leaves)."\t".join(',',@leaves)."\n",'.reject');
		}
	}
	if ($ppnb > 0){
		my $ppcode=$runID."_PP$ppnb";
		if (! -f "$opt_d/$opt_t.".$ppcode."_arch.tar.gz"){
			do_cmd("tar cfz $opt_d/$opt_t.".$ppcode."_arch.tar.gz $opt_d/$opt_t.$ppcode.*");
			do_cmd("rm $opt_d/$opt_t.$ppcode.*");
		}
	}

}
############################################################################################################################
sub compute_primers_Primer3Plus{
	my $tree = shift;
	my $current_node = shift;
	my $domst_ref = shift;
	my $domen_ref = shift;
	my (%primer_pairs,%primers,$pf);
	$primers{'left'}='';$primers{'right'}='';$primers{'pairs'}='';

	#unless(-e "$tree.primer3.is.done" && $opt_c =~/m/){
		#unlink("$tree.primer3.is.done") if -e "$tree.primer3.is.done";
		foreach my $dom (sort keys %$domst_ref){
			$pf = write_param_file("$tree.cons",$dom,$domst_ref->{$dom},$domen_ref->{$dom});
			#$debug && print STDERR "Forking $primer3_exec $pf > $pf.out 2>&1\n";
			#$manager->start and next;
	      		do_cmd("$primer3_exec $pf > $pf.out 2>&1");
	      		#$manager->finish;
		}
		#$debug && print STDERR "Waiting for forked processes to finish..........\n";
		#$manager->wait_all_children;
		#do_cmd("touch $tree.primer3.is.done");
	#}
	foreach my $dom (sort keys %$domst_ref){
		parse_primer3("$tree.cons.$dom.p3.out",$current_node,$dom,\%primer_pairs,\%primers,scalar(keys %$domst_ref));
		do_cmd("rm $tree.cons.$dom.p3.out");
	}
	return (\%primer_pairs,\%primers);
}
############################################################################################################################
sub nuc_degeneracy{
	my $sequence = shift;
	my @nucs=split //, $sequence;
	my %degenerate = ('N' => 4, 'B' => 3,'D' => 3,'H' => 3,'V' => 3,'K' => 2,'Y' => 2,
			  'S' => 2,'W' => 2,'R' => 2,'M' => 2,'A' => 1,'T' => 1,'G' => 1,'C' => 1);
	my $degeneracy = 1;
	
	foreach my $nuc(@nucs){
		die "Nuc '$nuc'=$degenerate{$nuc} is unknown in sequence '$sequence'" unless exists( $degenerate{$nuc});
		$degeneracy *= $degenerate{$nuc} if $degeneracy < 10000000;
	}
	return $degeneracy;
}
############################################################################################################################
sub nuc_composition{
	my $sequence = shift;
	my $fseparator= shift;
	my $rseparator= shift;
	$sequence=uc($sequence);
	#my %composition_hash;
	my $tot_length=length($sequence);
	my $seq = $sequence;
	my $tot_Ns    = ($seq =~ s/N//g);
	my $tot_NonAmb= ($seq =~ s/[CATG]//g);
	my $tot_Amb=length($seq)+$tot_Ns;
	my $degeneracy = nuc_degeneracy($sequence);
	my $out = ($degeneracy < 10000000?"Degeneracy:$fseparator$degeneracy$rseparator":'');
	$out .= 'Length:'.$fseparator.$tot_length.' bp'.$rseparator;
	#foreach my $nuc(('A','T','G','C'), (grep { $_ ne 'A' && $_ ne 'T' && $_ ne 'G' && $_ ne 'C' } sort keys %composition_hash)){
	#	next unless exists($composition_hash{$nuc});
	#	$out .= sprintf("$nuc:$fseparator%i$fseparator%5.2f%%$rseparator", $composition_hash{$nuc}, 100.0*$composition_hash{$nuc}/$tot_length);
	#}
	$out .= sprintf("ATGC:$fseparator%i$fseparator%5.2f%%$rseparator", $tot_NonAmb, 100.0*$tot_NonAmb/$tot_length);
	$out .= sprintf("Ambg:$fseparator%i$fseparator%5.2f%%$rseparator" , $tot_Amb, 100.0*$tot_Amb/$tot_length);
	$out .= sprintf("Ns:$fseparator%i$fseparator%5.2f%%$rseparator", $tot_Ns, 100.0*$tot_Ns/$tot_length);
	$max_ambig = $tot_Amb if $tot_Amb > $max_ambig;
#	$max_degen = $degeneracy if $degeneracy > $max_degen;
	$out .= $oligo_tm_cache{$sequence} if ((defined $opt_y && $degeneracy <= $opt_y) || (! defined $opt_y && length($sequence) < 35));
	return $out;
}
############################################################################################################################
sub inflate_degen{
	my $oligo_prev = shift;
	my $oligo_next = shift;
	my $file = shift;
	my $nb = 0;
	ATGC: for (my $i=0; $i < length($oligo_next); $i++){
		my $nuci = substr($oligo_next,$i,1);
		if ($nuci !~ /^[ATGC]$/o){
			my $inflated = $amb_tbl{$nuci};
			for ( my $j=0; $j < length($inflated); $j++ ){
				my $nucj = substr($inflated,$j,1);
				$nb += inflate_degen($oligo_prev.$nucj,substr($oligo_next,$i+1),$file);
			}
			return $nb;
		}else{
			$oligo_prev .= $nuci;
		}
	}
	print $file "$oligo_prev\n";# if $i == length($oligo_next);
	$nb++;
	return $nb;
}

############################################################################################################################
sub melting_point{
	my $seq = shift;
	my $i = 0;
	my $max = 0;
	my $min = 10000;
	my $mean = 0;
	foreach my $line (split(/\n/,$seq)){
#Num	Sequence                                	Length	[Oligo]	 [Salt]	%CG	Bas	Sal	Bre	San	Sug	Con	Class	Comments
#1	GCTCAACTAGAACCCCGAGGTAA                 	23	5.0e-07	 0.050	52.2	57.1	64.6	61.9	52.7	56.6	54.6	1	OK
#2	GCTCAACTAGAACCCCGAGGTAG                 	23	5.0e-07	 0.050	56.5	58.8	66.4	61.6	52.9	57.5	55.2	1	NCF
		my @cols = split(/\s+/,$line);
		$i++;
		my $tm;
		if ($cols[13] eq 'OK' || $cols[13] eq 'NCF'){
			$tm = $cols[11];
		}elsif($cols[13] eq 'LEN'){
			$tm = $cols[9];
		}else{
			print STDERR "Failed Tm calculation for inflated $seq ($cols[13])\n";
			$tm = 0;
		}
		$min = $tm if $tm < $min;
		$max = $tm if $tm > $max;
		$mean += $tm;
	}
	if ($i) {$mean = $mean / $i;}else{$mean = 0;}
	return "Tm: mean=".sprintf("%5.2f",$mean)." min=$min max=$max";
}
############################################################################################################################
sub calculate_tm{
	foreach my $file (($opt_f, $opt_r)){
		do_cmd("$seqkit seq -s $file > $opt_d/oligos.raw");
		pre_calculate_tm("$opt_d/oligos.raw");
		open (IN , $file) or die "$!";
		while (my $line = <IN>){
			chomp($line);
			if ($line =~ /^>(\S+)/){
				my $id = $1;
				$line = <IN>;
				chomp($line);
				$debug == 2 && print STDERR "Compute Tm for $id: $line\n";
				$line = "$id\t".nuc_composition($line," ","\t")."\n";
				print STDERR $line;
				print $line;
			}
		}
		close IN;
	}
}

############################################################################################################################
sub pre_calculate_tm{
	my $cfg_file = shift;
	my %oligo_inflated;
	my $oligo_list;
	if ($opt_k ne 'degentm'){
		$oligo_list = `egrep "^PRIMER_(LEFT|RIGHT)_[0-9]+_SEQUENCE" $cfg_file > $cfg_file.oligo_list`;
	}else{
		$oligo_list = `cat $cfg_file > $cfg_file.oligo_list`;
	}
	$debug && say STDERR "Result of pre_calculate_tm oligo extraction from $cfg_file : $oligo_list";
	open (P3, "$cfg_file.oligo_list") || die "$!";
	#PRIMER_LEFT_993_SEQUENCE=GTTGATGAAACAAGGGCRGCAT
	#PRIMER_RIGHT_993_SEQUENCE=ATTTGCTGAGTAAGGGTCGCT
	open my $file,">$cfg_file.oligo_list.inflated";
	my $tot = 0;
	my @seq_order;
	OLIGO: while (my $line = <P3>){
		my $seq;
		if ($opt_k ne 'degentm'){
			($seq) = ($line =~ /=([A-Z\-]+)$/);
			$seq =~ s/\-//g;
		}else{
			$line = uc $line;
			($seq) = ($line =~ /(^[A-Z\-]+)$/);
			$debug && print STDERR "Oligo $seq from $line\n";
		}
		if(defined $seq){
			if (exists $oligo_tm_cache{$seq}){
				$debug && print STDERR "Using cached Tm for $seq = $oligo_tm_cache{$seq}°C\n";
			}else{
				my $degen = nuc_degeneracy($seq);
				if ( !defined $opt_y || (defined $opt_y && $degen <= $opt_y)  ){
					$oligo_inflated{$seq} = inflate_degen('',$seq,$file);
					$tot += $oligo_inflated{$seq};
					push @seq_order, $seq;
					$debug && print STDERR "Inflated $seq into $oligo_inflated{$seq} oligos ($tot total)\n";
				}else{
					$debug && print STDERR "Oligo $seq high degeneracy ($degen) so not calculating Tm...\n";
					$oligo_tm_cache{$seq} = "Tm: mean=0 min=0 max=0";
				}
			}
		}else{
			$debug && print STDERR "No oligo sequence so not calculating Tm: $line\n";
		}
	}
	close P3;
	close $file;
	if ($tot){
		$debug && say STDERR "sum of $tot versus Total Inflated in file = ".`wc -l $cfg_file.oligo_list.inflated`;
		do_cmd("$dnaMate $cfg_file.oligo_list.inflated > $cfg_file.oligo_list.inflated.tm");
		open TM, "$cfg_file.oligo_list.inflated.tm";
		<TM>;
		foreach my $oligo (@seq_order){
			my $lines = '';
			for (my $i=0; $i < $oligo_inflated{$oligo}; $i++){
				$lines .= <TM>;
			}
			$oligo_tm_cache{$oligo} = melting_point($lines);
			$debug && print STDERR "Computed Tm for $oligo : $oligo_inflated{$oligo} inflated oligos = $oligo_tm_cache{$oligo}\n";
		}
		close TM;
		do_cmd("rm $cfg_file.oligo_list.inflated.tm");
	}
	do_cmd("rm $cfg_file.oligo_list.inflated");
	do_cmd("rm $cfg_file.oligo_list");
}

############################################################################################################################
sub parse_primer3{
	my $cfg_file = shift;
	my $current_node = shift;
	my $domain = shift;
	my $primer_pairs_ref = shift;
	my $primers_ref = shift;
	my $nb_domains = shift;
	my $prefix = ($nb_domains>1?"$domain:":'');
	my %passed;
	my $pindex=0;
	my $nb_passable = 0;
	my $nb_unpassable = 0;
	my $ldegenmin = 1000000000;
	my $rdegenmin = 1000000000;
	my $lTmmin = 1000000000;
	my $rTmmin = 1000000000;
	my $lTmmax = 0;
	my $rTmmax = 0;
	my $codegenmin = 100000000000;
	my $codegenpair = '';
	my $ldegen; my $winindex = 0; my $wincodegen = 100000000000; my $lTm;
	#my $passed_right=0;my$passed_left=0;
	pre_calculate_tm($cfg_file);
	open (P3, $cfg_file) || die "$!";
	PAIR: while (my $line = <P3>){
		if ($line =~ /^(PRIMER_WARNING|PRIMER_LEFT_EXPLAIN|PRIMER_RIGHT_EXPLAIN|PRIMER_PAIR_EXPLAIN)=(.+)$/){
			print STDERR "$1=$2\n";
			$primers_ref->{'left'}.=(length($primers_ref->{'left'})?',':'').
				"$prefix$1" if $line =~ /^PRIMER_LEFT_EXPLAIN=.+ok (\d+)\s*\n$/;
			$primers_ref->{'right'}.=(length($primers_ref->{'right'})?',':'').
				"$prefix$1" if $line =~ /^PRIMER_RIGHT_EXPLAIN=.+ok (\d+)\s*\n$/;
			$primers_ref->{'pairs'}.=(length($primers_ref->{'pairs'})?',':'').
				"$prefix$1" if $line =~ /^PRIMER_PAIR_EXPLAIN=.+ok (\d+)\s*\n$/;
		}elsif ($line =~ /^PRIMER_PAIR_NUM_RETURNED=(.+)\n$/){
			print STDERR "PRIMER_PAIR_NUM_RETURNED=$1\n";
		#}elsif (my ($dirt,$ppindext,$tm) = ($line =~ /^PRIMER_(LEFT|RIGHT)_(\d+)_TM=(.+)\n$/)){
		#	$ppindext++;
		#	$pindex=$ppindext;
		#	$passed{"TM_$dirt".'_'.$ppindext}=$tm;
		#	$debug == 2 && print STDERR " =>Primer $dirt $ppindext Tm = $tm\n";
		}elsif (my ($ppindexs,$ps) = ($line =~ /^PRIMER_PAIR_(\d+)_PRODUCT_SIZE=(.+)\n$/)){
			$ppindexs++;$pindex=$ppindexs;
			$passed{"PRODSIZE_$ppindexs"}=$ps;
			$debug == 2 && print STDERR " =>Product $ppindexs size = $ps\n";
		}elsif ($line =~ /^PRIMER_ERROR=(.+)\n$/){
			print STDERR "+++++++++++++++++++++++++++++++++++++++++++++++++++\n".
				     "CRITICAL ERROR! PLEASE CHECK PRIMER3 CONFIG!!!\nPRIMER_ERROR=$1\n".
				     "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		}elsif( my ($dir,$ppindex,$seq) = ($line =~ /^PRIMER_(LEFT|RIGHT)_(\d+)_SEQUENCE=([A-Za-z]+)\n$/)){
			$ppindex++;$pindex=$ppindex;
			my $ncomp = nuc_composition($seq," ",", ");
			(my $degen) = ($ncomp =~ /Degeneracy: (\d+)/);
			(my $Tm) = ($ncomp =~ /Tm: mean=(\d+\.?\d*)/);#TODO maybe also allow thresholding on *min* Tm value ?
			$Tm = 0 unless $Tm;
			$ldegen = $degen if $dir eq 'LEFT';
			$lTm = $Tm if $dir eq 'LEFT';
			$passed{"TM_$dir".'_'.$pindex}=$Tm;
			if ($Tm){
				$rTmmin = $Tm if ($rTmmin > $Tm && $dir eq 'RIGHT');
				$lTmmin = $Tm if ($lTmmin > $Tm && $dir eq 'LEFT');
				$rTmmax = $Tm if ($rTmmax < $Tm && $dir eq 'RIGHT');
				$lTmmax = $Tm if ($lTmmax < $Tm && $dir eq 'LEFT');
			}
			if ($degen){
				$rdegenmin = $degen if ($rdegenmin > $degen && $dir eq 'RIGHT');
				$ldegenmin = $degen if ($ldegenmin > $degen && $dir eq 'LEFT');
			}
			my $codegen;
			if ($dir eq 'RIGHT' && defined($ldegen) && defined($degen)){
				$codegen = ($opt_D eq '*'?$ldegen * $degen:$ldegen + $degen);
				if ($codegen < $codegenmin){
					$codegenmin = $codegen;
					$codegenpair = "(LEFT DEGEN=$ldegen, RIGHT DEGEN=$degen, L.Tm=$lTm, R.Tm=$Tm, pair N°$pindex)";
				}
			}
			$passed{"$dir.DEGEN_".$pindex} = $degen;
			$passed{"$dir.SEQ_".$pindex} = $seq;
			if (defined($degen) && $degen <= $opt_y && (! defined($opt_e) || ($Tm >= $opt_e && $Tm <= $opt_E) ) ) {
				#$passed{"$dir.NCOMP"} = $ncomp;
				$debug == 2 && print STDERR " =>Primer $dir $ppindex acceptable because of low degeneracy ($degen) & decent Tm ($Tm) '$seq' !\n\t$ncomp\n";
			}else{
				$degen = 'undef' unless defined $degen;
				$debug == 2 && print STDERR " =>Primer $dir $ppindex rejected because of high degeneracy ($degen) or low Tm ($Tm) '$seq'!\n\t$ncomp\n";
			}
		}else{next PAIR;}
		if (defined($passed{'RIGHT.DEGEN_'.$pindex}) && $passed{'RIGHT.DEGEN_'.$pindex} <= $opt_y && defined($passed{'LEFT.DEGEN_'.$pindex}) && $passed{'LEFT.DEGEN_'.$pindex} <= $opt_y &&
		    exists($passed{'PRODSIZE_'.$pindex}) && 
		    ( ! defined($opt_e) || 
		    	$passed{"TM_RIGHT_".$pindex} >= $opt_e && $passed{"TM_LEFT_" .$pindex} >= $opt_e  &&
		    	$passed{"TM_RIGHT_".$pindex} <= $opt_E && $passed{"TM_LEFT_" .$pindex} <= $opt_E  &&
		    abs( $passed{"TM_RIGHT_".$pindex} - $passed{"TM_LEFT_" .$pindex} ) <= $opt_T 
		    	) ) {
			$nb_passable++;
			print STDERR " =>Primer pair N°".$pindex." passes all thresholds (DEGEN_LEFT=".$prefix.$passed{'LEFT.DEGEN_'.$pindex}.", DEGEN_RIGHT=".$prefix.$passed{'RIGHT.DEGEN_'.$pindex}.", Tm.L=".$passed{"TM_LEFT_".$pindex}."°C Tm.R=".$passed{"TM_RIGHT_".$pindex}."°C)\n" if $nb_passable <= 10;
			my $codegen = ($opt_D eq '*'?$passed{'RIGHT.DEGEN_'.$pindex} * $passed{'LEFT.DEGEN_'.$pindex}:$passed{'RIGHT.DEGEN_'.$pindex} + $passed{'LEFT.DEGEN_'.$pindex});
			if ( $wincodegen > $codegen ) {
				$winindex = $pindex;
				$wincodegen = $codegen;
			}
		}elsif( defined($passed{'RIGHT.DEGEN_'.$pindex}) && defined($passed{'LEFT.DEGEN_'.$pindex}) && exists($passed{'PRODSIZE_'.$pindex}) ){
			$nb_unpassable++;
			print STDERR " =>Primer pair N°".$pindex." FAILS thresholds (DEGEN_LEFT=".$prefix.$passed{'LEFT.DEGEN_'.$pindex}.", DEGEN_RIGHT=".$prefix.$passed{'RIGHT.DEGEN_'.$pindex}.", Tm.L=".$passed{"TM_LEFT_".$pindex}."°C Tm.R=".$passed{"TM_RIGHT_".$pindex}."°C)\n" if $nb_unpassable <= 10;
		}
	}
	if ($nb_passable && $winindex){
		$primers_ref->{"PRIMER_RIGHT.TM"}   .=(length($primers_ref->{"PRIMER_RIGHT.TM"})?',':'')   .$prefix.$passed{"TM_RIGHT_" .$winindex};
		$primers_ref->{"PRIMER_LEFT.TM"}    .=(length($primers_ref->{"PRIMER_LEFT.TM" })?',':'')   .$prefix.$passed{"TM_LEFT_"  .$winindex};
		$primers_ref->{'amplicon_length'}   .=(length($primers_ref->{'amplicon_length'})?',':'')   .$prefix.$passed{"PRODSIZE_".$winindex};
		$primers_ref->{"PRIMER_LEFT.DEGEN"} .=(length($primers_ref->{"PRIMER_LEFT.DEGEN"})?',':'') .$prefix.$passed{'LEFT.DEGEN_'.$winindex};
		$primers_ref->{"PRIMER_RIGHT.DEGEN"}.=(length($primers_ref->{"PRIMER_RIGHT.DEGEN"})?',':'').$prefix.$passed{'RIGHT.DEGEN_'.$winindex};
		$primers_ref->{"PRIMER_RIGHT.SEQ"}  .=(length($primers_ref->{"PRIMER_RIGHT.SEQ"})?',':'')  .$prefix.$passed{'RIGHT.SEQ_'.$winindex};
		$primers_ref->{"PRIMER_LEFT.SEQ"}   .=(length($primers_ref->{"PRIMER_LEFT.SEQ"})?',':'')   .$prefix.$passed{'LEFT.SEQ_'.$winindex};
		$primer_pairs_ref->{$domain}=1;#$1 if $1 > 0;
		print STDERR "  Selected primer pair: N°$winindex (DEGEN=$passed{'LEFT.DEGEN_'.$winindex} $opt_D $passed{'RIGHT.DEGEN_'.$winindex}, Tm.L=".$passed{"TM_LEFT_".$winindex}."°C Tm.R=".$passed{"TM_RIGHT_".$winindex}."°C)\n";
	}
	print STDERR "  Min LEFT degen = $ldegenmin\n  Min RIGHT degen = $rdegenmin\n";
	print STDERR "  Min combined degen = $codegenmin $codegenpair\n";
	print STDERR "  Min LEFT Tm = $lTmmin\n  Min RIGHT Tm = $rTmmin\n";
	print STDERR "  Max LEFT Tm = $lTmmax\n  Max RIGHT Tm = $rTmmax\n";
	close P3;
}

############################################################################################################################
sub write_param_file{
	my $seq_file = shift;
	my $domain = shift;
	my $start = shift;
	my $end = shift;
	open (SEQ, $seq_file) || die "$!";
	my $seq = '';
	while (my $line = <SEQ>){
		next if $line =~ /^>/;
		chomp($line);
		$seq =~ s/[^A-Za-z\-]//g;
		$seq .= $line;
	}
	die "BAD consensus start/stop: $start-$end for seq length =".length($seq) 
		if ($start > length($seq) or $start > $end or $end > length($seq) or $start < 1 or $end < 1);
	close SEQ;
	open (CFG, $opt_m) || die "$opt_m : $!";
	open (RUN_CFG, ">$seq_file.$domain.p3") || die "$!";
	while (my $line = <CFG>){
		if ($line =~ /^SEQUENCE_ID=/){
			$line = "SEQUENCE_ID=$seq_file.$domain\n";
		}elsif ($line =~ /^SEQUENCE_TEMPLATE=/){
			$line = "SEQUENCE_TEMPLATE=$seq\n";
		}elsif ($line =~ /^SEQUENCE_INCLUDED_REGION=/){
			if ($opt_s ne 'include') {
				$line = '';
			}else{
				$line = "SEQUENCE_INCLUDED_REGION=$start,".($end-$start+1)."\n";
			}
		}elsif ($line =~ /^SEQUENCE_TARGET=/){
			if ($opt_s ne 'target') {
				$line = '';
			}else{
				$line = "SEQUENCE_TARGET=$start,".($end-$start+1)."\n";
			}
		}elsif ($line =~ /^PRIMER_NUM_RETURN=/ && $opt_x){
			$line = "PRIMER_NUM_RETURN=$opt_x\n";
		}elsif ($line =~ /^SEQUENCE_EXCLUDED_REGION=/ || $line =~ /^=/){
			$line = '';
		}
		print RUN_CFG $line;
	}
	open EX,"$seq_file.mask" or die "MASK file $seq_file.mask: $!";
	my $prev = - 10;my $mstart=-10;my $len=0;my $excludes = '';
	while (my $line = <EX>){
		#51	51	1	INDEL
		(my $pos) = ($line =~ /^(\d+)/);
		if ($pos == $prev + 1){
			$len++;
		}else{
			if ($mstart >= 0){
				$excludes.=($excludes?' ':'')."$mstart,$len";
			}
			$mstart=$pos;
			$len=1;
		}
		$prev = $pos;
	}
	$excludes.=($excludes?' ':'')."$mstart,$len" if $mstart >= 0;
	close EX;
	if ($excludes) {
		print RUN_CFG "SEQUENCE_EXCLUDED_REGION=$excludes\n";
		$debug && print STDERR "SEQUENCE_EXCLUDED_REGION=$excludes\n";
	}
	print RUN_CFG "=\n";
	close RUN_CFG;
	close CFG;
	return "$seq_file.$domain.p3";
}

############################################################################################################################
sub stats{
	my @leaves = @ingroupNames;
	my $t = 
	"\tTotal nb of sequences\t".scalar(@leaves)."\n".
	"\tTotal nb of amplifiable clusters\t".($stats{'clade_positives'})."\n".
	"\tTotal nb of seqs in ampl. clades\t".($stats{'clade_positives_leaves_total'})."\n".
	"\t% of seqs in ampl. clades\t". ( scalar(@leaves)?nearest(.1, ($stats{'clade_positives_leaves_total'} / scalar(@leaves)*100) )."\%\n" : "0\%\n" ).
	"\tAvg nb sequences in ampl. clades\t".( $stats{'clade_positives'}?nearest(.01, ($stats{'clade_positives_leaves_total'} / $stats{'clade_positives'}) ) : '0' )."\n".
	"\tNb leaves in positives list\t".scalar(split(/\|/,$stats{'clade_positives_list'}) )."\n".
	"\tBottom-up clade compression ratio\t".($stats{'clade_positives'}?nearest(.1,scalar(@leaves) / $stats{'clade_positives'})."\n":"0\%\n").
	"\tMaximum number of ambiguous nucleotides per primer:\t$max_ambig\n";
	#"\tAverage nucleotide composition of the positive primer set:\n".
	#nuc_composition($compiled_primers,"\t","\n");
	#my $tt = 
	#"QC: List of positive sequences\t\t\t".$stats{'clade_positives_list'}."\n".
	#"QC: List of negative sequences\t\t\t".$stats{'clade_negatives_list'}."\n";
	summary($t);my $cnt=0;
	foreach my $leaf (sort @leaves){
		$cnt++;
		if ($stats{'clade_positives_list'} !~ /$leaf/ && $stats{'clade_negatives_list'} !~ /$leaf/) {
			print STDERR "$leaf is neither pos nor neg!!!!!!!!!!!!!!\n";
		}elsif($stats{'clade_positives_list'} =~ /$leaf/ && $stats{'clade_negatives_list'} =~ /$leaf/){
			print STDERR "$leaf is both pos and neg!!!!!!!!!!!!\n";
		}elsif($stats{'clade_positives_list'} =~ /$leaf/){
			#print STDERR "$leaf is +\n";
		}elsif($stats{'clade_negatives_list'} =~ /$leaf/){
			print STDERR "$leaf is - !!!!!!!!!!!\n";
		}else{
			print STDERR "$leaf is gobledeegook!!!!!!!!!!!!!!\n";
		}
	}
	print STDERR "##################################################\n$t\n";
	print STDOUT "##################################################\n$t\n";
}
sub wrap_up{
	open OUT, ">$opt_d/$opt_t.css";
	print OUT $css;
	close OUT;
	open OUT, ">$opt_d/$opt_t.ornaments";
	print OUT $orn;
	close OUT;
	if (-f "$opt_d/$opt_t.primers.R.fasta" && -f "$opt_d/$opt_t.primers.F.fasta"){
		do_cmd("$nw_utils/nw_display -sr -w 1200 -b 'opacity:0' -i 'opacity:0' -l 'font-size:4px;font-family:sans;' -v 6 -c $opt_d/$opt_t.css -o $opt_d/$opt_t.ornaments $opt_d/$opt_t.rerooted.renamed > $opt_d/$opt_t.tree.svg");
		do_cmd_cached("revseq -sequence $opt_d/$opt_t.primers.R.fasta -outseq $opt_d/$opt_t.primers.R.rc.fasta","$opt_d/$opt_t.primers.R.rc.fasta",'m');
		do_cmd_cached("cat $opt_d/$opt_t.primers.R.rc.fasta $opt_d/$opt_t.primers.F.fasta > $opt_d/$opt_t.primers.FRrc.fasta","$opt_d/$opt_t.primers.FRrc.fasta",'m');
		if (! -f "$opt_d/$opt_t.ingroup.fna.aln"){
			do_cmd("cp $opt_n $opt_d/$opt_t.ingroup.fna.aln");
		}
		do_cmd_cached("$seqret $opt_d/$opt_t.ingroup.fna.aln -outseq $opt_d/$opt_t.ingroup.fna.fasta.aln","$opt_d/$opt_t.ingroup.fna.fasta.aln",'m');
		do_cmd_cached("$mafft --mapout --addfragments $opt_d/$opt_t.primers.FRrc.fasta $opt_d/$opt_t.ingroup.fna.fasta.aln > $opt_d/$opt_t.primers.aln","$opt_d/$opt_t.primers.aln",'m');
	}
	print STDERR "/////////////////////////\nAll said & done.\n/////////////////////////\n\n";
	print STDOUT "///////////////////////// we're done ! /////////////////////////\n\n";
	
}
############################################################################################################################
sub add_css{
	my $type = shift;
	my $node = shift;
	#$col+=4;$col=0 if $col > 15;
	$css .= "\"stroke-width:2; stroke:#".$colors[$col]."\" $type $node\n";
}
############################################################################################################################
sub add_ornaments_REFS{
	my @nodes = @_;
	foreach my $node (@nodes){
		$orn .= "\"<circle style='fill:red;stroke:none' r='4'/>\" Individual $node\n" if $node =~ /^REF/;
	}
}
sub add_ornaments_pp{#add_ornaments_pp(\@pps,$cid);
	my $pps_ref = shift;
	my $node = shift;
	#$orn .= "\"<circle style='fill:blue;stroke:none' r='4'/>\" Individual $node\n";
	#$orn .= "\"<rect style='fill:#$colors[$col];stroke:none' width='5' height='4'/>\" Individual $node\n";
	return unless scalar @$pps_ref;
	$orn .= "\"";
	my $i = 70;
	foreach my $pid (sort @$pps_ref){#$debug && print STDERR "$node:$pid\n";
		if (! defined $pp2color{$pid}){
			$col+=1;$col=0 if $col >= scalar(@colors);
			$pp2color{$pid}=$col;
		}
		$orn .= "<rect style='fill:#$colors[$pp2color{$pid}];stroke:none' x='+$i' y='-3' width='5' height='3'/>";
		$i+= 5;
	}
	$orn .= "\" Individual $node\n"; 
}
############################################################################################################################
sub do_cmd{
	my $cmd = shift;
	$debug && print STDERR "\nExecute: $cmd\n";
	system ($cmd)  == 0 || die "COMMAND '$cmd' FAILED: $!";
}
############################################################################################################################
sub do_cmd_cached{
	my $cmd = shift;
	my $file = shift;
	my $code = shift;
	if (cached($file,$code)){
		$debug && print STDERR "\nUsing cached: $file\n";
	}else{
		do_cmd($cmd);
		#do_cmd("touch $file.$code.cached");
	}
}
############################################################################################################################
sub cached{
	my $file = shift;
	my $code = shift;
	if (-e $file && -e "$file" && $opt_c && $opt_c =~/$code/){
		return 1;
	}else{
		unlink "$file.$code.cached" if -e "$file.$code.cached";
		unlink "$file" if -e "$file";
		return 0;
	}
}

############################################################################################################################
sub sub_fasta{
	my $tree = shift;
	my $fasta_file = shift;
	my $out_file = shift;
	my $id_list_ref = shift;
	unless(cached($out_file,'f')){
		my @ids = @$id_list_ref;#get_leaves($tree);
		open IDS, ">$tree.ids" || die "$!";
		print IDS join("\n",@ids);
		close IDS;
		#do_cmd_cached("$fasta_get -f $fasta_file -l $tree.ids > $out_file","$out_file",'f');
		extract_fasta("$tree.ids","$fasta_file","$out_file");
	}else{
		$debug && print STDERR "\nUsing cached $out_file\n";
	}
	return $out_file;
}

############################################################################################################################
sub get_leaves{
	my $tree = shift;
	my $cmd = "$nw_utils/nw_labels -It $tree";
	$debug && print STDERR "\n$cmd\n";
	my $o = `$cmd`;
	#print STDERR $o;
	chomp($o);
	return (split /\t/, $o);
}

############################################################################################################################
sub summary{
	my $text = shift;
	#the summary.csv will report outcome of primer search for each subclade
	if (!$sm && (-e "$opt_d/$opt_t.summary.csv" )){
		unlink "$opt_d/$opt_t.summary.csv" if -e "$opt_d/$opt_t.summary.csv";
	}
	open (SUMMARY, ">>$opt_d/$opt_t.summary.csv") || die "$!";
	print SUMMARY "NB\tCLUSTER\tSEQS_IN_CLUSTER\tNB_SEQS_HITS\tSUM_SEQ_ABUNDANCE\tNB_DOMAINS_WITH_PCR\tDOMAINS_WITH_PCR_LIST\tAMPLICON_SIZE\tPRIMER_LEFT_TM\tPRIMER_LEFT_DEGENERACY\tPRIMER_RIGHT_TM\tPRIMER_RIGHT_DEGENERACY\tPRIMER_LEFT\tPRIMER_RIGHT\tTOTAL_NB_LEFT_PRIMERS\tTOTAL_NB_RIGHT_PRIMERS\tTOTAL_NB_PRIMER_PAIRS\tSEQS_LIST\tOUTCLUSTER_IMPORTED_SEQS_LIST\n" if ( !$sm);
	print SUMMARY $text;
	$sm++;
	close SUMMARY;
}

############################################################################################################################
sub primers{
	my $text = shift;
	my $dir = shift;
	if (!$pm && (-e "$opt_d/$opt_t.primers.F.fasta" || -e "$opt_d/$opt_t.primers.R.fasta")){
		unlink "$opt_d/$opt_t.primers.F.fasta";
		unlink "$opt_d/$opt_t.primers.R.fasta";
	}
	#the primers.csv will list all found primers
	open (PRIMERS, ">>$opt_d/$opt_t.primers.$dir.fasta") || die "$!";
	#print PRIMERS "SOME_HEADER\n" unless $pm;
	print PRIMERS $text;
	$pm++;
	close PRIMERS;
}
############################################################################################################################
sub consensus{
	my $aln_file = shift;
	my $threshold = shift;
	my $output_file = shift;
	my $mask_file = shift;
	my $ignore_gaps = $opt_w;

	my $min_repr_threshold=(100.0-$threshold)/100.0;
	$debug && print STDERR "Min representation threshold: $min_repr_threshold\n";
	$debug && print STDERR "Ignoring gaps: $ignore_gaps\n";

	my ($seqname, $path)=fileparse($aln_file);
	$seqname=~s/\.aln//;

	###############################################################################

	$debug && print STDERR "Reading ALN file...\n";

	open(ALN_FH, "<$aln_file") || die "Could not open $aln_file\n";

	my %seq_hash;
	my $clustal_file=0;
	while(<ALN_FH>){
		chomp;
		if($_=~/CLUSTAL/ || $_=~/MUSCLE/){
			$clustal_file=1;
			next;	
		}elsif($_=~/^\s+[\*\:\.]/){
			next;
		}elsif($_=~/^\s*$/){
			next;
		}else{
			my ($id, $sequence)=split /\s+/, $_;
			if($id ne ""){
				$seq_hash{$id}.=uc($sequence);
			}
		}
	}
	close(ALN_FH);

	if($clustal_file!=1){
		print STDERR "Are you sure this is clustal file?\n";
	}

	###############################################################################

	$debug && print STDERR "Computing nucleotide profiles...\n";
	my @nuc_at_pos;
	my $length=-1;
	SAL: foreach my $id(keys %seq_hash){
		my @nucs=split //, $seq_hash{$id};
	
		# Assume length of consensus is the length of the first sequence, but die if not all
		#  lengths are the same length.
		my $curlength=$#nucs+1;
		if($length==-1){
			$length=$curlength;
		}elsif($length != $curlength){
			print STDERR "Error, for some reason the length of alignment for $id is not the same as all the other sequences ($length != $curlength)\n";
			delete($seq_hash{$id});
			next SAL;
		}
		# Compute the profile at each position
		#	Each position is represented in an array
		#	Each nucleotide code is summed up in a hash
		for(my $i=0; $i<$length; $i++){
			${$nuc_at_pos[$i]}{$nucs[$i]}++;
		}
	}

	###############################################################################

	$debug && print STDERR "Normalizing profiles and eliminating nucs that don't pass threshold...\n";
	open(FH, ">$mask_file") || die "Could not open $mask_file\n";
	for(my $i=0; $i<$length; $i++){

		# Sum up the total number of nucleotides at each position
		my $tot=0;
		my $gap=0;
		foreach my $nuc (keys %{$nuc_at_pos[$i]}){
			$tot+=${$nuc_at_pos[$i]}{$nuc} ;
			$gap++ if ($nuc eq '-');#if there is a gap at pos $i
		}
		if ($gap && !$ignore_gaps){
			foreach my $nuc (keys %{$nuc_at_pos[$i]}){
				delete ${$nuc_at_pos[$i]}{$nuc};
			}
			${$nuc_at_pos[$i]}{'N'}=1;
			print FH "$i\t$i\t1\tINDEL\n";
		}else{
			# Eliminated the variation at each position if it doesn't exceed the min_repr_threshold
			foreach my $nuc (keys %{$nuc_at_pos[$i]}){
				${$nuc_at_pos[$i]}{$nuc}/=$tot;
				if( $min_repr_threshold > ${$nuc_at_pos[$i]}{$nuc} ){
					delete ${$nuc_at_pos[$i]}{$nuc};
				}
			}
		}
	}
	close FH;
	###############################################################################

	$debug && print STDERR "Building consensus sequence...\n";

	my $consensus="";
	for(my $i=0; $i<$length; $i++){
		$consensus.=get_amb($nuc_at_pos[$i]);
	}

	###############################################################################

	$debug && print STDERR "Writing Consensus: $output_file\n";

	output_fasta($output_file, $consensus,$seqname);

}
###############################################################################

sub get_amb{
	my $nuc_hash_ref=shift;
	my $code='N';

	# The order of these if statements are important.

	if(has($nuc_hash_ref, "A")){
		$code="A";
	}
	if(has($nuc_hash_ref, "C")){
		$code="C";
	}
	if(has($nuc_hash_ref, "G")){
		$code="G";
	}
	if(has($nuc_hash_ref, "T")){
		$code="T";
	}
	if(has($nuc_hash_ref, "AC")){
		$code="M";
	}
	if(has($nuc_hash_ref, "AG")){
		$code="R";
	}
	if(has($nuc_hash_ref, "AT")){
		$code="W";
	}
	if(has($nuc_hash_ref, "CG")){
		$code="S";
	}
	if(has($nuc_hash_ref, "CT")){
		$code="Y";
	}
	if(has($nuc_hash_ref, "GT")){
		$code="K";
	}
	if(has($nuc_hash_ref, "ACG")){
		$code="V";
	}
	if(has($nuc_hash_ref, "ACT")){
		$code="H";
	}
	if(has($nuc_hash_ref, "AGT")){
		$code="D";
	}
	if(has($nuc_hash_ref, "CGT")){
		$code="B";
	}
	if(has($nuc_hash_ref, "GATC")){
		$code="N";
	}
	if(has($nuc_hash_ref, "N")){
		$code="N";
	}
	return($code);
}

###############################################################################

sub has{
	# See if all the members in the string, have a representive in the hash

	my $hash_ref=shift;
	my %hash=%{$hash_ref};	# Hash of nucleotide counts
	my $str=shift;		# String containing representatives of ambiguity code

	# Test all the nucs in the hash to see if they exist in the hash
	my @nucs=split //, $str;
	my $num_hits=0;
	foreach my $nuc(@nucs){
		if(defined($hash{$nuc})){
			$num_hits++;
		}	
	}

	# Only return true if all the members in the string are represented in the hash
	if($num_hits==($#nucs+1)){
		return(1);
	}else{
		return(0);
	}
}

###############################################################################

sub output_fasta{
	my $filename=shift;
	my $sequence=shift;
	my $seqname=shift;

	open(FH, ">$filename") || die "Could not open $filename\n";

	print FH ">$seqname\n";
	my $length=length($sequence);
	my $width=60;
	my $pos=0;
	do{
		my $out_width=($width>$length)?$length:$width;
		print FH substr($sequence, $pos, $width) . "\n";
		$pos+=$width;
		$length-=$width;
	}while($length>0);

	close(FH);
}

sub evaluate_primers{
	unless (-f "$opt_d/$opt_t.primers.F.fasta" && -f "$opt_d/$opt_t.primers.R.fasta"){
		return;
	}
	my $targets = $opt_z;
	my $run = "$opt_d/primer_evaluation";
	#my $blastx = $opt_e;

	my (%targetids,%pr,%matef,%mater,
	    %seen,%single,%double,$zeid,$i,%percentId,%nb100,%alLength,
	    %subjects,%primers,%ctgMinLeft,%ctgMaxRight,%ctgExtSubjStart,
	    %ctgExtSubjEnd,%isRef,%ctgExtQueryStart,%ctgExtQueryEnd,
	    %truncatedLeft,%truncatedRight,%primersghit,%primerdbhit,
	    %primerdbcnts,%primersgcnts,%primerallcnts,%ampliconSize,
	    %ampliconMaxSize,%ampliconMinSize,%ppSeqHits);
	say STDERR '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@';
	say STDERR '@             Evaluate primers                           @';
	say STDERR '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@';
	do_cmd("$seqkit seq -ni $targets > $run.$targets.ids") if ! -f "$run.$targets.ids";
	#####################################################
	open(IN, "$run.$targets.ids") or die "$!";
	$i=0;
	while (my $line = <IN>){
		chomp($line);
		$targetids{$line}=1;
		$i++;
	}
	close IN;
	print STDERR "Nb sequences in target file: $i\n";
	my $blast_res = '';
	if($opt_h){
		#####################################################
		# if in prealigned mode, no blastx file is provided, so
		# prepare blastx of target nuc seqs versus the -h region
		# first extract the -h region of the nuc aln:
		my ($rst,$ren) = ($opt_h =~ /^(\d+)\-(\d+)/);
		$ren++;
		$rst = int($rst / 3);
		$rst = 1 unless $rst;
		$ren = int($ren / 3);
		my $outfile = "$opt_n.subregions.faa";
		$outfile =~ s/^.+\///;
		if (! -f "$opt_d/$outfile.blastx"){
			do_cmd("$t_coffee -other_pg seq_reformat -in $opt_p -action +extract_block cons $rst $ren -output fasta_seq > $opt_d/$outfile");
			do_cmd("makeblastdb -in $opt_d/$outfile -parse_seqids -dbtype prot -blastdb_version 5 >/dev/null");
			#makeblastdb -in /shared/bank/pdb/current/derived_data/pdb_seqres.txt -parse_seqids -title "PDB" -dbtype prot -out pdb_blast -blastdb_version 5
			do_cmd("$blastx -max_target_seqs 10 -best_hit_overhang 0.1 -max_intron_length 500 -outfmt 6 -evalue 1e-10 -query $targets -db $opt_d/$outfile  > $opt_d/$outfile.blastx");
		}
		$blast_res = "$opt_d/$outfile.blastx";
	}else{$blast_res = $opt_b}
	open (IN,$blast_res) or die "$!";
	$i=0;
#scaffold174051_1_80	scaffold174051_1_80	100.00	164	0	0	3856	4347	1	164	2e-107	 	333
#alpha_NC_009657	alpha_NC_009657		100.000	295	0	0	1	885	1	295	0.0		618
#gamma_FN430415		alpha_NC_038861		54.045	309	122	7	1	909	1	295	3.78e-108	328
#gamma_FN430415		alpha_HQ392470		53.571	308	125	6	1	909	1	295	7.74e-108	327

	while (my $line = <IN>){
		chomp($line);
		my ($queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore) = split (/\t/, $line);
		#$debug && print STDERR "$queryId \t";
		if ($queryId eq $subjectId){
			$isRef{$queryId}=1;
			$debug ==2 && print STDERR "$queryId \t is a ref!\n"
		}#else{$debug && print STDERR "\n";}
		#$debug && print STDERR "More than one hit for $queryId, was:$percentId{$queryId}% now $percIdentity - target was $subjects{$queryId}\n" if exists($percentId{$queryId});
		if ( exists($percentId{$queryId}) && $percentId{$queryId} > $percIdentity && $queryId ne $subjectId) {
			next
		}
		$percentId{$queryId}=$percIdentity;
		$nb100{$queryId}++ if $percIdentity == 100;
		$alLength{$queryId}=$alnLength;
		$subjects{$queryId}=$subjectId;
		$ctgExtSubjStart{$queryId}=$subjectStart;
		$ctgExtSubjEnd{$queryId}=$subjectEnd;
		$ctgExtQueryStart{$queryId}=$queryStart;
		$ctgExtQueryEnd{$queryId}=$queryEnd;
		$i++;
	}
	close IN;
	print STDERR "Nb lines in blastx file: $i - nb contigs with a blast hit: ".scalar(keys %percentId)." - nb 100% id hits: ".scalar(keys %nb100)."\n";
	#####################################################
	open IN , "$opt_d/$opt_t.primers.F.fasta";
	while (my $line = <IN>){
		chomp($line);
		if ($line =~ /^>(\S+)/){
			my $id = $1;
			$line = <IN>;
			chomp($line);
			$pr{$id}=$line;
		}
	}
	close IN;
	open IN , "$opt_d/$opt_t.primers.R.fasta";
	while (my $line = <IN>){
		chomp($line);
		if ($line =~ /^>(\S+)/){
			my $id = $1;
			$line = <IN>;
			chomp($line);
			$pr{$id}=$line;
		}
	}
	close IN;
	#####################################################
	#unless (-f "$run.dreg.eval.out"){
		my $out = '';
		foreach my $id (keys %pr){
			my $prc = `echo \"$pr{$id}\" | revseq -filter`;
			($prc) = ($prc =~ /:\n(\S+)/);
			my $primer_f = pattern($pr{$id});
			$prc = pattern($prc);
			$out .= ">$id\n".   `$dreg  -sequence  $targets -pattern \"$primer_f\" -stdout 2>/dev/null`;
			$out .= ">$id rc\n".`$dreg  -sequence  $targets -pattern \"$prc\" -stdout 2>/dev/null`;
		}
		open OUT, ">$run.dreg.eval.out";
		print OUT $out ;
		close OUT;
	#}
	#####################################################
	open IN,"$run.dreg.eval.out" or die "Arggghhhhh $!";
	while (my $line = <IN>){
		chomp($line);
		$zeid = $1 if $line =~ /^>(\S+)/;
		next if $line =~ /^SeqName/ || $line =~ /^>(\S+)/;
		$debug == 2 && print STDERR "=>primer $zeid\n";
		if ($line =~ /^(\S+)\s+\d+/){
	#TARA_168_SRF_lt-0d22_G_scaffold53300_1_gene140902	2416	2434	0.000	+	regex:GG[CT]GACCG[CG]GT[CG]C[AC][AG]T[AT]T[ACG]
			my @cols = split (/\t/, $line);
			my $cid = $cols[0];
			$debug == 2 && print STDERR "-$zeid|$cols[0]\n";
			my $pid = $zeid;
			$pid =~ s/_(L|R)$//;
			if ($zeid =~ /_L$/){
				$matef{$pid}{$cid}=$cols[1];
				#if (exists($isRef{$cid})) {
					$ctgMinLeft{$cid}{$pid}=9999999999 if !exists($ctgMinLeft{$cid}{$pid});
					$ctgMinLeft{$cid}{$pid}=$cols[1] if $cols[1] < $ctgMinLeft{$cid}{$pid};
					$debug == 2 && print STDERR "ctgMinLeft{$cid}{$pid}=$ctgMinLeft{$cid}{$pid}\n";
				#}
			}
			if ($zeid =~ /_R$/){
				$mater{$pid}{$cid}=$cols[1];
				#if (exists($isRef{$cid})) {
					$ctgMaxRight{$cid}{$pid}=0 if !exists($ctgMaxRight{$cid}{$pid});
					$ctgMaxRight{$cid}{$pid}=$cols[2] if $cols[2] > $ctgMaxRight{$cid}{$pid};
					$debug == 2 && print STDERR "ctgMaxRight{$cid}{$pid}=$ctgMaxRight{$cid}{$pid}\n";
				#}
			}
			$seen{$cid}++;
			$primers{$pid}++;
		}
	}
	close IN;
	foreach my $cid (sort keys %isRef){
		foreach my $pid (sort keys %primers){
			$debug == 2 && print STDERR "$cid: -----------No left but right for $pid\n" if (!exists($ctgMinLeft{$cid}{$pid}) &&  exists($ctgMaxRight{$subjects{$cid}}{$pid}));
			$debug == 2 && print STDERR "$cid: left and right for $pid\n" if (exists($ctgMinLeft{$cid}{$pid}) &&  exists($ctgMaxRight{$subjects{$cid}}{$pid}));
			if (!exists($ctgMinLeft{$cid}{$pid}) ){
				$ctgMinLeft{$cid}{$pid} = $ctgExtQueryStart{$cid};
			}
			if (!exists($ctgMaxRight{$cid}{$pid}) ){
				$ctgMaxRight{$cid}{$pid} = $ctgExtQueryEnd{$cid};
			}
		}
	}
	#####################################################
	foreach my $cid (sort keys %targetids){
		$debug == 2 && print STDERR "Contig $cid \t";
		$debug == 2 && print STDERR "\n" if (!exists($isRef{$cid}));
		$debug == 2 && print STDERR " ===>  is a ref ctg\n" if (exists($isRef{$cid}));
		$ampliconSize{$cid} = "" unless exists($ampliconSize{$cid});
		foreach my $pid (sort keys %primers){
			$ampliconMaxSize{$pid} = 0 unless exists($ampliconMaxSize{$pid});
			$ampliconMinSize{$pid} = 1000000000 unless exists($ampliconMinSize{$pid});
			$primerallcnts{$pid} = 0 unless exists($primerallcnts{$pid});
			if (exists($matef{$pid}{$cid}) && exists($mater{$pid}{$cid})){
				$ppSeqHits{$pid} = [] if ! defined $ppSeqHits{$pid};
				push(@{$ppSeqHits{$pid}},$cid);
				$double{$cid}++;
				$primerdbhit{$pid}++;
				$primerdbcnts{$pid} += $abundances{$cid} if exists($abundances{$cid});
				$primerallcnts{$pid} += $abundances{$cid} if exists($abundances{$cid});
				$primerallcnts{$pid}++ if ($opt_k eq 'evaluate' && !exists($abundances{$cid}));
				$debug == 2 && print STDERR "\t\t$pid is double pos F:$matef{$pid}{$cid} R:$mater{$pid}{$cid}!\n";
				my $length = 0;
				if ($ctgMaxRight{$cid}{$pid} && $ctgMinLeft{$cid}{$pid}) {
					$length = $ctgMaxRight{$cid}{$pid} - $ctgMinLeft{$cid}{$pid} + 1 ;
				}
				$ampliconSize{$cid} .= ($ampliconSize{$cid}?'|':'')."$pid:".($length);
				$ampliconMaxSize{$pid} = $length if $length > $ampliconMaxSize{$pid};
				$ampliconMinSize{$pid} = $length if $length < $ampliconMinSize{$pid};
			}elsif(exists($matef{$pid}{$cid}) || exists($mater{$pid}{$cid})){
				$single{$cid}++;
				$primersghit{$pid}++;
				$primersgcnts{$pid} += $abundances{$cid} if exists($abundances{$cid});
				$primerallcnts{$pid} += $abundances{$cid} if exists($abundances{$cid});
				$primerallcnts{$pid}++ if ($opt_k eq 'evaluate' && !exists($abundances{$cid}));
				$debug == 2 && print STDERR "\t$pid is single pos F:".( exists($matef{$pid}{$cid})?$matef{$pid}{$cid}:'' ).
					" R:".( exists($mater{$pid}{$cid})?$mater{$pid}{$cid}:'' ) ."!\n";
			}
			#if (!exists($isRef{$cid}) ) {
			if ( exists($subjects{$cid}) && exists($ctgExtSubjEnd{$cid}) && exists($ctgExtSubjEnd{$subjects{$cid}}) && exists($ctgMaxRight{$subjects{$cid}}{$pid}) && exists($ctgMinLeft{$subjects{$cid}}{$pid}) && exists($ctgExtQueryStart{$subjects{$cid}}) ){
				if ($ctgExtSubjEnd{$cid} < $ctgExtSubjEnd{$subjects{$cid}}-
					($ctgExtQueryEnd{$subjects{$cid}}-$ctgMaxRight{$subjects{$cid}}{$pid})/3 ){#+85?????
						$truncatedRight{$cid}=1;
					$debug == 2  && warn "truncatedRight IMPOSSIBLE !!!" if exists($isRef{$cid});
					$debug == 2 && warn "no ctgMaxRight{$subjects{$cid}}{$pid}" if !exists($ctgMaxRight{$subjects{$cid}}{$pid});
					$debug == 2 && print STDERR "\t*********truncatedRight:$ctgExtSubjEnd{$cid} < $ctgExtSubjEnd{$subjects{$cid}}-($ctgExtQueryEnd{$subjects{$cid}}-$ctgMaxRight{$subjects{$cid}}{$pid})/3=".($ctgExtSubjEnd{$subjects{$cid}}-($ctgExtQueryEnd{$subjects{$cid}}-$ctgMaxRight{$subjects{$cid}}{$pid})/3+1)."\n";
				}
				if ($ctgExtSubjStart{$cid} > $ctgExtSubjStart{$subjects{$cid}}+
					($ctgMinLeft{$subjects{$cid}}{$pid} - $ctgExtQueryStart{$subjects{$cid}})/3 ){#+75???
						$truncatedLeft{$cid}=1;
					$debug == 2 && print STDERR "\t*********truncatedLEFT: $ctgExtSubjStart{$cid} > $ctgExtSubjStart{$subjects{$cid}}+($ctgMinLeft{$subjects{$cid}}{$pid} - $ctgExtQueryStart{$subjects{$cid}})/3=".($ctgExtSubjStart{$subjects{$cid}}+($ctgMinLeft{$subjects{$cid}}{$pid} - $ctgExtQueryStart{$subjects{$cid}})/3)."\n";
					$debug == 2  && warn "truncatedLeft IMPOSSIBLE !!!" if exists($isRef{$cid});
					$debug == 2 && warn "no ctgMinLeft{$subjects{$cid}}{$pid} !!!" if !exists($ctgMinLeft{$subjects{$cid}}{$pid});
				}
			}
			#}
		}
		delete($single{$cid}) if $double{$cid};
		#last;
	}
	print STDERR "Nb contigs: ".scalar(keys %targetids)." - nb double ctg hits: ".scalar(keys %double).' ('.nearest(.1,scalar(keys %double)/scalar(keys %targetids)*100).'%) - single: '.scalar(keys %single).' ('.nearest(.1,scalar(keys %single)/scalar(keys %targetids)*100)."%) - total hit: ".(scalar(keys %single)+scalar(keys %double))." (".nearest(.1,(scalar(keys %single)+scalar(keys %double))/scalar(keys %targetids)*100)."%)\n".
	"Truncated Left: ".scalar(keys %truncatedLeft)." - Truncated Right: ".scalar(keys %truncatedRight)."\n";
	############################################################################################################################
	my $quant_pos=0;
	my $quant_neg=0;
	open PP,">$run.dreg.eval.contig.stats.tsv" or die "$!";
	print PP "ID                     \tNbPrimerPairPos\tAmpliconSizes\tDoublePrimers\tSinglePrimers\tPercentId\tAlnLength\tTruncatedLeft\tTruncatedRight\tReadCounts\n";
	foreach my $ctg (keys %targetids){
		$double{$ctg} = 0 unless exists($double{$ctg});
		$single{$ctg} = 0 unless exists($single{$ctg});
		$truncatedLeft{$ctg} = 0 unless exists($truncatedLeft{$ctg});
		$truncatedRight{$ctg} = 0 unless exists($truncatedRight{$ctg});
		$abundances{$ctg} = '' unless exists($abundances{$ctg});
		$percentId{$ctg} = '' unless exists($percentId{$ctg});
		$alLength{$ctg} = '' unless exists($alLength{$ctg});
		$debug && print STDERR "$ctg is double primed but truncated $truncatedLeft{$ctg} $truncatedRight{$ctg}\n" if ($double{$ctg} && ($truncatedLeft{$ctg} + $truncatedRight{$ctg})>0 );
		my $pppos = scalar(split /\|/,$ampliconSize{$ctg});
		print PP "$ctg\t$pppos\t$ampliconSize{$ctg}\t$double{$ctg}\t$single{$ctg}\t$percentId{$ctg}\t$alLength{$ctg}\t$truncatedLeft{$ctg}\t$truncatedRight{$ctg}\t$abundances{$ctg}\n";
		if ($abundances{$ctg}){
			$quant_pos += $abundances{$ctg} if ($double{$ctg} || ($single{$ctg} && ($truncatedLeft{$ctg} || $truncatedRight{$ctg}) ) );
			$quant_neg += $abundances{$ctg} if (!$double{$ctg} || ($single{$ctg} && !$truncatedLeft{$ctg} && !$truncatedRight{$ctg} )  );
		}
	}
	close PP;
	print STDERR "Abundances POS:$quant_pos NEG:$quant_neg\n";
	my $sumFile = 0;
	my (%sum,%sic);
	open PP,">$run.dreg.eval.primer.stats.tsv" or die "$!";
	print PP "PrimerPairID\tPrimerLID\tPrimerLSeq\tPrimerLComposition\tPrimerRID\tPrimerRSeq\tPrimerRComposition\tAmpliconMinSize\tAmpliconMaxSize\tNbDoubleHits\tNbUniqueHits\tNbSingleHits\tAbundanceDoubleHits\tAbundanceSingleHits\n";
	if (-f "$opt_d/$opt_t.summary.csv"){
		open IN, "$opt_d/$opt_t.summary.csv";
		while (<IN>){
			my($id,$cluster,$seqsincl,$nbtarahits,$sumtarareads,@rest) = split /\t/;
			$sum{$cluster}=join ("\t",@rest);
			$sic{$cluster}=$seqsincl;
		}
		close IN;
		$sumFile = 1;
	}

	if ($sumFile) {
		open (SM,">$run.summary.csv") or die "$!";
		print SM "CLUSTER\tSEQS_IN_CLUSTER\tNB_SEQ_AMPLIFIED\tNB_SEQ_UNIQUE\tSUM_SEQ_ABUNDANCE\t".$sum{'CLUSTER'};
	}
	%primer_listing = %primers if $opt_k eq 'evaluate';
	foreach my $pid (sort {$primerallcnts{$b} <=> $primerallcnts{$a}} keys %primer_listing){
		$primersghit{$pid} = 0 unless exists($primersghit{$pid});
		$primersgcnts{$pid} = 0 unless exists($primersgcnts{$pid});
		$primerdbhit{$pid} = 0 unless exists($primerdbhit{$pid});
		$primerdbcnts{$pid} = 0 unless exists($primerdbcnts{$pid});
		my $uniqueHits = 0;
		if (defined($pid) && exists( $ppSeqHits{$pid} )){
			$uniqueHits = scalar @{$ppSeqHits{$pid}};
			foreach my $cid (@{$ppSeqHits{$pid}}){
				my $isHit = 0;
				PPS: foreach my $pid2 (keys %primer_listing){
					next PPS if $pid2 eq $pid;
					if (grep /^\Q$cid\E$/, @{$ppSeqHits{$pid2}}){
						$isHit = 1;
						last PPS;
					}
				}
				$uniqueHits-- if $isHit;
			}
		}
		my $pAnnot = $pid."_L\t".$pr{$pid.'_L'}."\t". nuc_composition($pr{$pid.'_L'}," ","|")."\t".$pid."_R\t".$pr{$pid.'_R'}."\t". nuc_composition($pr{$pid.'_R'}," ","|");
		print PP "$pid\t$pAnnot\t$ampliconMinSize{$pid}\t$ampliconMaxSize{$pid}\t$primerdbhit{$pid}\t$uniqueHits\t$primersghit{$pid}\t$primerdbcnts{$pid}\t$primersgcnts{$pid}\n";
		$sumFile && print SM "$pid\t$sic{$pid}\t$primerdbhit{$pid}\t$uniqueHits\t$primerdbcnts{$pid}\t".$sum{$pid};
	}
	$sumFile && close SM;
	close PP;
	foreach my $cid (@ingroupNames){
		my @pps;
		foreach my $pid (keys %primer_listing){
			push(@pps, $pid) if grep(/^\Q$cid\E$/, @{$ppSeqHits{$pid}});
		}
		add_ornaments_pp(\@pps,$cid);
	}
}
############################################################################################################################
sub pattern{
	my $seq = shift;
	#foreach my $amb (keys %amb_tbl){
	#	$seq =~ s/$amb/\[$amb_tbl{$amb}\]/g;
	#}
	return $seq;
}
############################################################################################################################
sub init{
	getopts("t:o:p:n:d:m:b:c:s:g:e:a:z:y:u:k:x:i:j:f:r:l:w:q:h:S:T:E:D:O:");
	my $usage = "usage:
	$0
		-t <file: phylogenetic tree of the sequences in NEWICK format>
		-O <file: pairwise distance matrix of sequences>, default is internally computed by nw_distances
		-o <name of outgroup in above tree> leave out to use tree rooted as is
		-k [topdown|bottomup|evaluate|degentm]
		-S [preOrdered|cladeSize|abundance] topdown clade priority, default cladeSize
		-p <file: protein multi FASTA>
		-n <file: nucleotide multi FASTA>
		-z <file: extra seqs multi FASTA>
		-d <output directory name>
		-m <file: primer3 parameter template>
		-b <files: blast DB of target/include regions> (or else see -h)
		-h <start-end> only with -q prealigned: boundaries of region to search for primers, nuc aln coordinates (or else see -b)
		-s [include|target] primer3 search strategy 
		-g [genafpair|globalpair] mafft params, default is genafpair
		-a <file: abundances of each input sequence>
		-y <max degeneracy of primers>
		-D [*|+] operation to compute combined degeneracy
		-e <Min Tm threshold> primer3 can't compute Tm for degenerate oligos, UjiLity does and can apply a threshold
		-E <Max Tm threshold> primer3 can't compute Tm for degenerate oligos, UjiLity does and can apply a threshold
		-T <max Tm difference> between primers of a pair, default 4
		-u <number of unsuccessful sequences to try before closing primer pair> bottom up only, default is 50
		-x <max number of primer3 primer pairs proposals to evaluate> default 200
		-i <number of allowed missmatches when evaluating primers> default 0
		-j <threshold for consensus building> default 100
		-w [0|1] ignore columns with gaps for consensus building, default 0
		-f <file of forward primers to evaluate> must also specify -r -d -k degentm
		-r <file of reverse primers to evaluate> must also specify -f -d -k degentm
		-l [0|1|2] debug mode (beware, very verbose!), default 0
		-q [align|prealigned] either align prots & thread nuc (-p needed), or nuc file is prealigned, default is align
		-c [ifpncbmwzy] use cached results files for:
			i initialisations
			f fasta file extractions
			p protein aln
			n nucleotide threading
			c consensus building
			b blastx searching
			m primer3 runs
			w dreg searches
			zy extra seq computes
";
	$opt_j = 100 		unless defined($opt_j);
	$opt_i = 0 		unless defined($opt_i);
	$opt_c = '' 		unless defined($opt_c);
	#$opt_k = 'bottomup' 	unless defined($opt_k);
	$opt_u = 50 		unless defined($opt_u);
	$opt_T = 4 		unless defined($opt_T);
	$opt_x = 200 		unless defined($opt_x);
	$opt_o = '' 		unless defined($opt_o);
	$opt_O = '' 		unless defined($opt_O);
	$opt_h = '' 		unless defined($opt_h);
	$opt_g = 'genafpair' 	unless defined($opt_g);
	$opt_S = 'cladeSize' 	unless defined($opt_S);
	$opt_q = 'align' 	unless defined($opt_q);
	$opt_w = 0 		unless defined($opt_w);
	$opt_D = '*' 		unless defined($opt_D);
	$debug = $opt_l if (defined $opt_l && $opt_l);	
	if ( !(
		defined($opt_t) && defined($opt_k) &&
		defined($opt_n) && defined($opt_p) && 
		defined($opt_d) && defined($opt_m) &&
		defined($opt_h) && #defined($opt_c) && 
		defined($opt_s) && #defined($opt_a) && 
		defined($opt_z) && defined($opt_y) 
	 ) and !(
		defined($opt_z) && defined($opt_h) &&
		defined($opt_r) && defined($opt_f) &&
		defined($opt_d) && defined($opt_k) &&
		defined($opt_n) && defined($opt_p) &&
		defined($opt_t)
	   ) and !(
	   	defined($opt_r) && defined($opt_f) &&
	   	defined($opt_d) 
	   ) ) {
			print STDERR $usage;
			exit();
	}
	($runID) = ($opt_d =~ /([a-zA-Z0-9\-\._]+)$/);
	#work directory for all output
	do_cmd("mkdir -p $opt_d") unless -d $opt_d;
	# redirect STDERR to log file
	open (STDERR, ">$opt_d/UjiLity_run_output.txt") or die "OOOOOps $!";
	say STDERR "UjiLity version ".VERSION;
	say STDERR "Parameters used:";
	foreach my $arg ('t', 'O', 'o', 'p', 'n', 'd', 'b', 'c', 's', 'g', 'e', 'E', 'a', 'z', 'y', 'D', 'u', 'k', 'm', 'x', 'i', 'j', 'r', 'f', 'l','w','q','h', 'S', 'T') {
		say STDERR "$arg\t: ",(defined(${"opt_$arg"})?${"opt_$arg"}:'');
	}
	say STDERR "runID\t: $runID";
	if ($opt_k ne 'degentm'){
		if ($opt_i) {
			$dreg .= " -pmismatch $opt_i";
		}
		$stats{'clade_total'}=0;
		$stats{'clade_positives'}=0;
		$stats{'clade_positives_leaves_total'}=0;
		$stats{'clade_negatives'}=0;
		$stats{'clade_negatives_leaves_total'}=0;
		$stats{'clade_positives_list'}='';
		$stats{'clade_negatives_list'}='';
		write_lua_file();
		if (! -d dirname("$opt_d/$opt_t.rerooted")){
			do_cmd("mkdir -p ".dirname("$opt_d/$opt_t.rerooted"));
		}
		if ( defined($opt_o) && $opt_o ) {
			#reroot the tree 
			do_cmd_cached("$nw_utils/nw_reroot $opt_t $opt_o > $opt_d/$opt_t.rerooted","$opt_d/$opt_t.rerooted",'i');
		}else{
			do_cmd_cached("cp $opt_t $opt_d/$opt_t.rerooted","$opt_d/$opt_t.rerooted",'i');
		}	
		#rename the internal nodes
		do_cmd_cached("$nw_utils/nw_luaed -f $opt_d/number_inodes.lua $opt_d/$opt_t.rerooted > $opt_d/$opt_t.rerooted.renamed","$opt_d/$opt_t.rerooted.renamed",'i');
		if ($opt_k ne 'evaluate') {
			#$manager = new Parallel::ForkManager( 8 );
			#transform into cladogram to resolve issues with zero length branches (~trifurcation pb for binary division)
			do_cmd_cached("$nw_utils/nw_topology $opt_d/$opt_t.rerooted.renamed > $opt_d/$opt_t.rerooted.renamed.cladogram","$opt_d/$opt_t.rerooted.renamed.cladogram",'i');
			if ( defined($opt_o) && $opt_o ) {
				#reroot the tree 
				do_cmd_cached("$nw_utils/nw_reroot $opt_t $opt_o > $opt_d/$opt_t.rerooted","$opt_d/$opt_t.rerooted",'i');
				#get ingroup = sister clade to outgroup:
				do_cmd_cached("$nw_utils/nw_clade -s $opt_d/$opt_t.rerooted.renamed.cladogram $opt_o > $opt_d/$opt_t.ingroup","$opt_d/$opt_t.ingroup",'i');
				do_cmd_cached("$nw_utils/nw_clade -s $opt_d/$opt_t.rerooted.renamed $opt_o > $opt_d/$opt_t.ingroup_phylogram","$opt_d/$opt_t.ingroup_phylogram",'i');
			}else{
				#tree is midpoint rooted, so no outgroup to exclude
				
				do_cmd_cached("cp $opt_d/$opt_t.rerooted.renamed.cladogram $opt_d/$opt_t.ingroup","$opt_d/$opt_t.ingroup",'i');
				do_cmd_cached("cp $opt_d/$opt_t.rerooted.renamed $opt_d/$opt_t.ingroup_phylogram","$opt_d/$opt_t.ingroup_phylogram",'i');
			}
			if ($opt_O && -f $opt_O){
				do_cmd("cp $opt_O $opt_d/$opt_t.ingroup_phylogram.distances");
			}
			unless (-f "$opt_d/$opt_t.ingroup_phylogram.distances"){
				do_cmd_cached("$nw_utils/nw_distance -s f -n -mm $opt_d/$opt_t.ingroup_phylogram > $opt_d/$opt_t.ingroup_phylogram.distances","$opt_d/$opt_t.ingroup_phylogram.distances",'i');
			}
		#				SEQ_1	SEQ_2	SEQ_3
		#SEQ_1	0				2.32945			2.32945
		#SEQ_2			2.32945				0			0.00000007
		#SEQ_3				2.32945				0.00000007		0

			open (IN,"$opt_d/$opt_t.ingroup_phylogram.distances") || die "$!" ;
			my $lnb = 0;my @heads;
			while (my $line = <IN>){
				chomp($line);
				$lnb++;
				if ($lnb == 1){
					@heads = split /\t/,$line;
				}else{
					my ($id, @rows) = split /\t/,$line;
					my $rownb = 0;
					foreach my $d (@rows){
						$rownb++;
						$dist{$id}{$heads[$rownb]}=$d;
						#print STDERR "dist{$id}{$heads[$rownb]}=$d\n";
					}
				}
			}
			close IN;	
			@ingroupNames = get_leaves("$opt_d/$opt_t.ingroup");
			#do_cmd("cp $opt_d/$opt_t.ingroup $opt_d/$opt_t.annotated.nw");
			if (-e "$opt_d/$opt_t.ep.dreg.eval.primer.stats"){
				print STDERR "\nReading in the external environmental assessment counts\n";
				open (IN,"$opt_d/$opt_t.ep.dreg.eval.primer.stats") || die "$!" ;
				while (my $line = <IN>){
					chomp($line);
					#N102	3	0	551	0
					my ($primerid,$hits,$null,$counts,$nulll) = split /\t/,$line;
					$envHits{$primerid}=$hits;
					$envCounts{$primerid}=$counts;
				}
				close IN;
			}
		}else{
			#$opt_t = 'no_tree';
			if (defined($opt_f) && -e $opt_f){
				do_cmd("cp $opt_f $opt_d/$opt_t.primers.F.fasta") 
			}elsif( defined($opt_f) ){
				die "$opt_f not found: $!";
			}
			if (defined($opt_r) && -e $opt_r){
				do_cmd("cp $opt_r $opt_d/$opt_t.primers.R.fasta") 
			}elsif( defined($opt_r) ){
				die "$opt_r not found: $!";
			}
			#do_cmd("$t_coffee -other_pg seq_reformat -in $opt_n -output code_name | cut -f1 -d\" \" > "$opt_d/$opt_t.ingroup"");
			if (-f $opt_z){
				do_cmd("$seqkit seq -ni $opt_z > $opt_d/primer_evaluation.$opt_z.ids");
				#####################################################
				open(IND, "$opt_d/primer_evaluation.$opt_z.ids") or die "$!";
				while (my $line = <IND>){
					chomp($line);
					push(@ingroupNames,$line);
				}
				close IND;
			}
		}
		if (defined($opt_a)){
			die "Where is the abundance file $opt_a : $!" unless -e $opt_a;
			print STDERR "\nReading in the abundances/weights of the DNA sequences\n";
			open IN,$opt_a ;
			#my $min = 1e12;
			#TARA_004_DCM_0.22-1.6_G_scaffold136330_1_gene107070	367.413047879777
			while (my $line = <IN>){
				chomp($line);
				my ($ctg,$abun) = split /\s+/,$line;
				$abundances{$ctg} = $abun;
				#$min = $abun if $min > $abun;
			}
			close IN;
			# Take care of missing abundances
			my $miss = 0;
			foreach my $id (@ingroupNames){
				unless ( exists($abundances{$id}) ) {
					$abundances{$id} = 0;#$min/2;
					$miss++;
				}
			}
			print STDERR "\tSet $miss missing abundances/weights to zero.\n" if $miss;
		}else{
			print STDERR "\nNo file provided for abundances/weights of the DNA sequences\n";
			my $miss = 0;
			foreach my $id (@ingroupNames){
				$abundances{$id} = 1;
				$miss++;
			}
			print STDERR "\tSet $miss abundances/weights to 1.\n" if $miss;
		}

		#my $color_scheme = Color::Scheme->new();
		#$color_scheme->variation('hard');#soft, light,hard,pale
		#$color_scheme->scheme('analogic');#contrast, tetrade, analogic,mono,triade
		#$color_scheme->distance(1);
	    	#$color_scheme->web_safe(1);
	    	#$color_scheme->from_hue(0);
	    	#$color_scheme->add_complement( 0 );
		#@colors = $color_scheme->colors();
		#add_ornaments_REFS(@ingroupNames);
		@colors = ("FFFF00", "1CE6FF", "FF34FF", "FF4A46", "008941", "006FA6", "A30059", "7A4900", "63FFAC", "B79762", "004D43", "8FB0FF", "997D87",
"5A0007", "809693", "FEFFE6", "1B4400", "4FC601", "3B5DFF", "FF2F80",
"61615A", "BA0900", "6B7900", "00C2A0", "FFAA92", "FF90C9", "B903AA", "D16100",
"DDEFFF", "000035", "7B4F4B", "A1C299", "300018", "0AA6D8", "013349", "00846F",
"372101", "FFB500", "C2FFED", "A079BF", "CC0744", "C0B9B2", "C2FF99", "001E09",
"00489C", "6F0062", "0CBD66", "EEC3FF", "456D75", "B77B68", "7A87A1", "788D66",
"885578", "FAD09F", "FF8A9A", "D157A0", "BEC459", "456648", "0086ED", "886F4C",
"34362D", "B4A8BD", "00A6AA", "452C2C", "636375", "A3C8C9", "FF913F", "938A81",
"575329", "00FECF", "B05B6F", "8CD0FF", "3B9700", "04F757", "C8A1A1", "1E6E00",
"7900D7", "A77500", "6367A9", "A05837", "6B002C", "772600", "D790FF", "9B9700",
"549E79", "FFF69F", "201625", "72418F", "BC23FF", "99ADC0", "3A2465", "922329",
"5B4534", "FDE8DC", "404E55", "0089A3", "CB7E98", "A4E804", "324E72", "6A3A4C");
	}
	$amb_tbl{"M"}= "AC"; 	
	$amb_tbl{"R"}= "AG"; 	
	$amb_tbl{"W"}= "AT"; 	
	$amb_tbl{"S"}= "CG"; 	
	$amb_tbl{"Y"}= "CT"; 	
	$amb_tbl{"K"}= "GT"; 	
	$amb_tbl{"V"}= "ACG"; 	
	$amb_tbl{"H"}= "ACT"; 	
	$amb_tbl{"D"}= "AGT"; 	
	$amb_tbl{"B"}= "CGT"; 	
	$amb_tbl{"N"}= "ATGC";
}


sub pal2nal{
#
#    pal2nal.pl  (v14)                                      Mikita Suyama
#	adapted here by Pascal Hingamp for UjiLity
#    Usage:  pal2nal.pl  pep.aln  nuc.fasta  [nuc.fasta...]  [options] 
#
#        pep.aln:    protein alignment either in CLUSTAL or FASTA format
#
#        nuc.fasta:  DNA sequences (single multi-fasta or separated files)
#
#        Options:  -output (clustal(default)|paml|fasta|codon)
#		   -oufile FILENAME
#                  -nogap
#                  -html
#                  -nostderr
#                  -blockonly
#                  -nomismatch
#                  -codontable N (default=1(universal))
#
#        - IDs in pep.aln are used in the output.
#
#        - Sequence order is automatically checked (see the comment of v12).
#
#        - In-frame stop -> "*" or "_"
#
#        - Frameshift
#
#            Gene    HCDG
#            Pseudo  HC1G     2 del in pseudo
#
#            Gene    HCDG
#            Pseudo  HC2G     1 del in pseudo
#
#            Gene    ND-TY
#            Pseudo  ND1TY    1 ins in pseudo
#
#            Gene    ND-TY
#            Pseudo  ND2TY    2 ins in pseudo
#
#            Gene    EREQK
#            Pseudo  EK4QK    1 ins in pseudo
#

	no strict "vars";
	#$| = 1;
	my $recSep= $/;
	undef $/;
	$myos = $^O;
	my @ARGV = @_;
	if ($#ARGV < 1) {
	    #& showhelp();
	    exit;
	} else {
	    $getoutform = 0;
	    $nogap = 0;
	    $html = 0;
	    $nostderr = 0;
	    #$fasta = 0;
	    undef($alnfile);
	    undef(@nucfiles);
	    $outform = "clustal";
	    $blockonly = 0;
	    $nomismatch = 0;
	    $getcodontable = 0;
	    $codontable = 1;
	    foreach $i (0..$#ARGV) {
		if ($ARGV[$i] eq "-h") {
		    & showhelp();
		} elsif ($ARGV[$i] eq "-output") {
		    $getoutform = 1;
		} elsif ($getoutform) {
		    $outform = $ARGV[$i];
		    if ($outform ne "clustal" &&
		        $outform ne "paml" &&
		        $outform ne "fasta" &&
		        $outform ne "codon") {
		        print STDERR "\nERROR:  valid output format: clustal, paml, fasta, or codon\n\n";
		        exit;
		    }
		    $getoutform = 0;
		} elsif ($ARGV[$i] eq "-outfile") {
		    $getoutfile = 1;
		} elsif ($getoutfile) {
		    $outfilename = $ARGV[$i];
		    if (-f $outfilename) {
		        print STDERR "\nWARNING:  outfile $outfilename exists!\n";
		        #exit;
		    }
		    $getoutfile = 0;
		} elsif ($ARGV[$i] eq "-blockonly") {
		    $blockonly = 1;
		} elsif ($ARGV[$i] eq "-nogap") {
		    $nogap = 1;
		} elsif ($ARGV[$i] eq "-nomismatch") {
		    $nomismatch = 1;
		} elsif ($ARGV[$i] eq "-codontable") {
		    $getcodontable = 1;
		} elsif ($getcodontable) {
		    $codontable = $ARGV[$i];
		    if ($codontable != 1 && $codontable != 2 && $codontable != 3 &&
		        $codontable != 4 && $codontable != 5 && $codontable != 6 &&
		        $codontable != 9 && $codontable != 10 && $codontable != 11 &&
		        $codontable != 12 && $codontable != 13 && $codontable != 14 &&
		        $codontable != 15 && $codontable != 16 && $codontable != 21 &&
		        $codontable != 22 && $codontable != 23) {
		        print STDERR "\nERROR:  invalid codontable number, $codontable!!\n\n";
		        exit;
		    }
		    $getcodontable = 0;
		} elsif ($ARGV[$i] eq "-html") {
		    $html = 1;
		} elsif ($ARGV[$i] eq "-nostderr") {
		    $nostderr = 1;
		} elsif (!$alnfile) {
		    $alnfile = $ARGV[$i];
		} else {
		    push(@nucfiles, $ARGV[$i]);
		}
	    }
	}
	open (my $OUT, '>', "$outfilename");
	select $OUT;
	if ($html) {
	    print "<pre>\n";
	}

	#---------------#
	# check options  ("-outform codon" is not valid with -blockonly, -nogap, -nomismatch.)
	#---------------#

	if ($outform eq "codon" &&
	    ($blockonly || $nogap || $nomismatch)) {
	    if ($html) {
		print "\nERROR:  if \"codon(Output format)\" is selected, don't use \"Remove gaps, inframe stop codons\" or \"Remove mismatches\" or \"Use only selected positions\".\n\n";
	    } else {
		print STDERR "\nERROR:  \"-outform codon\" is not valid with -blockonly, -nogap, -nomismatch\n\n";
	    }
	    exit;
	}


	#---------------------#
	#  Get nuc sequences
	#---------------------#

	undef(@nucid);
	undef(%id2nucseq);
	$nseq = -1;
	foreach $i (0..$#nucfiles) {
	    open(NUCFILE, "< $nucfiles[$i]") || die "Can't open $nucfiles[$i]";
	    $nucfiledata = <NUCFILE>;
	    close(NUCFILE);
	    $nucfiledata =~ s/\x0D\x0A|\x0D|\x0A/\n/g;
	    foreach (split(/\n/, $nucfiledata)) {
		if (!/^#/ && /\S+/) {
		    if (/^>(\S+)/) {
		        ++$nseq;
		        $tmpnucid = $1;
		        push(@nucid, $tmpnucid);
		    } else {
		        s/[^a-zA-Z]//g;
		        $id2nucseq{$tmpnucid} .= $_;
		    }
		}
	    }
	}


	#-------------------#
	#  Get aa alignemt
	#-------------------#

	undef(@aaid);
	undef(%id2aaaln);
	undef(%aaidcnt);
	undef(@aaseq);
	undef($gblockseq);

	$gettype = 1;
	open(ALNFILE, "< $alnfile") || die "Can't open $alnfile";
	while (<ALNFILE>) {
	    chomp;
	    if ($gettype && !/^#/ && /\S+/) {
		if (/^CLUSTAL/) {
		    $inalntype = "clustal";
		} elsif (/^>/) {
		    $inalntype = "fasta";
		} elsif (/^Gblocks/) {
		    $inalntype = "gblocks";
		} else {
		    $inalntype = "clustal";
		}

		$gettype = 0;
	    }
	}
	close(ALNFILE);

	open(ALNFILE, "< $alnfile") || die "Can't open $alnfile";
	$getblock = 0;
	$aafiledata = <ALNFILE>;
	close(ALNFILE);
	$aafiledata =~ s/\x0D\x0A|\x0D|\x0A/\n/g;
	foreach (split(/\n/, $aafiledata)) {
	    if ($inalntype eq "clustal") {
		if (!/^CLUSTAL/ && /^\S+/ && !/^#/) {
		    s/\s+$//;
		    @dat = split(/\s+/, $_);
		    ++$aaidcnt{$dat[0]};
		    push(@aaid, $dat[0]) if ($aaidcnt{$dat[0]} == 1);
		    $dat[1] =~ tr/a-z/A-Z/;
		    $id2aaaln{$dat[0]} .= $dat[1];

		    $tmplen = length($_);

		    /^\S+\s+/;
		    $idspc = length($&);
		    $subalnlen = length($dat[1]);

		    $getblock = 1;
		} elsif ($getblock) {
		    $_ .= ' ' x ($tmplen - length($_));
		    $gblockseq .= substr($_, $idspc, $subalnlen);

		    $getblock = 0;
		}
	    } elsif ($inalntype eq "fasta") {
		if (/^>(\S+)/) {
		    $tmpid = $1;
		    push(@aaid, $tmpid);
		} else {
		    s/\s+//g;
		    tr/a-z/A-Z/;
		    $id2aaaln{$tmpid} .= $_;
		}
	    } elsif ($inalntype eq "gblocks") {
		$getaln = 0;
		if (/^\s+\=/) {
		    $getaln = 1;
		} elsif (/^Parameters/) {
		    $getaln = 0;
		} elsif ($getaln) {
		    @dat = split(/\s+/, $_);
		    if (/^Gblocks/) {
		        $gblockseq .= $dat[1];
		    } elsif (/^\S+/) {
		        ++$aaidcnt{$dat[0]};
		        push(@aaid, $dat[0]) if ($aaidcnt{$dat[0]} == 1);
		        $dat[1] =~ tr/a-z/A-Z/;
		        $id2aaaln{$dat[0]} .= $dat[1];
		    }
		}
	    }
	}

	foreach $i (0..$#aaid) {
	    push(@aaseq, $id2aaaln{$aaid[$i]});
	}


	#------------------------------------#
	#  For frameshifts in tfastx/tfasty
	#------------------------------------#

	foreach $i (0..$#aaid) {
	    $aaseq[$i] =~ s/\\/1/g;
	    $aaseq[$i] =~ s/\/(-*)[A-Z\*]/-${1}2/g;
	}


	#-------------------------------------#
	#  Check the input seqs (pep <=> nuc)
	#-------------------------------------#

	if ($#aaid != $#nucid) {
	    $naa = $#aaid + 1;
	    $nnuc = $#nucid + 1;
	    if ($html) {
		print "\nERROR: number of input seqs differ (aa: $naa;  nuc: $nnuc)!!\n\n";
	    } else {
		print STDERR "\nERROR: number of input seqs differ (aa: $naa;  nuc: $nnuc)!!\n\n";
		print STDERR "   aa  '@aaid'\n";
		print STDERR "   nuc '@nucid'\n";
	    }
	    exit;
	}


	#---------------------------------#
	# corresponding IDs or same order
	#---------------------------------#

	@commonids = & common_elem(*aaid, *nucid);

	sub common_elem {
	    local(*aarr, *barr) = @_;

	    local(%mark, @comarr);
	    grep($mark{$_}++, @aarr);
	    @comarr = grep($mark{$_}, @barr);
	}
		
	if ($#commonids == $#aaid) {
	    $idcorrespondence = "sameID";
	} else {
	    $idcorrespondence = "ordered";
	}


	#-------------------#
	#  Codon sequences
	#-------------------#

	undef(@codonseq);
	undef(%aaidpos2mismatch);
	undef(@outmessage);
	foreach $i (0..$#aaid) {
	    if ($idcorrespondence eq "sameID") {
		$tmpnucid = $aaid[$i];
	    } elsif ($idcorrespondence eq "ordered") {
		$tmpnucid = $nucid[$i];
	    } else {
		print STDERR "\nERROR in ID correspondence.\n\n";
		exit;
	    }

	    %codonout = & pn2codon($aaseq[$i], $id2nucseq{$tmpnucid}, $codontable);
	    @message = @{$codonout{'message'}} if defined $codonout{'message'};

	    foreach $j (0..$#message) {
		$outl = "WARNING: $aaid[$i] $message[$j]";
		push(@outmessage, $outl);
		@dat = split(/\s+/, $message[$j]);
		$dat[1] =~ s/:$//;
		$aaidpos2mismatch{"$aaid[$i] $dat[1]"} = 1;
	    }

	    if ($codonout{'result'} == 1 || $codonout{'result'} == 2) {
		push(@codonseq, $codonout{'codonseq'});
	    } else {
		$erraa = $aaseq[$i];
		$erraa =~ s/-//g;
		# 1 while $erraa =~ s/(.{60})(.+)/$1\n$2/;
		@frags = & my1while($erraa, 60);
		$erraa = join("\n", @frags);

		$errnuc = $id2nucseq{$tmpnucid};
		$errnuc =~ s/-//g;
		# 1 while $errnuc =~ s/(.{60})(.+)/$1\n$2/;
		@frags = & my1while($errnuc, 60);
		$errnuc = join("\n", @frags);

		if ($html) {
		    print "</pre>\n";
		    print "<H1>ERROR in your input files!</H1>\n";
		    print "<pre>\n";
		    print "#---  ERROR: inconsistency between the following pep and nuc seqs  ---#\n";

		    print ">$aaid[$i]\n";
		    print "$erraa\n";

		    print ">$tmpnucid\n";
		    print "$errnuc\n";

		    print "</pre>\n";
		    $tmpwhich = "tmp/tmpwhich.$$";
		} else {
		    print STDERR "#---  ERROR: inconsistency between the following pep and nuc seqs  ---#\n";

		    print STDERR ">$aaid[$i]\n";
		    print STDERR "$erraa\n";

		    print STDERR ">$tmpnucid\n";
		    print STDERR "$errnuc\n";

		    $tmpwhich = "tmpwhich.$$";
		}

		if ($myos eq "linux") {
		    system("which bl2seq > $tmpwhich");

		    $foundbl2seq = 0;
		    open(TMPWHICH, "< $tmpwhich") || die "Can't open $tmpwhich";
		    while (<TMPWHICH>) {
		        chomp;
		        $foundbl2seq = 1 if (/bl2seq$/);
		    }
		    close(TMPWHICH);
		    unlink($tmpwhich);

		    if ($foundbl2seq) {


		        #------------------------------------------#
		        # run bl2seq (if the command is available)
		        #------------------------------------------#

		        if ($html) {
		            $erraafile = "tmp/erraafile.$$";
		            $errnucfile = "tmp/errnucfile.$$";
		            $tblnout = "tmp/tbln.out.$$";
		        } else {
		            $erraafile = "erraafile.$$";
		            $errnucfile = "errnucfile.$$";
		            $tblnout = "tbln.out.$$";
		        }

		        open(ERRAAFILE, "> $erraafile") || die "Can't open $erraafile";
		        print ERRAAFILE ">$aaid[$i]\n";
		        print ERRAAFILE "$erraa\n";
		        close(ERRAAFILE);

		        open(ERRNUCFILE, "> $errnucfile") || die "Can't open $errnucfile";
		        print ERRNUCFILE ">$tmpnucid\n";
		        print ERRNUCFILE "$errnuc\n";
		        close(ERRNUCFILE);

		        system("bl2seq -p tblastn -F F -i $erraafile -j $errnucfile -o $tblnout");

		        if ($html) {
		            print "<BR>\n";
		            print "<BR>\n";
		            print "<H1>Check the following TBLASTN output.</H1><BR>\n";
		            print "<pre>\n";
		            print "      your pep -> 'Query'\n";
		            print "      your nuc -> 'Sbjct'\n";
		            print "<BR>\n";
		        } else {
		            print STDERR "\n";
		            print STDERR "\n";
		            print STDERR "        ###-----   Check the following TBLASTN output:           -----###\n";
		            print STDERR "        ###-----       your pep -> 'Query'                       -----###\n";
		            print STDERR "        ###-----       your nuc -> 'Sbjct'                       -----###\n";
		            print STDERR "\n";
		        }

		        open(TBLNOUT, "< $tblnout") || die "Can't open $tblnout";
		        $tblnoutdata = <TBLNOUT>;
		        close(TBLNOUT);
		        $tblnoutdata =~ s/\x0D\x0A|\x0D|\x0A/\n/g;
		        foreach (split(/\n/, $tblnoutdata)) {
		            if ($html) {
		                print "$_\n";
		            } else {
		                print STDERR "$_\n";
		            }
		        }

		        if ($html) {
		            print "</pre>\n";
		        }

		        unlink($erraafile);
		        unlink($errnucfile);
		        unlink($tblnout);
		    } else {
		        if ($html) {
		            print "<BR>\n";
		            print "<BR>\n";
		            print "Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.<BR>\n";
		            print "<BR>\n";
		        } else {
		            print STDERR "\n";
		            print STDERR "\n";
		            print STDERR "Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.\n";
		            print STDERR "\n";
		        }
		    }
		} else {    ####    non-linux environment
		    print STDERR "\n";
		    print STDERR "\n";
		    print STDERR "Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.\n";
		    print STDERR "\n";
		}
		exit;
	    }
	}


	#-------------------#
	#  Warning in '#'?
	#-------------------#

	if (defined $gblockseq && $gblockseq =~ /#/ && $blockonly) {
	    undef(@newoutmessage);
	    foreach $i (0..$#outmessage) {
		$mpos = (split(/\s+/, $outmessage[$i]))[3];
		$mpos =~ s/://;

		if (substr($gblockseq, $mpos - 1, 1) eq "#") {
		    push(@newoutmessage, $outmessage[$i]);
		}
	    }
	    @outmessage = @newoutmessage;
	}


	#--------------------#
	#  Warning messages
	#--------------------#

	if ($codontable != 1) {
	    # push(@outmessage, "Codon table $codontable is used");
	    splice(@outmessage, 0, 0, "Codontable $codontable is used");
	}

	if (!$nostderr) {
	    foreach $j (0..$#outmessage) {
		if ($html) {
		    if ($j == 0 && !$nomismatch) {
		        print "#------------------------------------------------------------------------#\n";
		    }
		    if ($nomismatch) {
		        # print "#  $outmessage[$j]  (excluded from the output)\n";
		    } else {
		        print "#  $outmessage[$j]\n";
		    }
		    if ($j == $#outmessage && !$nomismatch) {
		        print "#------------------------------------------------------------------------#\n\n";
		    }
		} else {
		    if ($j == 0 && !$nomismatch) {
		        print STDERR "#------------------------------------------------------------------------#\n";
		        print STDERR "#  Input files:  $alnfile @nucfiles\n";
		    }
		    if ($nomismatch) {
		        # print STDERR "#  $outmessage[$j]  (exlucded from the output)\n";
		    } else {
		        print STDERR "#  $outmessage[$j]\n";
		    }
		    if ($j == $#outmessage && !$nomismatch) {
		        print STDERR "#------------------------------------------------------------------------#\n\n";
		    }
		}
	    }
	}

	undef(%errorpos);
	foreach $j (0..$#outmessage) {
	    @dat = split(/\s+/, $outmessage[$j]);
	    $tmperrpos = $dat[3];
	    $tmperrpos =~ s/:$//;
	    $tmperrpos -= 1;
	    $errorpos{$tmperrpos} = 1;
	}


	#----------------------------------#
	#  Make an AA-based NUC alignment
	#----------------------------------#

	undef(@tmppos);
	undef(@codonaln);
	undef(@coloraln);
	undef($maskseq);
	foreach $i (0..length($aaseq[0]) - 1) {
	    $tmpmax = 0;
	    $apos = $i + 1;

	    #-----------#
	    # gblocks ?
	    #-----------#

	    $putcodon = 1;
	    if (defined $gblockseq && $gblockseq =~ /#/ && substr($gblockseq, $i, 1) ne "#" && $blockonly) {
		$putcodon = 0;
	    }
	    if ($nomismatch && $errorpos{$i}) {
		$putcodon = 0;
	    }

	    foreach $k (0..$#aaid) {
		$tmpaa = substr($aaseq[$k], $i, 1);
		if ($tmpaa !~ /\d/) {
		    $tmplen = 3;
		} else {
		    $tmplen = (int(($tmpaa - 1) / 3) + 1) * 3;  # 1, 2, 3 -> 3
		                                                # 4, 5, 6 -> 6
		                                                # 7, 8, 9 -> 9
		}
		$tmpmax = $tmplen if ($tmpmax < $tmplen);
	    }
	    foreach $k (0..$#aaid) {
		$tmpaa = substr($aaseq[$k], $i, 1);
		$codonaln[$k] = '' unless defined $codonaln[$k];
		if ($tmpaa !~ /\d/) {
		    if ($tmpaa eq '-') {
		        # if ($putcodon || (!$blockonly && !$nomismatch)) {
		        if ($putcodon) {
		            $codonaln[$k] .= '-' x $tmpmax;
		            if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
		                $coloraln[$k] .= 'R' x $tmpmax;
		            } else {
		                $coloraln[$k] .= '-' x $tmpmax;
		            }
		        }
		    } elsif ($tmpaa =~ /[A-Z]/ || $tmpaa eq '*') {
		        # if ($putcodon || (!$blockonly && !$nomismatch)) {
		        if ($putcodon) {
		            if (defined $codonseq[$k] && defined $tmppos[$k]){
		            	$codonaln[$k] .= substr($codonseq[$k], $tmppos[$k], 3);
		            }
		            if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
		                $coloraln[$k] .= 'RRR';
		            } else {
		                $coloraln[$k] .= '---';
		            }
		        }
		        $tmppos[$k] += 3;
		        # if ($putcodon || (!$blockonly && !$nomismatch)) {
		        if ($putcodon) {
		            $codonaln[$k] .= '-' x ($tmpmax - 3) if ($tmpmax - 3 > 0);
		            if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
		                $coloraln[$k] .= 'R' x ($tmpmax - 3) if ($tmpmax - 3 > 0);
		            } else {
		                $coloraln[$k] .= '-' x ($tmpmax - 3) if ($tmpmax - 3 > 0);
		            }
		        }
		    }
		} elsif ($tmpaa =~ /\d/) {
		    # if ($putcodon || (!$blockonly && !$nomismatch)) {
		    if ($putcodon) {
		        $codonaln[$k] .= substr($codonseq[$k], $tmppos[$k], $tmpaa);
		        if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
		            $coloraln[$k] .= 'R' x $tmpaa;
		        } else {
		            $coloraln[$k] .= '-' x $tmpaa;
		        }
		    }
		    $tmppos[$k] += $tmpaa;
		    # if ($putcodon || (!$blockonly &&!$nomismatch)) {
		    if ($putcodon) {
		        $codonaln[$k] .= '-' x ($tmpmax - $tmpaa);
		        if ($aaidpos2mismatch{"$aaid[$k] $apos"}) {
		            $coloraln[$k] .= 'R' x ($tmpmax - $tmpaa);
		        } else {
		            $coloraln[$k] .= '-' x ($tmpmax - $tmpaa);
		        }
		    }
		}
	    }
	    if (!$blockonly) {
		# if ($putcodon && substr($gblockseq, $i, 1) eq "#") {
		#     $maskseq .= '#' x $tmpmax;
		# } else {
		#     $maskseq .= ' ' x $tmpmax;
		# }
		if ($putcodon) {
		    if (defined $gblockseq && substr($gblockseq, $i, 1) eq "#") {
		        $maskseq .= "#" x $tmpmax;
		    } else {
		        $maskseq .= " " x $tmpmax;
		    }
		}
	    }
	}


	#-----------#
	#  -nogap?
	#-----------#

	$alilen = length($codonaln[0]);

	if ($nogap) {
	    $tmppos = 0;
	    undef(@nogapaln);
	    undef(@nogapcoloraln);
	    undef($nogapmaskseq);
	    while ($tmppos < $alilen) {
		$outok = 1;
		foreach $i (0..$#codonaln) {
		    $tmpcodon = substr($codonaln[$i], $tmppos, 3);
		    $outok = 0 if ($tmpcodon =~ /-/);
		    $outok = 0 if ($tmpcodon =~ /(((U|T)A(A|G|R))|((T|U)GA))/);
		}
		if ($outok) {
		    foreach $i (0..$#codonaln) {
		        $tmpcodon = substr($codonaln[$i], $tmppos, 3);
		        $nogapaln[$i] .= $tmpcodon;
		        $tmpcolorcodon = substr($coloraln[$i], $tmppos, 3);
		        $nogapcoloraln[$i] .= $tmpcolorcodon;
		    }
		    $nogapmaskseq .= substr($maskseq, $tmppos, 3);
		}
		$tmppos += 3;
	    }
	    @codonaln = @nogapaln;
	    @coloraln = @nogapcoloraln;
	    $maskseq = $nogapmaskseq;
	}


	#----------#
	#  Output
	#----------#

	$maxn = 0;
	foreach $i (0..$#aaid) {
	    $maxn = length($aaid[$i]) if ($maxn < length($aaid[$i]));
	}
	$maxn = 10 if ($maxn < 10);

	$alilen = length($codonaln[0]);

	# foreach $i (0..$#aaid) {
	#     1 while $codonaln[$i] =~ s/(.{60})(.+)/$1\n$2/;
	#     1 while $coloraln[$i] =~ s/(.{60})(.+)/$1\n$2/;
	# }
	# 1 while $maskseq =~ s/(.{60})(.+)/$1\n$2/;

	undef(%aaid2codonarr);
	undef(%aaid2colorarr);
	foreach $i (0..$#aaid) {
	    @{$aaid2codonarr{$aaid[$i]}} = & my1while($codonaln[$i], 60);
	    @{$aaid2colorarr{$aaid[$i]}} = & my1while($coloraln[$i], 60);
	}
	@maskarr = & my1while($maskseq, 60);


	if ($outform eq "clustal") {

	    #-----------#
	    #  clustal
	    #-----------#

	    print "CLUSTAL W multiple sequence alignment\n";
	    print "\n";

	    @output1 = @{$aaid2codonarr{$aaid[0]}};
	    if ($html) {
		foreach $i (0..$#output1) {
		    foreach $k (0..$#aaid) {
		        printf "%-${maxn}s    ", $aaid[$k];
		        $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
		        $outr = ${$aaid2colorarr{$aaid[$k]}}[$i];
		        $rlen = length($outf);
		        if ($outr =~ /R/) {
		            foreach $l (0..$rlen - 1) {
		                $tmpnuc = substr($outf, $l, 1);
		                $tmpr = substr($outr, $l, 1);
		                if ($tmpr eq "R") {
		                    print "<FONT color='red'>$tmpnuc</FONT>";
		                } else {
		                    print "$tmpnuc";
		                }
		            }
		            print "\n";
		        } else {
		            print "$outf\n";
		        }
		    }
		    $outmask = $maskarr[$i];
		    printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
		    print "\n";
		}
	    } else {
		foreach $i (0..$#output1) {
		    foreach $k (0..$#aaid) {
		        $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
		        printf "%-${maxn}s    $outf\n", $aaid[$k];
		    }
		    $outmask = $maskarr[$i];
		    printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && defined $gblockseq && $gblockseq =~ /#/);
		    print "\n";
		}
	    }

	} elsif ($outform eq "codon") {

	    #-------#
	    # codon
	    #-------#

	    # @outaa = @aaseq;

	    undef(@outaa);

	    $withn = 0;
	    foreach $j (0..$#aaid) {
		$withn = 1 if ($aaseq[$j] =~ /\d/);
	    }

	    if ($withn) {
		$alnlen = length($aaseq[0]);
		foreach $i (0..$alnlen - 1) {
		    $maxaan = 0;
		    foreach $j (0..$#aaid) {
		        $tmpaa = substr($aaseq[$j], $i, 1);
		        if ($tmpaa =~ /\d/ && $tmpaa > $maxaan) {
		            $maxaan = $tmpaa;
		        }
		    }
		    if ($maxaan >= 4) {
		        $tmplen = int(($tmpaa - 1) / 3) + 1;  # 4, 5, 6 -> 2
		                                              # 7, 8, 9 -> 3
		    } else {
		        $tmplen = 1;
		    }
		    foreach $j (0..$#aaid) {
		        $pushaa = '-' x $tmplen;
		        $tmpaa = substr($aaseq[$j], $i, 1);
		        substr($pushaa, 0, 1) = $tmpaa;
		        $outaa[$j] .= $pushaa;
		    }
		}
	    } else {
		@outaa = @aaseq;
	    }

	    undef(%aaid2aaarr);
	    foreach $i (0..$#aaid) {
		# 1 while $outaa[$i] =~ s/(.{20})(.+)/$1\n$2/;
		@{$aaid2aaarr{$aaid[$i]}} = & my1while($outaa[$i], 20);
	    }

	    @output1 = @{$aaid2codonarr{$aaid[0]}};
	    if ($html) {
		foreach $i (0..$#output1) {
		    foreach $k (0..$#aaid) {
		        printf "%-${maxn}s    ", '';
		        $outa = ${$aaid2aaarr{$aaid[$k]}}[$i];
		        # 1 while $outa =~ s/(\S)(\S+)/$1   $2/;
		        @frags = & my1while($outa, 1);
		        $outa = join("   ", @frags);

		        $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
		        # 1 while $outf =~ s/(\S{3})(\S+)/$1 $2/;
		        @frags = & my1while($outf, 3);
		        $outf = join(" ", @frags);

		        $outr = ${$aaid2colorarr{$aaid[$k]}}[$i];
		        # 1 while $outr =~ s/(\S{3})(\S+)/$1 $2/;
		        @frags = & my1while($outr, 3);
		        $outr = join(" ", @frags);

		        if ($outr =~ /R/) {
		            $rlen = length($outf);
		            foreach $l (0..$rlen - 1) {
		                $tmppep = substr($outa, $l, 1);
		                $tmpr = substr($outr, $l, 1);
		                if ($tmpr eq "R") {
		                    print "<FONT color='red'>$tmppep</FONT>";
		                } else {
		                    print "$tmppep";
		                }
		            }
		            print "\n";
		            printf "%-${maxn}s    ", $aaid[$k];
		            foreach $l (0..$rlen - 1) {
		                $tmpnuc = substr($outf, $l, 1);
		                $tmpr = substr($outr, $l, 1);
		                if ($tmpr eq "R") {
		                    print "<FONT color='red'>$tmpnuc</FONT>";
		                } else {
		                    print "$tmpnuc";
		                }
		            }
		            print "\n";
		        } else {
		            print "$outa\n";
		            printf "%-${maxn}s    ", $aaid[$k];
		            print "$outf\n";
		        }
		    }
		    $outmask = $maskarr[$i];
		    $outmask =~ s/\s/-/g;
		    # 1 while $outmask =~ s/(\S{3})(\S+)/$1 $2/;
		    @frags = & my1while($outmask, 3);
		    $outmask = join(" ", @frags);
		    $outmask =~ s/-/ /g;
		    printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
		    print "\n";
		}
	    } else {
		foreach $i (0..$#output1) {
		    foreach $k (0..$#aaid) {
		        $outa = ${$aaid2aaarr{$aaid[$k]}}[$i];
		        # 1 while $outa =~ s/(\S)(\S+)/$1   $2/;
		        @frags = & my1while($outa, 1);
		        $outa = join("   ", @frags);
		        printf "%-${maxn}s    $outa\n", '';

		        $outf = ${$aaid2codonarr{$aaid[$k]}}[$i];
		        # 1 while $outf =~ s/(\S{3})(\S+)/$1 $2/;
		        @frags = & my1while($outf, 3);
		        $outf = join(" ", @frags);
		        printf "%-${maxn}s    $outf\n", $aaid[$k];
		    }
		    $outmask = $maskarr[$i];
		    $outmask =~ s/\s/-/g;
		    # 1 while $outmask =~ s/(\S{3})(\S+)/$1 $2/;
		    @frags = & my1while($outmask, 3);
		    $outmask = join(" ", @frags);
		    $outmask =~ s/-/ /g;
		    printf "%-${maxn}s    $outmask\n", ' ' if (!$blockonly && $gblockseq =~ /#/);
		    print "\n";
		}
	    }
	} elsif ($outform eq "paml") {

	    #--------#
	    #  paml
	    #--------#

	    $nseq = $#aaid + 1;
	    printf " %3d %6d\n", $nseq, $alilen;
	    foreach $i (0..$#aaid) {
		print  "$aaid[$i]\n";
		if ($html) {
		    @outf = @{$aaid2codonarr{$aaid[$i]}};
		    @outr = @{$aaid2colorarr{$aaid[$i]}};
		    foreach $k (0..$#outf) {
		        if ($outr[$k] =~ /R/) {
		            $lenf = length($outf[$k]);
		            foreach $l (0..$lenf - 1) {
		                $tmpnuc = substr($outf[$k], $l, 1);
		                $tmpr = substr($outr[$k], $l, 1);
		                if ($tmpr eq "R") {
		                    print "<FONT color='red'>$tmpnuc</FONT>";
		                } else {
		                    print "$tmpnuc";
		                }
		            }
		            print "\n";
		        } else {
		            print "$outf[$k]\n";
		        }
		    }
		} else {
		    $outcodon = join("\n", @{$aaid2codonarr{$aaid[$i]}});
		    print  "$outcodon\n";
		}
	    }
	} elsif ($outform eq "fasta") {

	    #---------#
	    #  fasta
	    #---------#

	    foreach $i (0..$#aaid) {
		print  ">$aaid[$i]\n";
		if ($html) {
		    @outf = @{$aaid2codonarr{$aaid[$i]}};
		    @outr = @{$aaid2colorarr{$aaid[$i]}};
		    foreach $k (0..$#outf) {
		        if ($outr[$k] =~ /R/) {
		            $lenf = length($outf[$k]);
		            foreach $l (0..$lenf - 1) {
		                $tmpnuc = substr($outf[$k], $l, 1);
		                $tmpr = substr($outr[$k], $l, 1);
		                if ($tmpr eq "R") {
		                    print "<FONT color='red'>$tmpnuc</FONT>";
		                } else {
		                    print "$tmpnuc";
		                }
		            }
		            print "\n";
		        } else {
		            print "$outf[$k]\n";
		        }
		    }
		} else {
		    $outcodon = join("\n", @{$aaid2codonarr{$aaid[$i]}});
		    print  "$outcodon\n";
		}
	    }
	}

	if ($html) {
	    print "</pre>\n";
	}
	$/ = $recSep;
	use strict "vars";
	select STDOUT;
}

#-----------------------------------------------------------------------

sub pn2codon {
    #    pn2codon v6
    #
    #    input:   $pep    protein sequence
    #                         termination -> "_" or "*";
    #                         frameshift  -> digit
    #                         "-" or "."  -> gap
    #             $nuc    DNA or RNA sequence (lower/upper case letters)
    #
    #             $codontable  (corresponds to codon tables used in GenBank
    #                1  Universal code
    #                2  Vertebrate mitochondrial code
    #                3  Yeast mitochondrial code
    #                4  Mold, Protozoan, and Coelenterate Mitochondrial code
    #                   and Mycoplasma/Spiroplasma code
    #                5  Invertebrate mitochondrial
    #                6  Ciliate, Dasycladacean and Hexamita nuclear code
    #                9  Echinoderm and Flatworm mitochondrial code
    #               10  Euplotid nuclear code
    #               11  Bacterial, archaeal and plant plastid code
    #               12  Alternative yeast nuclear code
    #               13  Ascidian mitochondrial code
    #               14  Alternative flatworm mitochondrial code
    #               15  Blepharisma nuclear code
    #               16  Chlorophycean mitochondrial code
    #               21  Trematode mitochondrial code
    #               22  Scenedesmus obliquus mitochondrial code
    #               23  Thraustochytrium mitochondrial code
    #
    #    return:  hash
    #                $retval{'codonseq'}: codon seq (w/o gap)
    #                $retval{'message'}:  error/warning messame (array)
    #                $retval{'result'}:   1: OK, 2: mismatch, -1: no match found
    #
    #
    #                                      05/05/2002    Mikita Suyama
    #                                      12/07/2004
    #


    my($pep, $nuc, $codontable) = @_;


    my(%p2c);
    my($peplen, $qcodon, $codon);
    my($tmpaa, $message, $modpep, @qcodon, @fncodon, $wholecodon);
    my($i, $j, $anclen, @anchor, $peppos, $tmpcodon, $codonpos);

    my(%retval);


    if ($codontable == 1) {
        #-----------#
        # Universal Code (transl_table=1)
        #-----------#
        %p2c = (
            "B" => "((U|T|C|Y|A)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 2) {
        #--------------------------#
        # Vertebrate Mitochondrial Code (transl_table=2)
        #--------------------------#
        %p2c = (
            "B" => "((A(U|T).)|G(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
            "*" => "(((U|T)A(A|G|R))|(AG(A|G|R)))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 3) {
        #--------------------------#
        # Yeast Mitochondrial Code (transl_table=3)
        #--------------------------#
        %p2c = (
            "B" => "(A(U|T)(A|G|R))",
            "L" => "((U|T)(U|T)(A|G|R))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "((AC.)|(C(U|T).))",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 4) {
        #-----------#
        # Mold, Protozoan, and Coelenterate Mitochondrial Code
        # and Mycoplasma/Spiroplasma Code (transl_table=4)
        #-----------#
        %p2c = (
            "B" => "((A(U|T).)|((U|T)(U|T)(A|G|R))|(C(U|T)G)|(G(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 5) {
        #-----------#
        # Invertebrate Mitochondrial Code (transl_table=5)
        #-----------#
        %p2c = (
            "B" => "((A(U|T).)|((U|T|A|G|R)(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 6) {
        #-----------#
        # Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((T|U)GA)",
            "*" => "((T|U)GA)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "((CA(A|G|R))|((U|T)A(A|G|R)))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 9) {
        #-----------#
        # Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
        #-----------#
        %p2c = (
            "B" => "((A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 10) {
        #-----------#
        # Euplotid Nuclear Code (transl_table=10)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 11) {
        #-----------#
        # Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
        #-----------#
        %p2c = (
            "B" => "((A(U|T)G)|(.(U|T)G))",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 12) {
        #-----------#
        # Alternative Yeast Nuclear Code (transl_table=12)
        #-----------#
        %p2c = (
            "B" => "((A|C)(U|T)G)",
            "L" => "((C(U|T)(U|T|C|Y|A))|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y))|(C(U|T)G))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 13) {
        #-----------#
        # Ascidian Mitochondrial Code (transl_table=13)
        #-----------#
        %p2c = (
            "B" => "((T|U|A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "((GG.)|(AG(A|G|R)))",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 14) {
        #-----------#
        # Alternative Flatworm Mitochondrial Code (transl_table=14)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)AG)",
            "*" => "((U|T)AG)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y|A))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)G(G|A|R))",
            "X" => "...",
        );
    } elsif ($codontable == 15) {
        #-----------#
        # Blepharisma Nuclear Code (transl_table=15)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)AA)|((T|U)GA))",
            "*" => "(((U|T)AA)|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "((CA(A|G|R))|((U|T)AG))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 16) {
        #-----------#
        # Chlorophycean MitochondrialCode (transl_table=16)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|((U|T)AG))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)AA)|((T|U)GA))",
            "*" => "(((U|T)AA)|((T|U)GA))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 21) {
        #-----------#
        # Trematode Mitochondrial Code (transl_table=21)
        #-----------#
        %p2c = (
            "B" => "((A|G|R)(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R)))",
            "R" => "(CG.)",
            "S" => "(((U|T)C.)|(AG.))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y))",
            "_" => "((U|T)A(A|G|R))",
            "*" => "((U|T)A(A|G|R))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AAG)",
            "N" => "(AA(U|T|C|Y|A))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)(A|G|R))",
            "W" => "((U|T)G(A|G|R))",
            "X" => "...",
        );
    } elsif ($codontable == 22) {
        #-----------#
        # Scenedesmus obliquus mitochondrial Code (transl_table=22)
        #-----------#
        %p2c = (
            "B" => "(A(U|T)G)",
            "L" => "((C(U|T).)|((U|T)(U|T)(A|G|R))|((T|U)AG))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C(U|T|C|Y|G))|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "((U|T)(C|A|G|R)A)",
            "*" => "((U|T)(C|A|G|R)A)",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    } elsif ($codontable == 23) {
        #-----------#
        # Thraustochytrium Mitochondrial Code (transl_table=23)
        #-----------#
        %p2c = (
            "B" => "(((A|G|R)(U|T)G)|(A(U|T)(U|T)))",
            "L" => "((C(U|T).)|((U|T)(U|T)G))",
            "R" => "((CG.)|(AG(A|G|R)))",
            "S" => "(((U|T)C.)|(AG(U|T|C|Y)))",
            "A" => "(GC.)",
            "G" => "(GG.)",
            "P" => "(CC.)",
            "T" => "(AC.)",
            "V" => "(G(U|T).)",
            "I" => "(A(U|T)(U|T|C|Y|A))",
            "_" => "(((U|T)A(A|G|R))|((T|U)GA)|((T|U)(T|U)A))",
            "*" => "(((U|T)A(A|G|R))|((T|U)GA)|((T|U)(T|U)A))",
            "C" => "((U|T)G(U|T|C|Y))",
            "D" => "(GA(U|T|C|Y))",
            "E" => "(GA(A|G|R))",
            "F" => "((U|T)(U|T)(U|T|C|Y))",
            "H" => "(CA(U|T|C|Y))",
            "K" => "(AA(A|G|R))",
            "N" => "(AA(U|T|C|Y))",
            "Q" => "(CA(A|G|R))",
            "Y" => "((U|T)A(U|T|C|Y))",
            "M" => "(A(U|T)G)",
            "W" => "((U|T)GG)",
            "X" => "...",
        );
    }


    #---------------------------------------------------------------#
    # make codon sequence, $qcodon,  with all possible combinations
    #---------------------------------------------------------------#

    $peplen = length($pep);
    $qcodon = '';
    foreach $i (0..$peplen - 1) {
        $peppos = $i + 1;
        $tmpaa = substr($pep, $i, 1);
        if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) {
            if ($qcodon !~ /\w/ && substr($pep, $i, 1) eq "M") {

                #---------------#
                # the first Met
                #---------------#

                $qcodon .= $p2c{"B"};
            } else {
                $qcodon .= $p2c{substr($pep, $i, 1)};
            }
        } elsif ($tmpaa =~ /\d/) {
            $qcodon .= "." x $tmpaa;
        } elsif ($tmpaa =~ /[-\.]/) {
            # nothing to do
        } else {
            $message = "pepAlnPos $peppos: $tmpaa unknown AA type. Taken as 'X'";
            push(@{$retval{'message'}}, $message);
            $qcodon .= $p2c{'X'};
        }
    }
    # print "$qcodon\n";


    #-----------------------------#
    # does $nuc contain $qcodon ?
    #-----------------------------#

    if ($nuc =~ /$qcodon/i) {
        $codon = $&;

        $retval{'codonseq'} = $codon;
        $retval{'result'} = 1;

    } else {
        #-------------------#
        # make 10 aa anchor
        #-------------------#

#        undef(@{$retval{'message'}});

#        $modpep = $pep;
#        1 while $modpep =~ s/(.{10})(.{10,})/$1\n$2/;
#        @anchor = split(/\n/, $modpep);

        my (@preanchor) = ();
        my ($tmpanc) = undef;
        my $nanc = 0;
        foreach $i (0..$peplen - 1) {
            $tmpaa = substr($pep, $i, 1);
            $tmpanc .= $tmpaa;
            ++$nanc if ($tmpaa !~ /-/);
            if ($nanc == 10 || $i == $peplen - 1) {
                push(@preanchor, $tmpanc);
                undef($tmpanc);
                $nanc = 0;
            }
        }
        undef(@anchor);
        my $lastanchorlen = length($preanchor[-1]);
        if ($lastanchorlen < 10) {
            foreach $i (0..$#preanchor - 1) {
                if ($i < $#preanchor - 1) {
                    push(@anchor, $preanchor[$i]);
                } else {
                    push(@anchor, "$preanchor[$i]$preanchor[$i + 1]");
                }
            }
        } else {
            @anchor = @preanchor;
        }

        undef($wholecodon);
        foreach $i (0..$#anchor) {
            # print "    $anchor[$i]\n";
            $anclen = length($anchor[$i]);
            undef(@qcodon);
            undef(@fncodon);
            foreach $j (0..$anclen - 1) {
                $peppos = $i * 10 + $j + 1;
                $tmpaa = substr($anchor[$i], $j, 1);
                if ($tmpaa =~ /[ACDEFGHIKLMNPQRSTVWY_\*XU]/) {
                    if ($i == 0 && $qcodon[$i] !~ /\w/ && $tmpaa eq "M") {

                        #----------------------------------#
                        # the first Met can be AGT|GTG|CTG
                        #----------------------------------#

                        $qcodon[$i] .= "((A|C|G|R)TG)";
                    } else {
                        $qcodon[$i] .= $p2c{$tmpaa};
                    }
                    $fncodon[$i] .= $p2c{'X'};
                } elsif ($tmpaa =~ /\d/) {
                    $qcodon[$i] .= "." x $tmpaa;
                    $fncodon[$i] .= "." x $tmpaa;
                } elsif ($tmpaa =~ /[-\.]/) {
                    # nothing to do
                } else {
                    #del $message = "pepAlnPos $peppos: $tmpaa unknown AA type. Replaced by 'X'";
                    #del push(@{$retval{'message'}}, $message);
                    $qcodon[$i] .= $p2c{'X'};
                    $fncodon[$i] .= $p2c{'X'};
                }
            }
            if ($nuc =~ /$qcodon[$i]/i) {
                $wholecodon .= $qcodon[$i];
            } else {
                $wholecodon .= $fncodon[$i];
            }
        }

        if ($nuc =~ /$wholecodon/i) {
            $codon = $&;
            $codonpos = 0;
            my $tmpnaa = 0;
            foreach $i (0..$peplen - 1) {
                $peppos = $i + 1;
                $tmpaa = substr($pep, $i, 1);
                undef($tmpcodon);
                if ($tmpaa !~ /\d/ && $tmpaa !~ /-/) {
                    ++$tmpnaa;
                    $tmpcodon = substr($codon, $codonpos, 3);
                    $codonpos += 3;
                    if ($tmpnaa == 1 && $tmpaa eq "M") {
                        if ($tmpcodon !~ /((A|C|G|R)TG)/i) {
                            $message = "pepAlnPos $peppos: $tmpaa does not correspond to $tmpcodon";
                            push(@{$retval{'message'}}, $message);
                        }
                    } elsif ($tmpcodon !~ /$p2c{$tmpaa}/i) {
                        $message = "pepAlnPos $peppos: $tmpaa does not correspond to $tmpcodon";
                        push(@{$retval{'message'}}, $message);
                    }
                } elsif ($tmpaa =~ /\d/i) {
                    $tmpcodon = substr($codon, $codonpos, $tmpaa);
                    $codonpos += $tmpaa;
                }
                # print "$tmpaa    $tmpcodon\n";
            }
            #$codon;

            $retval{'codonseq'} = $codon;
            $retval{'result'} = 2;

        } else {

            $retval{'result'} = -1;

        }

    }

    return(%retval);
}


###############################################################################

sub my1while {
    my($seq, $wlen) = @_;

    my($seqlen, $tmppos, @frags);

    $seqlen = length($seq);
    $tmppos = 1;


    while ($tmppos <= $seqlen) {
        my $tmpfrag = substr($seq, $tmppos - 1, $wlen);
        push(@frags, $tmpfrag);
        $tmppos += $wlen;
    }

    return(@frags);

}
###############################################################################

sub write_lua_file{
	open LUA, ">$opt_d/number_inodes.lua" || die "$!";
	print LUA 'function start_tree() n = 0 end
function node()
if not l then n = n + 1; N.lbl = \'N\'..n end
end';
	close LUA;
}

###############################################################################

sub extract_fasta{
	# Make sure files open before wasting any time doing anything
	my $id_l = shift;
	my $fasta_f = shift;
	my $out_file = shift;
	$debug && print STDERR "extract_fasta: $id_l from  $fasta_f into $out_file\n";
	open(LIST_FH, "<$id_l") || die "Could not open $id_l\n";
	open(OUT_FH, ">$out_file") || die "Could not open $out_file\n";
	my %list_hash;
	my $list_length=0;
	my @id_order;
	while(<LIST_FH>){
		chomp;
		my ($id)=split /\t/, $_;
		$id =~ s/\'//g;
		$list_hash{$id}=1;
		$list_length++;
		push @id_order,$id;
	}
	close(LIST_FH);
	###############################################################################
	# Read in features

	my $num_found=0;

	open(FASTA_FH, "<$fasta_f") || die "Could not open $fasta_f\n";

	#print STDERR "Processing FASTA file...\n";
	my ($defline, $prev_defline, $sequence) = ('','','');
	while(<FASTA_FH>){
		chomp;
		
		if(/^>/){
			$defline=$_;
			if($sequence ne ""){
				$num_found += process_record($prev_defline, $sequence, \%list_hash);
				$sequence="";
			}
			$prev_defline=$defline;
		}else{
			$sequence.=$_;
		}
	}
	$num_found += process_record($prev_defline, $sequence, \%list_hash);

	close(FASTA_FH);
	foreach my $id (@id_order){
		print OUT_FH $list_hash{$id};
	}
	print STDERR "extract_fasta: $num_found out of $list_length sequences found\n";
	die "Some sequences were not found..." unless $num_found == $list_length;
	#print STDERR "Completed.\n";
	close(OUT_FH);
	###############################################################################


	sub process_record{
		my $defline = shift;
		my $sequence = shift;
		my $list_hash = shift;
		my $num_found = 0;
		$defline=~/^>(\S+)/;
		my $id=$1;
		
		if(defined($$list_hash{$id})){

			$$list_hash{$id} = "$defline\n";

			my $length=length($sequence);
			my $width=50;
			my $pos=0;
			do{
				my $out_width=($width>$length)?$length:$width;
				$$list_hash{$id} .= substr($sequence, $pos, $width) . "\n";
				$pos+=$width;
				$length-=$width;
			}while($length>0);
			
			$num_found++;
		}
	return $num_found;
	}
}
