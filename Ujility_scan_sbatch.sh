#!/bin/bash
#
#SBATCH --partition=fast             # partition
#SBATCH --ntasks=1                   # tasks, here 1
#SBATCH --nodes=1                    # nodes, here 1
#SBATCH --cpus-per-task=1		# single threaded, so 1
#SBATCH --mem=1GB			# ample RAM
#SBATCH -o outputs/slurm/%A.%a.out	# for STDOUT with job ID and ARRAY ID
#SBATCH -e outputs/slurm/%A.%a.err	# for STDERR with job ID and ARRAY ID
#SBATCH --time=1-00:00			# override at runtime if necessary
#SBATCH --array=1-6			# override at runtime if necessary
#SBATCH --job-name=Uji_v1		# override at runtime if necessary
#SBATCH --mail-type END 
#SBATCH --mail-user me@there.org 

#USAGE sbatch --job-name=Uji_v1 --array=1-50 --time=0-10:00 Ujility_scan_sbatch.sh job_name offset step degen mintm maxtm
#eg: sbatch --partition=fast --job-name=UjiLity_v2 --array=1-50 --time=0-10:00 Ujility_scan_sbatch.sh Uji_v2 14032 6 256 45 51
#Note: the directory outputs/slurm needs to have been created before submitting this batch

echo START: `date` 1>&2
based='/path/to/ujility/folder/with/data_and_scripts'
cd  ${based}
#the following software is required (eg conda or modules)
module load mafft blast primer3 emboss seqkit

# Variable defs
code=$1
offset=$2
step=$3
degen=$4
mintm=$5
maxtm=$6
jobs=$SLURM_ARRAY_TASK_ID

folder="$based/Uji_results/$code.$jobs"
mkdir -p "$folder/targets"
p3="primer3_run_params.p3"
msa_aa="MEGAPRIMER_v1_PolB_Refs_plus_Tara_plus_outgroup_OTV1_non_truncated2_prot.aln"
msa_nuc="MEGAPRIMER_v1_PolB_Refs_plus_Tara_plus_outgroup_OTV1_non_truncated2_nuc.aln"
distances="-O MEGAPRIMER_v1_PolB.nw.ingroup_phylogram.distances" # optional
tree="MEGAPRIMER_v1_PolB_2.nw"
extra_nucs="MEGAPRIMER_v1_PolB_Refs_plus_Tara_plus_outgroup_OTV1_all.fna" # sequences to test for e-PCR
weights="-a MEGAPRIMER_v1_PolB_Tara_contigs_abundances.csv" # optional
outgroup="" #"-o PHYCO_Ostreococcus_tauri_virus_1" #optional
cache="" #"-c ifpncbmwzy" #optional
tmdiff=$((maxtm-mintm))
unsuccess_seqs=20 # how many unsuccessful sequences to try before closing primer pair
max_primer_eval=500 # how many primer pairs to evaluate fro each new sequence added
width=270 # width of region to amplify
margin=80 # margins to above amplifiable region in which to design primers
pos=$(( (SLURM_ARRAY_TASK_ID-1)*step+offset ))
endpos=$((pos+width+2*margin))
rst=$((pos/3))
ren=$((endpos/3))
tot=$((margin+width))

echo Subset MSA 1>&2
srun t_coffee -other_pg seq_reformat -in $msa_aa  -action +extract_block cons $rst $ren     > $folder/targets/$code.$jobs.$pos.$endpos.faa
srun t_coffee -other_pg seq_reformat -in $msa_nuc -action +extract_block cons $pos $endpos  > $folder/targets/$code.$jobs.$pos.$endpos.fna

echo Start Ujility 1>&2
srun ./UjiLity.pl -k bottomup -d $folder -t $tree $outgroup $distances -p $folder/targets/$code.$jobs.$pos.$endpos.faa -n $folder/targets/$code.$jobs.$pos.$endpos.fna -m $p3 -z $extra_nucs -s target -y $degen -x $max_primer_eval -i 0 -j 100 -l 0 -q prealigned -h ${margin}-${tot} -w 0 -e $mintm -E $maxtm -u $unsuccess_seqs -D "+" -T $tmdiff $weights $cache > outputs/$code.$jobs.out

echo END: `date` 1>&2


