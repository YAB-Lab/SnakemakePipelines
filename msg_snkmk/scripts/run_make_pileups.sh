#!/usr/bin/env sh


usage () {
    echo ; echo ; echo usage: `basename $0` -o outdir -p parent1 -q parent2 -i indiv -S samtools_path -m msg_dir
    echo ; echo ; echo
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

# src=$(dirname $0)

while getopts "i:m:o:p:q:S:t:u:v:w:x:y:z:" opt
do
  case $opt in
      i) indiv=$OPTARG ;; #
      m) msg=$OPTARG ;; #
      o) outdir=$OPTARG ;; #
      p) parent1=$OPTARG ;; #
      q) parent2=$OPTARG ;; #
      S) samtools=$OPTARG ;; #

      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$outdir" ] && [ -n "$parent1" ] && [ -n "$parent2" ] && [ -n "$indiv" ] && [ -n "$samtools" ] || usage
[ -d $outdir ] || mkdir -p $outdir

indivdir=$outdir/$indiv
[ -d $indivdir ] || mkdir -p $indivdir

date
echo "version 0.0"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate samtools_0.1.9

echo "Creating pileup for $indiv"
echo "bash $msg/make-pileups.sh -i $indiv -d $indivdir -p $parent1 -q $parent2 -S ${samtools} 2>&1 | grep -vF 'deleted'"
bash $msg/make-pileups.sh -i $indiv -d $indivdir -p $parent1 -q $parent2 -S ${samtools} 2>&1 | grep -vF 'deleted'
