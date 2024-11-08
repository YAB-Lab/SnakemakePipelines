#!/usr/bin/env sh


usage () {
    echo ; echo ; echo usage: `basename $0` -o outdir -i indiv -m msg_dir -c all
    echo ; echo ; echo
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

# src=$(dirname $0)

while getopts "c:i:m:o:" opt
do
  case $opt in
      c) chroms=$OPTARG ;; #
      i) indiv=$OPTARG ;; #
      m) msg=$OPTARG ;; #
      o) outdir=$OPTARG ;; #

      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

# [ -n "$outdir" ] && [ -n "$indiv" ] && [ -n "$msg" ] && [ -n "$chrom" ] || usage
[ -d $outdir ] || mkdir -p $outdir

indivdir=$outdir/$indiv
[ -d $indivdir ] || mkdir -p $indivdir

date
echo "version 0.0"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate r_msg

echo "Writing HMM input data for $indiv"
cmd="Rscript $msg/write-hmm-data.R -i $indiv -d $indivdir -c $chroms"
echo $cmd
exec 3>&1; exec 1>&2; echo $cmd; exec 1>&3 3>&-
$cmd || {
    echo "Error during write-hmm-data.R for $indiv"
}
