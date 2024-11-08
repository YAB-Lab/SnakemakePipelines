#!/usr/bin/env sh


usage () {
    echo ; echo ; echo usage: `basename $0` -s samdir -o outdir -p parent1 -q parent2 -i indiv -c chroms -w bwaalg -k repeat_threshold -S samtools_path -m msg_dir
    echo ; echo ; echo
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

# src=$(dirname $0)

while getopts "a:c:e:f:g:h:i:j:k:l:m:n:o:p:q:r:R:s:S:t:u:v:w:x:y:z:" opt
do
  case $opt in
      c) chroms=$OPTARG ;; #
      i) indiv=$OPTARG ;; #
      k) repeatthresh=$OPTARG ;; #
      m) msg=$OPTARG ;; #
      o) outdir=$OPTARG ;; #
      p) parent1=$OPTARG ;; #
      q) parent2=$OPTARG ;; #
      s) samdir=$OPTARG ;; #
      S) samtools=$OPTARG ;; #
      w) bwaalg=$OPTARG ;; #

      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$samdir" ] && [ -n "$outdir" ] && [ -n "$parent1" ] && [ -n "$parent2" ] && [ -n "$indiv" ] && [ -n "$samtools" ] || usage
[ -d $outdir ] || mkdir -p $outdir

indivdir=$outdir/$indiv
[ -d $indivdir ] || mkdir -p $indivdir


date
echo "version 0.0"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate msg

echo "Extracting reference allele information from SAM files for $indiv ($parent1 and $parent2)"
echo "python $msg/extract-ref-alleles.py -i $indiv -d $samdir -o $indivdir --parent1 $parent1 --parent2 $parent2 --chroms $chroms --bwa_alg $bwaalg --repeat_threshold $repeatthresh"
python $msg/extract-ref-alleles.py -i $indiv -d $samdir -o $indivdir --parent1 $parent1 --parent2 $parent2 --chroms $chroms --bwa_alg $bwaalg --repeat_threshold $repeatthresh || {
    echo "Error during extract-ref-alleles.py for $indiv"
}
