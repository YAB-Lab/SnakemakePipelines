#!/usr/bin/env sh

usage () {
    echo ; echo usage: `basename $0` -c chroms -x xchrom -f diffac -n phenofile -p filepar1 -q filepar2 -k cross -t cfg -m msg_dir -s sex
    echo ;
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

src=$(dirname $0)

while getopts "a:b:c:d:f:g:h:i:j:k:l:m:n:o:p:q:r:R:s:S:t:u:v:w:x:y:z:" opt
do
  case $opt in
      c) chroms=$OPTARG ;; #
      x) xchrom=$OPTARG ;; #
      f) diffac=$OPTARG ;; #
      k) cross=$OPTARG ;; #
      m) msg=$OPTARG ;; #
      n) phenofile=$OPTARG ;; #
      p) filepar1=$OPTARG ;;
      q) filepar2=$OPTARG ;; #
      t) cfg=$OPTARG ;; #
      s) sex=$OPTARG ;; #

      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$phenofile" ] && [ -n "$chroms" ] && [ -n "$filepar1" ] && [ -n "$filepar2" ] || usage

date
echo "version 0.0"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate msg

echo "[Common]" > $cfg
echo "filepar2 = $filepar2" >> $cfg
echo "filepar1 = $filepar1" >> $cfg
echo "phenofile = $phenofile" >> $cfg
echo "difffac = $diffac" >> $cfg
echo "chroms = $chroms" >> $cfg
echo "xchroms = $xchrom" >> $cfg
echo "cross = $cross" >> $cfg
echo "autosome_prior = 0.5" >> $cfg
echo "X_prior = 0.5" >> $cfg
echo "sex_all = $sex" >> $cfg

echo ; echo "Pullin' dem tin blood!"
echo ;
cmd="python $msg/pull_thin_tsv.py $cfg"

exec 3>&1; exec 1>&2; echo $cmd; exec 1>&3 3>&-
echo $cmd
$cmd || {
    echo "Error during pull_thin"
}
