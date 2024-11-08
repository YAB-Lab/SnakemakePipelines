#!/usr/bin/env sh

usage () {
    echo ; echo usage: `basename $0` -c chroms -d hmm_fit_dir -f diffac -l chromlengths -m msg_dir -n pnathresh -p chroms2plot -t thinfac
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
      d) hmm_fit_dir=$OPTARG ;; #
      f) diffac=$OPTARG ;; #
      l) chromlengths=$OPTARG ;; #
      m) msg=$OPTARG ;; #
      n) pnathresh=$OPTARG ;; #
      p) chroms2plot=$OPTARG ;;
      t) thinfac=$OPTARG ;; #

      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$hmm_fit_dir" ] && [ -n "$chroms" ] && [ -n "$msg" ] && [ -n "$thinfac" ] || usage

date
echo "version 0.0"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate r_msg

echo "Outputting genotype tables and summary plots"
cmd="Rscript $msg/summaryPlots.R -c $chroms -p $chroms2plot -d $hmm_fit_dir -t $thinfac -f $difffac -n $pnathresh -l $chromlengths"

exec 3>&1; exec 1>&2; echo $cmd; exec 1>&3 3>&-
echo $cmd
$cmd || {
    echo "Error during ancestry probs"
}
