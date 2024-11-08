#!/usr/bin/env sh

usage () {
    echo usage: `basename $0` -o outdir -R Routdir -s sex -p deltaparent1 -q deltaparent2 -a recRate -r rfac -i indiv -c chroms -x sexchroms -z priors -t theta -u one_site_per_contig -j pepthresh -n msg_dir -m gff_thresh_conf -v use_filter_hmmdata_pl -l chromlengths -y chroms2plot
    exit 2
}

die () {
    echo "$1"
    exit ${2:-1}
}

src=$(dirname $0)

while getopts "a:b:c:e:f:g:h:i:j:k:l:m:n:o:p:q:r:R:s:S:t:u:v:w:x:y:z:" opt
do
  case $opt in
      a) recRate=$OPTARG ;; #
      c) chroms=$OPTARG ;; #
      p) deltapar1=$OPTARG ;; #
      q) deltapar2=$OPTARG ;; #
      i) indiv=$OPTARG ;; #
      j) pepthresh=$OPTARG ;; #
      l) chromlengths=$OPTARG ;; #
      m) gff_thresh_conf=$OPTARG ;; #
      n) msg=$OPTARG ;; #
      o) outdir=$OPTARG ;; #
      r) rfac=$OPTARG ;; #
      R) Routdir=$OPTARG ;; #
      s) sex=$OPTARG ;; #
      t) theta=$OPTARG ;; #
      u) one_site_per_contig=$OPTARG ;; #
      v) filter_hmmdata_pl=$OPTARG ;; #
      x) sexchroms=$OPTARG ;; #
      y) chroms2plot=$OPTARG ;;
      z) priors=$OPTARG ;; #

      *) usage ;;
  esac
done
shift $(($OPTIND - 1))

[ -n "$outdir" ] && [ -n "$indiv" ] || usage
[ -n "$Routdir" ] || usage
[ -d $outdir ] || mkdir -p $outdir

[ -n "$deltapar1" ] || deltapar1=.01
[ -n "$deltapar2" ] || deltapar2=$deltapar1
[ -n "$recRate" ] ||   recRate=0
[ -n "$rfac" ] ||      rfac=.000001

date
echo "version 0.0"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate r_msg

indivdir=$outdir/$indiv
[ -d $indivdir ] || mkdir -p $indivdir

echo "Fitting HMM for $indiv"
Rindivdir=$Routdir/$indiv
[ -d $Rindivdir ] || mkdir -p $Rindivdir
cmd="Rscript $msg/fit-hmm.R -d $outdir -i $indiv -s $sex -o $Routdir -p $deltapar1 -q $deltapar2 -a $recRate -r $rfac -c $chroms -x $sexchroms -z $priors -t $theta -g $gff_thresh_conf -u $one_site_per_contig -j $pepthresh -v $filter_hmmdata_pl -l $chromlengths -y $chroms2plot"

exec 3>&1; exec 1>&2; echo $cmd; exec 1>&3 3>&-
echo $cmd
$cmd || {
    echo "Error during fit-hmm.R for $indiv"
}
