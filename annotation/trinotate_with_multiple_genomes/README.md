The Trinotate conda package is no longer supported on conda so the package must be installed in another way. 

Following the Trinotate wiki, this can be achieved either through directly downloading the files or installing it with Singularity. 
For the uses in this Snakemake workflow to maintain reproducibility and reliability, I suggest using the Singularity method. For this pipeline, Trinotate version 4.0.2 was used following the Trinotate wiki setup instructions (https://github.com/Trinotate/Trinotate/wiki/TrinotateViaSingularity).

This pipeline involves the local downloading of additional scripts: HMMer (v3.4), RNAMMer (v1.2), SignalP (v4.1), and TMHMM (v2.0c). They will need to be initialized for use by Trinotate as described in the Trinotate wiki.
