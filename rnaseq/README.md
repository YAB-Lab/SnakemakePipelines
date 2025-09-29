### Instructions:

To run the RNA-Seq pipelines, edit the `config.yml` file with input genome/gtf files and a JSON file of the smaples. You can create the JSON samples file using the `make_json_sample` script like so:

```
make_json_sample.py /path/to/fastq/directory/*fq.gz
```

(`make_json_sample.py` can be found on the main repository page)


FastQ files must have the filename format:

<sample_name>.R1.fastq.gz or <sample_name>.R1.fq.gz

`<sample_name>` should not include any periods.

### Running the pipeline:

Run the pipeline from within a Conda environemnt in which Snakemake in installed (all other pipeline requirements are included as Conda envs or Snakemake wrappers). For example, if your `config.yml` file is in the current directory:

```
snakemake --use-conda -j 12
```

