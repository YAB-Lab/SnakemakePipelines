### Instructions:

To run the RNA-Seq pipelines, edit the `config.yml` file with input genome/gtf files and a JSON file of the smaples. You can create the JSON samples file using the `make_json_samples.py` script like so:

```make_json_samples.py /path/to/fastq/directory/*fq.gz```

FastQ files must have the filename format:

<sample_name>.R1.fastq.gz or <sample_name>.R1.fq.gz

`<sample_name>` should now include any periods.
