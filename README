# STAR task for RNA-seq

This task provides a convenience wrapper around the [STAR aligner](https://github.com/alexdobin/STAR). It aligns FASTQ reads to a given transcriptome, sorts the resulting BAMs for input into [RSEM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323), and performs [flagstat](http://www.htslib.org/doc/samtools-flagstat.html).

## Running

We encourage running this task as a Docker image, which is publicly available through GitHub packages. To pull the image, first [install Docker](https://docs.docker.com/engine/install/), then run
```
docker pull docker.pkg.github.com/krews-community/rnaseq-star-task/rnaseq-star:latest
```
The task requires a pre-built STAR index to be provided as a tarball. The STAR indexes available at ENCODE (for example [ENCFF912WCJ](https://www.encodeproject.org/files/ENCFF912WCJ/)) are compatible with this task. Alternatively, you can package a STAR index you have built into a compatible tarball by running
```
tar -cf index.tar -C /path/to/directory/with/index/files .
```
Then, to align FASTQs, simply run:
```
docker run \
    --volume /path/to/inputs:/input \
    --volume /path/to/outputs:/output \
    docker.pkg.github.com/krews-community/rnaseq-star-task/rnaseq-star:latest \
    java -jar /app/star.jar --r1 /input/r1.fastq.gz --r2 /input/r2.fastq.gz \
        --index /input/index.tar --output-directory /output
```
This will produce an output directory containing reads aligned to the genome with associated flagstat QC (`output_genome.bam`) and reads aligned to the transcriptome with associated flagstat QC (`output_anno.bam`).

### Parameters
This task supports several parameters:
|name|description|default|
|--|--|--|
|r1|FASTQ containing single-end reads or reads for pair 1|(required)|
|r2|FASTQ containing reads for pair 2 in a paired-end experiment|(none)|
|index|path to index tarball; may be gzipped or not|(required)|
|output-directory|path to output directory|(required)|
|library-id|name to assign to reads in the BAM file|(empty string)|
|output-prefix|prefix to use when naming output files|output|
|cores|number of cores available to the task|1|
|ram-gb|gigabytes of RAM available to the task|16|

## For developers

The task provides integrated unit testing, which runs a simple alignment of a small number of reads to a human mitochondrial index and checks that the outputs match expected values. To run the tests, first install Docker and Java, then clone this repo, then run
```
scripts/test.sh
```
Contributions to the code are welcome via pull request.
