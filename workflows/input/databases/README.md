# Reference Databases (Local Only)

This directory is intentionally kept out of git because reference genomes and indices are large.
Only this `README.md` is committed so the folder exists on GitHub.

## hg38 FASTA (UCSC)

From the repository root:

```bash
mkdir -p workflows/input/databases/hg38_ucsc
cd workflows/input/databases/hg38_ucsc

curl -L -o hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
```

## Bowtie2 Index For KneadData (hg38)

This repository expects the Bowtie2 index prefix at:

`workflows/input/databases/hg38_bowtie2/hg38`

Build it with:

```bash
mkdir -p workflows/input/databases/hg38_bowtie2

# Keep the gz file, but Bowtie2 needs the decompressed FASTA for indexing.
gunzip -k workflows/input/databases/hg38_bowtie2/hg38.fa.gz

bowtie2-build \
  workflows/input/databases/hg38_bowtie2/hg38.fa \
  workflows/input/databases/hg38_bowtie2/hg38
```

After indexing, you should see files like `hg38.1.bt2` (or `hg38.1.bt2l`) in `hg38_bowtie2/`.
