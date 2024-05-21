# Somatic Point Mutations

This repository provides a Nextflow pipeline for calling somatic point mutations from tumor/normal pairs using Whole Genome Sequencing (WGS) or Exome data.

## Getting Started

### Running the Pipeline

To run the pipeline, use the following command:

```sh
nextflow run digenoma-lab/somatic_point_mutations -r v1.1 --tn test.csv -params-file strelka-params.yml -profile kutral
```

### Input File: Tumor/Normal Pairs

Prepare a CSV file indicating the paths to CRAM or BAM files, including index and optional manta_indel VCF files. The CSV file should follow this format:

```csv
sampleId,normal,normal_index,tumor,tumor_index,manta_indel,manta_indel_index
A,A.cram,A.cram.crai,AT.cram,AT.cram.crai,,
B,B.cram,B.cram.crai,BT.cram,BT.cram.crai,,
C,C.cram,C.cram.crai,CT.cram,CT.cram.crai,,
D,D.cram,D.cram.crai,DT.cram,DT.cram.crai,DT.manta.vcf.gz,DT.manta.vcf.gz.tbi
```

## Pipeline Options

The `somatic_point_mutations` pipeline has several required and optional arguments.

### Required Arguments

- `--tn`: CSV file with tumor/normal pairs.
- `--fasta`: Reference genome file in FASTA format.
- `--fai`: Reference genome index file in FAI format.

### Optional Arguments

- `--outdir`: Directory for Nextflow results. Default: `./results`.
- `--exome`: Set if the data is Exome rather than WGS. Default: `false`.
- `--target_bed`: Target regions for Strelka in BED format for hg38. Default: `/somatic_point_mutations/auxfiles/hg38.bed.gz`.
- `--target_bed_index`: Index for target BED regions. Default: `/somatic_point_mutations/auxfiles/hg38.bed.gz.tbi`.

### Annovar Options

- `--annovar_bin`: Path to `annovar_table.pl` executable. Default: `/annovar/annovar/table_annovar.pl`.
- `--annovar_bd`: Path to Annovar database for hg38. Default: `/databases/annovar/hg38`.
- `--annovar_protocol`: Databases included in Annovar analysis. Default: `ensGene,clinvar_20220320,revel,dbnsfp42c,gnomad30_genome,avsnp150,icgc28`.
- `--annovar_operation`: Operations according to Annovar selected databases. Default: `g,f,f,f,f,f,f`.

## Example Usage

```sh
nextflow run digenoma-lab/somatic_point_mutations -r v1.1 \
  --tn test.csv \
  --fasta /path/to/reference.fasta \
  --fai /path/to/reference.fasta.fai \
  --outdir ./results \
  --exome true \
  --target_bed /path/to/target.bed.gz \
  --target_bed_index /path/to/target.bed.gz.tbi \
  --annovar_bin /path/to/annovar/table_annovar.pl \
  --annovar_bd /path/to/annovar/hg38 \
  --annovar_protocol ensGene,clinvar_20220320,revel,dbnsfp42c,gnomad30_genome,avsnp150,icgc28 \
  --annovar_operation g,f,f,f,f,f,f
   -profile kutral
```

## Support

If you encounter any issues or have questions, please open an issue on the [GitHub repository](https://github.com/digenoma-lab/somatic_point_mutations).


