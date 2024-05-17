# somatic_point_mutations
call somatic point mutations from tumor/normal pairs

## nextflow 
```
nextflow run somatic_point_mutations/main.nf --tn test.csv -params-file strelka-params.yml -profile kutral
```

### Tumor/Normal file

CSV file indicating the paths to CRAM or BAM files, including index and alternative manta_indel VCFs files.
```
sampleId,normal,normal_index,tumor,tumor_index,manta_indel,manta_indel_index
A,A.cram,A.cram.crai,AT.cram,AT.cram.crai,,
B,B.cram,B.cram.crai,BT.cram,BT.cram.crai,,
C,C.cram,C.cram.crai,CT.cram,CT.cram.crai,,
D,D.cram,D.cram.crai,DT.cram,DT.cram.crai,DT.manta.vcf.gz,DT.manta.vcf.gz.tbi
```

### Available Options

```
somatic_point_mutations.nf: somatic point mutation caller pipeline.
Required arguments:
  --tn  tumor normal pairs
                [default: false]
  --fasta  reference file in fasta format
                [default: false]
  --fai  reference index file in fai format
                [default: false]

Optional arguments:
  --outdir       The NextFlow result directory.
                [default: ./results]
  --exome       Data is Exome rather than WGS.
                [default: false]
  --target_bed  target region for strelka in bed format hg38
                [default: /somatic_point_mutations/auxfiles/hg38.bed.gz]
  --target_bed_index  bed index for target regions
                [default: /somatic_point_mutations/auxfiles/hg38.bed.gz.tbi]
```


