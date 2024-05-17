
nextflow.enable.dsl = 2

params.reference = ""
params.subreads = ""
params.aln_only = null


def help_function() {
    help = """somatic_point_mutations.nf: somatic point mutation caller pipeline.
             |Required arguments:
             |  --tn  tumor normal pairs
             |                [default: ${params.tn}]
             |  --fasta  reference file in fasta format
             |                [default: ${params.fasta}]
             |  --fai  reference index file in fai format
             |                [default: ${params.fai}]
             |
             |Optional arguments:
             |  --outdir       The NextFlow result directory.
             |                [default: ${params.outdir}]
             |  --exome       Data is Exome rather than WGS.
             |                [default: ${params.exome}]
             |  --target_bed  target region for strelka in bed formad
             |                [default: ${params.target_bed}]
             |  --target_bed_index  bed index 
             |                [default: ${params.target_bed_index}]
             |                """.stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

/*
The stub block can be defined before or after the script block. 
When the pipeline is executed with the -stub-run option and a process’s stub is not defined, the script block is executed.

This feature makes it easier to quickly prototype the workflow logic without using the real commands. 
The developer can use it to provide a dummy script that mimics the execution of the real one in a quicker manner. 
In other words, it is a way to perform a dry-run.
*/
process STRELKA_SOMATIC {
    tag "$meta-stk"

    publishDir "$params.outdir/strelka_somatic/${meta}", mode: "copy"

    //conda "${baseDir}/strelka_env.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
     //   'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
       // 'biocontainers/strelka:2.9.10--h9ee0642_1' }"

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor),  path(manta_candidate_small_indels), path(manta_candidate_small_indels_tbi)
    path target_bed
    path target_bed_index
    path fasta
    path fai

    output:
    tuple val(meta), path("*.somatic_indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    path "versions.yml"                                 , emit: versions

    
    script:
    def args = task.ext.args ?: ''
    def exome = params.exome ? " --exome ": ""
    def prefix = task.ext.prefix ?: "${meta}"
    def options_target_bed = target_bed ? "--callRegions ${target_bed}" : ""
    def options_manta = manta_candidate_small_indels ? "--indelCandidates ${manta_candidate_small_indels}" : ""
    """

    configureStrelkaSomaticWorkflow.py \\
        --tumor $input_tumor \\
        --normal $input_normal \\
        --referenceFasta $fasta \\
        ${options_target_bed} \\
        ${options_manta} \\
        $args \\
        $exome \\
        --runDir strelka


    python strelka/runWorkflow.py -m local -j $task.cpus
    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def exome = params.exome ? " --exome ": ""
    def prefix = task.ext.prefix ?: "${meta}"
    def options_target_bed = target_bed ? "--callRegions ${target_bed}" : ""
    def options_manta = manta_candidate_small_indels ? "--indelCandidates ${manta_candidate_small_indels}" : ""
    
    """
    echo configureStrelkaSomaticWorkflow.py \\
        --tumor $input_tumor \\
        --normal $input_normal \\
        --referenceFasta $fasta \\
        ${options_target_bed} \\
        ${options_manta} \\
        $args \\
        $exome \\
        --runDir strelka

    echo "" | gzip > ${prefix}.somatic_indels.vcf.gz
    touch ${prefix}.somatic_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.somatic_snvs.vcf.gz
    touch ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: 1.2
    END_VERSIONS
    """
}

//we index and filter PASS variants 
/*
process BCFTOOLS_PASS_INDEX {
    tag "$meta-bcf"
    publishDir "$params.outdir/VCFsfiltered/${meta}", mode: "copy"

    //conda "${moduleDir}/environment.yml"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //   'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
    //  'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_PASS.vcf.gz"), path("*_PASS.vcf.gz.csi"), emit: csi
    tuple val(meta), path("*.tbi"), optional:true, emit: tbi
    path "versions.yml"           , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def file_tag = vcf[0].name.replace(".vcf.gz","").replace(".vcf","")
   """
   bcftools view -f PASS -O z ${vcf} -o ${file_tag}_PASS.vcf.gz
   bcftools index -t ${file_tag}_PASS.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
   """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def file_tag = vcf[0].name.replace(".vcf.gz","").replace(".vcf","")
   """
   echo bcftools view -f PASS -O z ${vcf} -o ${file_tag}_PASS.vcf.gz
   echo bcftools index -t ${file_tag}_PASS.vcf.gz
   touch ${file_tag}_PASS.vcf.gz
   touch ${file_tag}_PASS.vcf.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.2 
    END_VERSIONS
   """
}*/



include { BCFTOOLS_PASS_INDEX as BFSNV } from './modules/bcftools'
include { BCFTOOLS_PASS_INDEX as BFINDEL } from './modules/bcftools'

workflow {
     if(params.help){
        help_function()
     }

    if(params.tn  == null  || params.fasta == null || params.fai == null ){
      //|| params.target_bed == null || params.target_bed_index == null){
        help_function()
    
    }
     //we read the ref fasta file
      ch_ref_fasta = Channel.value(file( "${params.fasta}" ))
      ch_ref_fai = Channel.value(file( "${params.fai}"))
      //optional bed file
      ch_bed_target = Channel.value(params.target_bed ? file( "${params.target_bed}"):[])
      ch_bed_target_index = Channel.value(params.target_bed_index ? file( "${params.target_bed_index}"):[])
     //we read the long-reads 
    //we reads pairs from csv
    pairs=Channel.fromPath(params.tn) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, file(row.normal), file(row.normal_index), file(row.tumor),file(row.tumor_index), row.manta_indel ? file(row.manta_indel):[],row.manta_indel_index ? file(row.manta_indel_index):[])}

//    pairs.view()
    STRELKA_SOMATIC(pairs,ch_bed_target,ch_bed_target_index,ch_ref_fasta,ch_ref_fai)
    BFSNV(STRELKA_SOMATIC.out.vcf_snvs)
    BFINDEL(STRELKA_SOMATIC.out.vcf_indels)



}



