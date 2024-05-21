
nextflow.enable.dsl = 2

def help_function() {
    help = """somatic_point_mutations.nf: somatic point mutation caller pipeline.
             |Required arguments:
             |  --tn  tumor normal pairs
             |                [default: ${params.tn}]
             |  --vcfs  VCFs of tumor normal pairs from streka 
             |                [default: ${params.vcfs}]
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
             |   Annovar options:
             |  --annovar_bin   path to annovar_table.pl exec 
             |              [default : ${params.annovar_bin}
             | --annovar_bd  path to annovar database hg38
             |               [default : ${params.annovar_bd}]
             | --annovar_protocol databases included in annovar analysis
             |              [default: ${params.annovar_protocol}]
             | --annovar_operation operation according to annovar selected databases
             |              [default: ${params.annovar_operation}]
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

    conda "${baseDir}/strelka_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'biocontainers/strelka:2.9.10--h9ee0642_1' }"

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
    # we add the genotype information to strelka files
    sh ${baseDir}/auxfiles/fixGT.sh *.vcf.gz 

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

//we run ANNOVAR annotation tool
process ANNOVAR{
    tag "${meta}-AV"
  publishDir "$params.outdir/annovar", mode: "copy"

  input:
    tuple val(meta), path(variants), path(vindex)
  output:
    tuple val (meta), path("*multianno.vcf"), path("*multianno.txt"), emit: annovar

  script:
   """
    ${params.annovar_bin} ${variants} \\
    ${params.annovar_bd} -out ${meta}_annovar  --thread  $task.cpus \\
    -nastring . -vcfinput --buildver hg38  --codingarg -includesnp --remove --onetranscript \\
    -protocol ${params.annovar_protocol} -operation ${params.annovar_operation} 
    """
  stub:
    """
   echo ${params.annovar_bin} ${variants} \\
        ${params.annovar_bd} -out ${meta}_annovar --thread  $task.cpus \\
    -nastring . -vcfinput --buildver hg38  --codingarg -includesnp --remove --onetranscript \\
    -protocol ${params.annovar_protocol} -operation ${params.annovar_operation} 
    touch ${meta}_annovar_multianno.vcf ${meta}_annovar_multianno.txt
    """
}


include { BCFTOOLS_PASS_INDEX as BFSNV } from './modules/bcftools'
include { BCFTOOLS_PASS_INDEX as BFINDEL } from './modules/bcftools'
include { BCFTOOLS_CONCAT as BCON} from './modules/bcftools'


workflow {
     if(params.help){
        help_function()
     }

     //we do add the vcfs options to run annovar directly from unmerged strelka files
    if( (params.tn  == null && params.vcfs  == null)  || params.fasta == null || params.fai == null ){
        help_function()
    }

    if(params.tn!=null){
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
//This part is for strelka somatic point mutation caller 
    STRELKA_SOMATIC(pairs,ch_bed_target,ch_bed_target_index,ch_ref_fasta,ch_ref_fai)
    BFSNV(STRELKA_SOMATIC.out.vcf_snvs)
    snps=BFSNV.out.tbi
    BFINDEL(STRELKA_SOMATIC.out.vcf_indels)
    indels=BFINDEL.out.tbi
    pm_merge=snps.join(indels)
    BCON(pm_merge)
    //BCON.out.tbi.view()
    ANNOVAR(BCON.out.tbi)  
    //we do call annovar to run the annotation of point mutations
  }else{
    if(params.vcfs!=null){
        vcfs=Channel.fromPath(params.vcfs) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sampleId, file(row.snpvcf), file(row.snpindex), file(row.indelvcf),file(row.indelindex))}
        BCON(vcfs)
        ANNOVAR(BCON.out.tbi)      
    }
  } 
}



