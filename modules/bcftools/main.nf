//we index and filter PASS variants 
process BCFTOOLS_PASS_INDEX {
    tag "$meta-bcf"
    publishDir "$params.outdir/VCFsfiltered/${meta}", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
     'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_PASS.vcf.gz"), path("*_PASS.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.csi"), optional:true, emit: csi
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
   touch ${file_tag}_PASS.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.2 
    END_VERSIONS
   """
}


//we index and filter PASS variants 
process BCFTOOLS_CONCAT{
    tag "$meta-bcf"
    publishDir "$params.outdir/VCFMerge/", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
       'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
     'biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), path(vcfsnp), path(indexsnp), path(vcfindel), path(indexindel)

    output:
    tuple val(meta), path("*_pm.vcf.gz"), path("*_pm.vcf.gz.tbi"), emit: tbi
    path "versions.yml"           , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def file_tag = vcfsnp[0].name.replace("_snvs_PASS.vcf.gz","").replace(".vcf","")
   """
   bcftools concat --allow-overlaps -O z -o ${file_tag}_pm.vcf.gz ${vcfsnp} ${vcfindel}
   bcftools index -t ${file_tag}_pm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
   """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def file_tag = vcfsnp[0].name.replace("_snvs_PASS.vcf.gz","").replace(".vcf","")
   """
   echo bcftools concat --allow-overlaps -O z -o ${file_tag}_pm.vcf.gz ${vcfsnp} ${vcfindel}
   echo bcftools index -t ${file_tag}_pm.vcf.gz
   touch ${file_tag}_pm.vcf.gz
   touch ${file_tag}_pm.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.2 
    END_VERSIONS
   """
}

