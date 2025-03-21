#!/usr/bin/env python3

process MINIMAP2 {
    tag "mapping of $meta assemblies"
    
  //  publishDir "${params.output_dir}", mode:'copy'
    
    conda "${projectDir}/conda_environment/minimap2_samtools.yml"
    
    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5
    
    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("${meta}.alignment.sorted.bam"),     emit: bam_ch
    path("${meta}.alignment.sorted.bam.bai"),                  emit: bai_ch
    path("${meta}.fasta.fai"),                                 emit: fai_ch
    path "versions.yml",                                       emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    # mapping
     minimap2 -ax map-ont -t 1 $assembly $reads | samtools sort > ${prefix}.alignment.sorted.bam
     samtools faidx $assembly
     samtools index -b ${prefix}.alignment.sorted.bam
     
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process MINIMAP2_SAM {
    tag "$meta"

   // publishDir "${params.output_dir}", mode:'copy'

    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5
    
    conda "${projectDir}/conda_environment/minimap2_samtools.yml"
    
    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.sam"), emit: sam_ch
    path "versions.yml" , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 -ax map-ont -t 1 $assembly $reads  > ${meta}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process SAM_SORT_AND_INDEX {
    tag "$meta"

    publishDir "${params.output_dir}/sorted_bam", mode:'copy'

    conda "${projectDir}/conda_environment/minimap2_samtools.yml"
    
    input:
    tuple val(meta), path(sam), path(assembly)

    output:
    tuple val(meta), path("*.bam"), path("${meta}.alignment.sorted.bam.bai"), path("${meta}.fasta.fai"), emit: bam_ch
    path "versions.yml",                                                                                 emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    samtools sort $sam > ${meta}.alignment.sorted.bam
    samtools faidx $assembly
    samtools index -b ${meta}.alignment.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	 END_VERSIONS
    """
}

process LONGSHOT {
    tag "$meta"
    
    publishDir "${params.output_dir}/longshot", mode:'copy'

    
    conda "${projectDir}/conda_environment/longshot_environment.yml"
    
    input:
    tuple val(meta), path(bam), path(bai), path(fai), path(assembly)
    
    output:
    tuple val(meta), path("*.vcf"), emit: vcf_ch
    tuple val(meta), path("*.floria_vcf_header"), emit: vcf_header_ch
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    longshot \\
    	-A \\
        --bam $bam \\
        --ref $assembly \\
        --out ${meta}.vcf \\
	--strand_bias_pvalue_cutoff 0.01 \\
	-F -n > ${meta}.log
        
    gzip -c ${meta}.vcf > ${meta}.vcf.gz
    
    floria_vcf_header.py ${meta}.vcf.gz ${meta}.floria_vcf_header
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longshot: \$(echo \$(longshot --version 2>&1) | cut -f9 -d' ' )
    END_VERSIONS
    """
}


process FLORIA {
    tag "Strain resolution of $meta assemblies"
    
    
    publishDir "${params.output_dir}/floria_out", mode:'copy'
    
    conda "${projectDir}/conda_environment/floria.yml"
    
    errorStrategy { task.attempt <= 3 ? "retry" : "ignore" }
    maxRetries 5
    
    input:
    tuple val(meta), path(vcf), path(bam), path(bai), path(fai), path(assembly)
    
    output:
    tuple val(meta), path("*"), emit: all_ch
    tuple val(meta), path("${meta}/${meta}.log"), emit: log_ch
    tuple val(meta), path("${meta}"), emit: floria_out_ch 
    path("${meta}/${meta}_contig_ploidy_info.tsv"),   emit: cpi_ch 
    path "versions.yml"                   , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    # strain resolution
     
    floria $args -b $bam -v $vcf -r $assembly -t 1 -o ${meta} --output-reads --gzip-reads > ${meta}.log    
    
    mv ${meta}.log ${meta}/${meta}.log
    cp ${meta}/contig_ploidy_info.tsv ${meta}/${meta}_contig_ploidy_info.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        floria: \$( floria --version 2>&1 | cut -f2 -d' ' )
    END_VERSIONS
    """
}


process COMBINE_CONTIG_PLOIDY_INFO {
    publishDir "${params.output_dir}/floria_out", mode:'copy'
    tag { 'combine contig ploidy info files'} 
    
    
    input:
    path(contig_ploidy_info_files)
    

    output:
    path("combined_contig_ploidy_info.txt"), emit: combined_cpi_ch

    
    script:
    """ 
    CONTIG_PLOIDY_INFO_FILES=(${contig_ploidy_info_files})
    
    for index in \${!CONTIG_PLOIDY_INFO_FILES[@]}; do
    CONTIG_PLOIDY_INFO_FILE=\${CONTIG_PLOIDY_INFO_FILES[\$index]}
    
    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "\$(head -1 \${CONTIG_PLOIDY_INFO_FILE})" >> combined_contig_ploidy_info.txt
    fi
    echo "\$(awk 'FNR>=2 {print}' \${CONTIG_PLOIDY_INFO_FILE})" >> combined_contig_ploidy_info.txt
    done

    """
}


process HAPLOTAG_FLORIA_OUTPUT {
    tag "$meta"
    
    publishDir "${params.output_dir}/haplotag_out", mode:'copy'
    
    conda "${projectDir}/conda_environment/haplotag.yml"

    input:
    tuple val(meta), path(bam), path(bai), path(fai), path(floria_out_dir)
    
    output:
    tuple val(meta), path("*.bam"), emit: bam_ch
    tuple val(meta), path("*.bai"), emit: bai_ch
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    haplotag_output_dir.py -d $floria_out_dir -b $bam -o $meta
    
    """
}


process FLORIA_SAMTOOLS_FASTQ {

   publishDir "${params.output_dir}/bamtofastq_floria", mode:'copy'
   
    conda "${projectDir}/conda_environment/minimap2_samtools.yml"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta}.fastq.gz"), emit: fastq_ch
    path "versions.yml",                   emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    samtools fastq -@ 4 $bam | gzip > ${meta}.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	 END_VERSIONS
    """
}


process FLORIA_FLYE {
    publishDir "${params.output_dir}/floria_flye/${meta}_FLYE", mode:'copy'
    tag "flye on $meta"
    memory { 4.GB * task.attempt }
    
     conda "${projectDir}/conda_environment/flye.yml"

    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("*.fasta")   , emit: fasta_ch
    tuple val(meta), path("*.gfa")     , emit: gfa_ch
    tuple val(meta), path("*.gv")      , emit: gv_ch
    tuple val(meta), path("*.txt")     , emit: txt_ch
    tuple val(meta), path("*.log")     , emit: log_ch
    tuple val(meta), path("*.json")    , emit: json_ch
    path "versions.yml"                , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }
    """

    flye $mode $reads --keep-haplotypes --meta --out-dir . --threads $task.cpus $args

    mv assembly.fasta ${prefix}.fasta
    mv flye.log ${prefix}.flye.log
    mv assembly_graph.gfa ${prefix}.assembly_graph.gfa
    mv assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv params.json ${prefix}.params.json
    sed -i "s/^>/>${prefix}_/g" ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

}

process FLORIA_STRAINER {
    tag "strain count & bam extraction from $meta"
    
    
    publishDir "${params.output_dir}/floria_strainer_out/$meta", mode:'copy'
    
    conda "${projectDir}/conda_environment/floria_strainer.yml"
    
    errorStrategy { task.attempt <= 2 ? "retry" : "ignore" }
    maxRetries 5
    
    input:
    tuple val(meta), path(bam), path(bai), path(fai), path(floria_out_dir)
    
    output:
    tuple val(meta), path("*"), emit: all_ch
    tuple val(meta), path("${meta}.0.bam"), emit: first_bam_ch
    tuple val(meta), path("${meta}.1.bam"), emit: second_bam_ch
    tuple val(meta), path("${meta}.2.bam"), emit: third_bam_ch, optional:true
    tuple val(meta), path("${meta}.3.bam"), emit: fourth_bam_ch, optional:true
    tuple val(meta), path("${meta}.4.bam"), emit: fifth_bam_ch, optional:true
    tuple val(meta), path("${meta}.log"), emit: log_ch
    path "versions.yml"                   , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    # strain count and extracting corresponding bam
     
    floria-strainer --bam $bam -m split $floria_out_dir --basename ${meta} &> ${meta}.log
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        floria-strainer: \$( floria-strainer --version 2>&1 | cut -f3 -d' ' )
    END_VERSIONS
    """
}

process FLORIA_STRAINER_SAMTOOLS_FASTQ {

   publishDir "${params.output_dir}/bamtofastq_floria_strainer/${meta}", mode:'copy'
   
   conda "${projectDir}/conda_environment/minimap2_samtools.yml"

    input:
    tuple val(meta), path(bam0), path(bam1)

    output:
    tuple val(meta), path("${meta}.0.fastq.gz"), emit: first_fastq_ch
    tuple val(meta), path("${meta}.1.fastq.gz"), emit: second_fastq_ch
    path "versions.yml",                   emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    samtools fastq -@ 4 $bam0 | gzip > ${meta}.0.fastq.gz
    samtools fastq -@ 4 $bam1 | gzip > ${meta}.1.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	 END_VERSIONS
    """
}

process FLORIA_STRAINER_FLYE_READ1 {
    publishDir "${params.output_dir}/floria_strainer_flye1/${meta}_FLYE", mode:'copy'
    tag "flye on $meta"
    memory { 4.GB * task.attempt }
	
    conda "${projectDir}/conda_environment/flye.yml"
    
    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5

    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("*.fasta")   , emit: fasta_ch
    tuple val(meta), path("*.gfa")     , emit: gfa_ch
    tuple val(meta), path("*.gv")      , emit: gv_ch
    tuple val(meta), path("*.txt")     , emit: txt_ch
    tuple val(meta), path("*.log")     , emit: log_ch
    tuple val(meta), path("*.json")    , emit: json_ch
    path "versions.yml"                , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }
    """

    flye $mode $reads --keep-haplotypes --meta --out-dir . --threads $task.cpus $args

    mv assembly.fasta ${prefix}.fasta
    mv flye.log ${prefix}.flye.log
    mv assembly_graph.gfa ${prefix}.assembly_graph.gfa
    mv assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv params.json ${prefix}.params.json
    sed -i "s/^>/>${prefix}_/g" ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

}


process FLORIA_STRAINER_FLYE_READ2 {
    publishDir "${params.output_dir}/floria_strainer_flye2/${meta}_FLYE", mode:'copy'
    tag "flye on $meta"
    memory { 4.GB * task.attempt }

    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5

    conda "${projectDir}/conda_environment/flye.yml"
    
    input:
    tuple val(meta), path(reads)
    val mode

    output:
    tuple val(meta), path("*.fasta")   , emit: fasta_ch
    tuple val(meta), path("*.gfa")     , emit: gfa_ch
    tuple val(meta), path("*.gv")      , emit: gv_ch
    tuple val(meta), path("*.txt")     , emit: txt_ch
    tuple val(meta), path("*.log")     , emit: log_ch
    tuple val(meta), path("*.json")    , emit: json_ch
    path "versions.yml"                , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    def valid_mode = ["--pacbio-raw", "--pacbio-corr", "--pacbio-hifi", "--nano-raw", "--nano-corr", "--nano-hq"]
    if ( !valid_mode.contains(mode) )  { error "Unrecognised mode to run Flye. Options: ${valid_mode.join(', ')}" }
    """

    flye $mode $reads --keep-haplotypes --meta --out-dir . --threads $task.cpus $args

    mv assembly.fasta ${prefix}.fasta
    mv flye.log ${prefix}.flye.log
    mv assembly_graph.gfa ${prefix}.assembly_graph.gfa
    mv assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv params.json ${prefix}.params.json
    sed -i "s/^>/>${prefix}_/g" ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """

}
