#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()

// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)

final_params = check_params(merged_params)

// starting pipeline
pipeline_start_message(version, final_params)



// include processes
include { MINIMAP2_SAM; SAM_SORT_AND_INDEX; LONGSHOT; FLORIA; COMBINE_CONTIG_PLOIDY_INFO; HAPLOTAG_FLORIA_OUTPUT; FLORIA_SAMTOOLS_FASTQ; FLORIA_FLYE; FLORIA_STRAINER; FLORIA_STRAINER_SAMTOOLS_FASTQ; FLORIA_STRAINER_FLYE_READ1; FLORIA_STRAINER_FLYE_READ2 } from './modules/processes.nf' addParams(final_params)



workflow  {
          read_ch = channel
                          .fromPath( final_params.reads )
                          .map { file -> tuple(file.simpleName, file) }
			  .ifEmpty { error "Cannot find any reads matching: ${final_params.reads}" }

          assemblies_ch = channel
                                .fromPath( final_params.assemblies, checkIfExists: true )
                                .map { file -> tuple(file.simpleName, file) }
				.ifEmpty { error "Cannot find any assemblies matching: ${final_params.assemblies}" }

         joined_ch = read_ch.join(assemblies_ch)

 
	 MINIMAP2_SAM ( joined_ch )
	 
	 joined_sam_ch = MINIMAP2_SAM.out.sam_ch.join(assemblies_ch)

         SAM_SORT_AND_INDEX(joined_sam_ch)
	 
	 joined_sam_assembly_ch = SAM_SORT_AND_INDEX.out.bam_ch.join(assemblies_ch)


         LONGSHOT(joined_sam_assembly_ch)
	 
	 joined_longshot_ch = LONGSHOT.out.vcf_header_ch.join(joined_sam_assembly_ch)
	 
	 FLORIA(joined_longshot_ch)
	 
	 collected_contig_ploidy_ch = FLORIA.out.cpi_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_CONTIG_PLOIDY_INFO(collected_contig_ploidy_ch)
	 
	 HAPLOTAG_FLORIA_OUTPUT(SAM_SORT_AND_INDEX.out.bam_ch.join(FLORIA.out.floria_out_ch))
	 
	 FLORIA_SAMTOOLS_FASTQ(HAPLOTAG_FLORIA_OUTPUT.out.bam_ch)
	 
	 FLORIA_FLYE(FLORIA_SAMTOOLS_FASTQ.out.fastq_ch, params.valid_mode)
	 
	 FLORIA_STRAINER(SAM_SORT_AND_INDEX.out.bam_ch.join(FLORIA.out.floria_out_ch))
	 
	 FLORIA_STRAINER_SAMTOOLS_FASTQ(FLORIA_STRAINER.out.first_bam_ch.join(FLORIA_STRAINER.out.second_bam_ch))
	 
	 FLORIA_STRAINER_FLYE_READ1(FLORIA_STRAINER_SAMTOOLS_FASTQ.out.first_fastq_ch, params.valid_mode)
	 
	 FLORIA_STRAINER_FLYE_READ2(FLORIA_STRAINER_SAMTOOLS_FASTQ.out.second_fastq_ch, params.valid_mode)

}

workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
