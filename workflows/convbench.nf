/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SPLIT_BY_ASC           } from '../modules/local/split_by_asc/main'
include { CDR3_SIMILARITY        } from '../modules/local/cdr3_similarity/main'
include { CDR3_SIMILARITY_ASC    } from '../modules/local/cdr3_similarity/main'
include { MILO                   } from '../modules/local/milo/main'
include { MILO_ASC               } from '../modules/local/milo/main'
include { DASEQ                  } from '../modules/local/daseq/main'
include { DASEQ_ASC              } from '../modules/local/daseq/main'
include { BCRDIST                } from '../modules/local/bcrdist/main'
include { BCRDIST_ASC            } from '../modules/local/bcrdist/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { GET_ASC_AUROC          } from '../modules/local/get_asc_auroc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_convbench_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CONVBENCH {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    ch_samplesheet.dump(tag: "samplesheet")

    //
    // MODULE: Run split_by_ASC
    //

    if (params.asc_mode && !params.asc_guide){
        error "When --asc_mode is enabled, you must provide --asc_guide."
    }

    if (params.asc_mode){

        // make channel for guide
        asc_guide = Channel.value(file(params.asc_guide))

        // make ASC level channel    
        asc_splitting = SPLIT_BY_ASC(ch_samplesheet, asc_guide)

        ch_file_pairs = asc_splitting.flatMap{ meta, md_files, emb_files ->
            
            def meta_id = meta.id

            def labeled_md_files = md_files
                .collect{ it ->
                    [meta_id, it.name.split('_')[0], it]
                }

            def labeled_emb_files = emb_files
                .collect{ it ->
                    [meta_id, it.name.split('_')[0], it]
                }
            
            return tuple(
                labeled_md_files + labeled_emb_files
            )
        }
        .groupTuple(by: [0, 1])
        .map{ meta_id, asc_id, files -> // ensure files are always in the right order
        
            def airr = files.find { it.name.endsWith('_md.tsv.gz') }
            def embedding = files.find { it.name.endsWith('_emb.tsv.gz')}

            tuple(meta_id, asc_id, airr, embedding)
        }
    
    }

    //
    // MODULE: Run CDR3_similarity
    //

    if (params.conv_tools && params.conv_tools.split(',').contains('cdr3_similarity')){
        if(params.asc_mode){

            cdr3_sim_asc_result = CDR3_SIMILARITY_ASC(ch_file_pairs)

            def tool_id = 'cdr3_similarity'

            cdr3_sim_asc_input = cdr3_sim_asc_result.auc_input
                .map{meta_id, seq_file ->
                    tuple(tool_id, meta_id, seq_file)
                }
                .groupTuple(by: [0, 1])

            GET_ASC_AUROC(cdr3_sim_asc_input)
        
        } else{
            CDR3_SIMILARITY(
                ch_samplesheet
            )
        }
    }

    //
    // MODULE: Run Milo
    //
    if (params.conv_tools && params.conv_tools.split(',').contains('milo')){
        if(params.asc_mode){
            
            milo_asc_result = MILO_ASC(ch_file_pairs)

            def tool_id = 'Milo'

            milo_asc_input = milo_asc_result.auc_input
                .map{meta_id, seq_file ->
                    tuple(tool_id, meta_id, seq_file)
                }
                .groupTuple(by: [0, 1])

            GET_ASC_AUROC(milo_asc_input)
            
        } else{
            MILO(
                ch_samplesheet
            )
        }
    }
    //
    // MODULE: Run DA-seq
    //
    if (params.conv_tools && params.conv_tools.split(',').contains('daseq')){
        if(params.asc_mode){

            daseq_asc_result = DASEQ_ASC(ch_file_pairs)

            def tool_id = 'DA-seq'

            daseq_asc_input = daseq_asc_result.auc_input
                .map{meta_id, seq_file ->
                    tuple(tool_id, meta_id, seq_file)
                }
                .groupTuple(by: [0, 1])

            GET_ASC_AUROC(daseq_asc_input)

        } else{
            DASEQ(
                ch_samplesheet
            )
        }
    }

    //
    // MODULE: Run BCRdist
    //
    if (params.conv_tools && params.conv_tools.split(',').contains('bcrdist')){
        if(params.asc_mode){

            bcrdist_asc_result = BCRDIST_ASC(ch_file_pairs)

            def tool_id = 'BCRdist'

            bcrdist_asc_input = bcrdist_asc_result.auc_input
                .map{meta_id, seq_file ->
                    tuple(tool_id, meta_id, seq_file)
                }
                .groupTuple(by: [0, 1])

            GET_ASC_AUROC(bcrdist_asc_input)

        } else{
            BCRDIST(
                ch_samplesheet
            )
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'convbench_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'convbench'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
