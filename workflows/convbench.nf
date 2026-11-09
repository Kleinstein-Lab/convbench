/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS                                              } from '../modules/local/preprocess/main'
include { CDR3_SIMILARITY                                         } from '../modules/local/cdr3_similarity/main'
include { CDR3_SIMILARITY as CDR3_SIMILARITY_FIRSTV               } from '../modules/local/cdr3_similarity/main'
include { MILO                                                    } from '../modules/local/milo/main'
include { BCRDIST                                                 } from '../modules/local/bcrdist/main'
include { EVALUATE                                                } from '../modules/local/evaluate/main'
include { MULTIQC                                                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                  } from '../subworkflows/local/utils_nfcore_convbench_pipeline'

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

    ch_samplesheet.dump(tag: "samplesheet") // Debug

    // get CDR3 similarity process-specific similarity parameters
    def linkage = params.linkage.toString().split(',').collect { it.trim() }
    ch_linkage = channel.from(linkage)
                        .dump(tag: 'linkage') // Debug

    def threshold = params.threshold.toString().split(',').collect { it.trim().toFloat() }
    ch_threshold = channel.from(threshold)
                          .map { it -> it as Float } // convert to Float
                          .dump(tag: 'threshold') // Debug
    
    // combine clustering parameters
    ch_clustering_params = ch_linkage.combine(ch_threshold)
                                     .dump(tag: 'cdr3_similarity_clustering_params')

    //
    // MODULE: Run preprocessing (library sizes for all runs, ASC splitting in asc_mode)
    //

    if (params.asc_mode && !params.asc_guide){
        error "When --asc_mode is enabled, you must provide --asc_guide."
    }

    // an empty list skips ASC splitting
    def asc_guide = params.asc_mode ? channel.value(file(params.asc_guide, checkIfExists: true)) : channel.value([])

    PREPROCESS(ch_samplesheet, asc_guide)

    // channel: [ dataset_id, library_sizes ]
    ch_library_sizes = PREPROCESS.out.library_sizes
        .map { meta, library_sizes -> [ meta.id, library_sizes ] }

    // tool inputs, one per dataset in standard mode and one per ASC in asc_mode
    // meta.group_id is the ASC (or 'all') and meta.n_groups the number of groups in the dataset
    // channel: [ meta, airr, embedding ]
    if (params.asc_mode){

        ch_tool_inputs = PREPROCESS.out.asc_files
            .flatMap { meta, md_files, emb_files ->

                def md_list = [md_files].flatten()
                def emb_list = [emb_files].flatten()
                def group_ids = md_list.collect { md -> md.name.replaceAll(/_md\.tsv\.gz$/, '') }

                group_ids.collect { group_id ->
                    [
                        meta + [sample_id: meta.id, group_id: group_id, n_groups: group_ids.size()],
                        md_list.find { md -> md.name == "${group_id}_md.tsv.gz" },
                        emb_list.find { emb -> emb.name == "${group_id}_emb.tsv.gz" }
                    ]
                }
            }
    } else{

        ch_tool_inputs = ch_samplesheet
            .map { meta, airr, embedding ->
                [ meta + [sample_id: meta.id, group_id: 'all', n_groups: 1], airr, embedding ]
            }
    }

    ch_tool_inputs.dump(tag: 'tool_inputs') // Debug

    // CDR3 similarity inputs, one per group and clustering parameter set
    ch_cdr3_inputs = ch_tool_inputs
        .combine(ch_clustering_params)
        .map { meta, airr, embedding, linkage_method, cluster_threshold ->
            def fmeta = meta + [
                id: "${meta.id}_${linkage_method}_t${cluster_threshold}".toString(),
                linkage: linkage_method,
                threshold: cluster_threshold
            ]
            [ fmeta, airr, embedding ] // channel: [ [meta, airr, embedding] ]
        }
        .dump(tag: 'cdr3_input_final')

    def conv_tools = params.conv_tools ? params.conv_tools.split(',').collect { tool -> tool.trim() } : []

    def ch_evaluation_inputs = channel.empty()

    //
    // MODULE: Run CDR3_similarity
    //
    if ('cdr3_similarity' in conv_tools){
        CDR3_SIMILARITY(ch_cdr3_inputs)

        ch_evaluation_inputs = ch_evaluation_inputs.mix(
            CDR3_SIMILARITY.out.eval_input.map { row -> labelEvalInput(row, 'cdr3_similarity', 'cluster', 'CDR3 Similarity') }
        )
    }

    //
    // MODULE: Run CDR3_similarity_firstv
    //
    if ('cdr3_similarity_firstv' in conv_tools){
        CDR3_SIMILARITY_FIRSTV(ch_cdr3_inputs)

        ch_evaluation_inputs = ch_evaluation_inputs.mix(
            CDR3_SIMILARITY_FIRSTV.out.eval_input.map { row -> labelEvalInput(row, 'cdr3_similarity_firstv', 'cluster', 'CDR3 Similarity - First V') }
        )
    }

    //
    // MODULE: Run Milo
    //
    if ('milo' in conv_tools){
        MILO(ch_tool_inputs)

        ch_evaluation_inputs = ch_evaluation_inputs.mix(
            MILO.out.eval_input.map { row -> labelEvalInput(row, 'Milo', 'milo', 'Milo') }
        )
    }

    //
    // MODULE: Run BCRdist
    //
    if ('bcrdist' in conv_tools){
        BCRDIST(ch_tool_inputs)

        ch_evaluation_inputs = ch_evaluation_inputs.mix(
            BCRDIST.out.eval_input.map { row -> labelEvalInput(row, 'BCRdist', 'cluster', 'BCRdist') }
        )
    }

    //
    // MODULE: Run evaluation
    //

    // gather the outputs from all groups (ASCs) of each tool run, then add the library sizes
    // channel: [ meta, counts, memberships, tool_info, milo_objects, library_sizes ]
    ch_evaluate = ch_evaluation_inputs
        .map { meta, counts, memberships, tool_info, milo_object ->
            [ groupKey([meta.tool_dir, meta.id], meta.n_groups), meta, counts, memberships, tool_info, milo_object ]
        }
        .groupTuple()
        .map { _key, metas, counts, memberships, tool_info, milo_objects ->
            def meta = metas[0].findAll { k, _v -> k != 'group_id' }
            [ meta.sample_id, meta, counts, memberships, tool_info, milo_objects.flatten() ]
        }
        .combine(ch_library_sizes, by: 0)
        .map { _sample_id, meta, counts, memberships, tool_info, milo_objects, library_sizes ->
            [ meta, counts, memberships, tool_info, milo_objects, library_sizes ]
        }
        .dump(tag: 'evaluate_input')

    EVALUATE(ch_evaluate)

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
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Label a tool's evaluation inputs with the tool information used by EVALUATE
// row: [ meta, counts, memberships, tool_info, (milo_object) ]
// Tools other than Milo have no Milo object, so it is replaced with an empty list
//
def labelEvalInput(row, tool_dir, tool_type, tool_name) {
    def meta = row[0] + [tool_dir: tool_dir, tool_type: tool_type, tool_name: tool_name]
    def milo_object = row.size() > 4 ? row[4] : []
    return [ meta, row[1], row[2], row[3], milo_object ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
