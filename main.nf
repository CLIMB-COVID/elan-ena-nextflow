#!/usr/bin/env nextflow

if( !params.study ) error "Missing ENA `study` param"
if( !params.manifest ) error "Missing ena.csv `manifest` param"
if( !params.webin_jar ) error "Missing `webin_jar` path param"
if( !params.out ) error "Missing `out` path param"

if( !System.getenv("WEBIN_USER") ) error '$WEBIN_USER unset'
if( !System.getenv("WEBIN_PASS") ) error '$WEBIN_PASS unset'

flag_ascp = ""
if( params.ascp ){
    flag_ascp = "-ascp"
}
flag_test = ""
if( params.test ){
    flag_test = "-test"
}
description_s = ""
if( params.description ){
    description_s = params.description
}

out_name = file(params.out).name
out_dir = file(params.out).parent


workflow_repo = "CLIMB-COVID/elan-ena-nextflow"
workflow_v = workflow.manifest.version
workflow_cid = ""
if( workflow.commitId ){
    workflow_repo = workflow.repository
    workflow_v = workflow.revision
    workflow_cid = workflow.commitId.substring(0, 7)
}

Channel
    .fromPath(params.manifest)
    .splitCsv(header:true, sep:'\t')
    .map{ it << [climb_fn: file(it.climb_fn), hoot:0] }
    .set{manifest_ch}

process prep_fasta {

    input:
    val row from manifest_ch

    output:
    tuple row, file("${row.climb_fn.baseName}.ena-a.fasta.gz") into chrlist_ch

    """
    elan_rehead.py ${row.climb_fn} '${row.published_name}' | gzip > ${row.climb_fn.baseName}.ena-a.fasta.gz
    """
}

process generate_chrlist {

    input:
    tuple row, file(ena_fasta) from chrlist_ch

    output:
    tuple row, file(ena_fasta), file("${row.climb_fn.baseName}.chr_list.txt.gz") into pyena_input_ch

    script:
    """
    echo "${row.published_name} 1 Monopartite" | gzip > ${row.climb_fn.baseName}.chr_list.txt.gz
    """
}

process pyena_submission {
    errorStrategy 'ignore'
    conda "environments/pyena.yaml"

    input:
    tuple row, file(ena_fasta), file(chr_list) from pyena_input_ch

    output:
    tuple row, file(ena_fasta), file(chr_list) into genmanifest_ch
    file("${row.central_sample_id}.${row.run_name}.pyena.txt") into dh_ocarina_report_ch


    script:
    """
    pyena --study-accession ${params.study} --sample-only --no-ftp \
          --sample-name COG-UK/${row.central_sample_id} \
          --sample-center-name "${row.center_name}" \
          --sample-taxon '2697049' \
          --sample-attr 'collector name' 'not provided' \
          --sample-attr 'collecting institution' "${row.center_name}" \
          --sample-attr 'collection date' ${row.collection_date} \
          --sample-attr 'geographic location (country and/or sea)' 'United Kingdom' \
          --sample-attr 'geographic location (region and locality)' '${row.adm1}' \
          --sample-attr 'definition for seropositive sample' 'not provided' \
          --sample-attr 'serotype (required for a seropositive sample)' 'not provided' \
          --sample-attr 'host common name' 'not provided' \
          --sample-attr 'host health state' 'not provided' \
          --sample-attr 'host scientific name' 'Human' \
          --sample-attr 'host sex' 'not provided' \
          --sample-attr 'host subject id' 'not provided' \
          --sample-attr 'isolate' 'not provided' \
          --sample-attr 'receipt date' '${row.received_date}' \
          --sample-attr 'sample capture status' 'active surveillance in response to outbreak' \
          --sample-attr 'virus identifier' 'not provided' \
          --sample-attr 'ENA-CHECKLIST' 'ERC000033' \
          --sample-attr 'min_cycle_threshold' '${row.min_ct}' \
          --sample-attr 'max_cycle_threshold' '${row.max_ct}' \
          --experiment-attr 'artic_primer_version' '${row.exp_primers}' \
          --experiment-attr 'artic_protocol_version' '${row.exp_protocol}' \
          --run-name ${row.published_name} > ${row.central_sample_id}.${row.run_name}.pyena.txt
    """
}

// dh_ocarina_report_ch
//     .splitCsv(header:['success', 'real', 'ena_sample_name', 'ena_run_name', 'bam', 'study_acc', 'sample_acc', 'exp_acc', 'run_acc'], sep:' ')
//     .map { row-> tuple(row.ena_run_name, row.sample_acc, row.run_acc) }
//     .set { dh_ocarina_report_ch_split }

// process tag_ocarina {
//     tag { bam }
//     label 'ocarina'
//     conda "../environments/ocarina.yaml"

//     input:
//     tuple ena_run_name, sample_acc, run_acc from dh_ocarina_report_ch_split

//     errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }
//     maxRetries 3

//     cpus 6 //# massively over-request local cores to prevent sending too much to API at once

//     script:
//     """
//     ocarina --oauth --env put publish --publish-group '${ena_run_name}' --service 'ENA-SAMPLE' --accession ${sample_acc} --public --submitted
//     ocarina --oauth --env put publish --publish-group '${ena_run_name}' --service 'ENA-RUN' --accession ${run_acc} --public --submitted
//     """
// }

process generate_manifest {
    input:
    tuple row, file(ena_fasta), file(chr_list) from genmanifest_ch

    output:
    tuple row, file(ena_fasta), file(chr_list), file("${row.climb_fn.baseName}.manifest.txt") into webin_validate_ch

    script:
    def engine = new groovy.text.SimpleTemplateEngine()
    this_description = engine.createTemplate(description_s).make(['row':row]).toString()
    """
    echo "STUDY ${params.study}
    SAMPLE ${row.ena_sample_id}
    RUN_REF ${row.pag_name}
    ASSEMBLYNAME ${row.assemblyname}
    DESCRIPTION """ << this_description << """
    ASSEMBLY_TYPE COVID-19 outbreak
    MOLECULETYPE genomic RNA
    COVERAGE ${row.mean_cov}
    PROGRAM ${row.program}
    PLATFORM ${row.platform}
    CHROMOSOME_LIST ${chr_list}
    FASTA ${ena_fasta}
    AUTHORS ${row.authors}
    ADDRESS ${row.address}
    SUBMISSION_TOOL ${workflow_repo}
    SUBMISSION_TOOL_VERSION ${workflow_v}@${workflow_cid}" > ${row.climb_fn.baseName}.manifest.txt
    """
}

process webin_validate {
    input:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) from webin_validate_ch

    errorStrategy 'ignore' //# Drop assemblies that fail to validate

    output:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) into webin_submit_ch

    script:
    """
    java -jar ${params.webin_jar} -context genome -userName \$WEBIN_USER -password \$WEBIN_PASS -manifest ${ena_manifest} -centerName '${row.center_name}' ${flag_ascp} -validate
    """
}

// process webin_submit {
//     errorStrategy 'ignore' //# Drop assemblies that fail to validate

//     input:
//     tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) from webin_submit_ch

//     output:
//     tuple row, file(ena_fasta), file(chr_list), file(ena_manifest), file("genome/${row.assemblyname.replaceAll('#', '_')}/submit/receipt.xml") into webin_parse_ch

//     script:
//     """
//     java -jar ${params.webin_jar} -context genome -userName \$WEBIN_USER -password \$WEBIN_PASS -manifest ${ena_manifest} -centerName '${row.center_name}' ${flag_ascp} -submit ${flag_test}
//     """
// }

// process receipt_parser {
//     conda "$baseDir/environments/receipt.yaml"

//     input:
//     tuple row, file(ena_fasta), file(chr_list), file(ena_manifest), file(ena_receipt) from webin_parse_ch

//     output:
//     file("${row.climb_fn.baseName}.accession.txt") into accession_report_ch

//     script:
//     """
//     parse_receipt.py ${ena_manifest} ${ena_receipt} ${row.published_name} > ${row.climb_fn.baseName}.accession.txt
//     """
// }

// accession_report_ch
//     .collectFile(keepHeader: true, name: "${out_name}", storeDir: "${out_dir}")
