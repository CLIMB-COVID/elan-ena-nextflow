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

    maxForks 4

    input:
    val row from manifest_ch

    output:
    tuple row, file("${row.climb_fn.baseName}.ena-a.fasta.gz") into chrlist_ch

    """
    elan_rehead.py ${row.climb_fn} '${row.assemblyname}' | gzip > ${row.climb_fn.baseName}.ena-a.fasta.gz
    """
}

process generate_chrlist {

    maxForks 4

    input:
    tuple row, file(ena_fasta) from chrlist_ch

    output:
    tuple row, file(ena_fasta), file("${row.climb_fn.baseName}.chr_list.txt.gz") into pyena_input_ch

    script:
    """
    echo "${row.assemblyname} 1 Monopartite" | gzip > ${row.climb_fn.baseName}.chr_list.txt.gz
    """
}

process pyena_submission {
    errorStrategy 'ignore'
    conda "${workflow.projectDir}/environments/pyena.yaml"

    maxForks 4

    input:
    tuple row, file(ena_fasta), file(chr_list) from pyena_input_ch

    output:
    tuple row, file(ena_fasta), file(chr_list) into genmanifest_ch
    file("${row.central_sample_id}.pyena.txt") into dh_ocarina_report_ch

    script:

    if (! params.test) {
        pyena_test_flag = "--my-data-is-ready"
    } else {
        pyena_test_flag = ""
    }

    """
    pyena --study-accession ${params.study} --no-ftp --sample-only ${pyena_test_flag}\
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
          --run-name ${row.published_name} > ${row.central_sample_id}.pyena.txt
    """
}

//   --experiment-attr 'artic_primer_version' '${row.library_primers}' \
//   --experiment-attr 'artic_protocol_version' '${row.library_protocol}' \
//   --run-name ${row.published_name} \
//   --run-file-path ${}
//   --run-center-name "${row.run_center_name}" \
//   --run-instrument '${row.instrument_model}' \
//   --run-lib-protocol '${row.library_seq_kit}|${row.library_seq_protocol}' \
//   --run-lib-source ${row.library_source} \
//   --run-lib-selection ${row.library_selection} \
//   --run-lib-strategy ${row.library_strategy} 

dh_ocarina_report_ch
    .splitCsv(header:['success', 'real', 'ena_sample_name', 'ena_run_name', 'bam', 'study_acc', 'sample_acc', 'exp_acc', 'run_acc'], sep:' ')
    .map { row-> tuple(row.ena_run_name, row.sample_acc, row.run_acc) }
    .into { dh_ocarina_report_ch_split; dh_accession_report_ch }

process tag_ocarina {
    tag { bam }
    label 'ocarina'
    conda "${workflow.projectDir}/environments/ocarina.yaml"

    maxForks 4 //# Limit number of forks to prevent sending too much to API at once
    maxRetries 3

    input:
    tuple ena_run_name, sample_acc, run_acc from dh_ocarina_report_ch_split

    errorStrategy { sleep(Math.pow(2, task.attempt) * 300 as long); return 'retry' }

    script:
    if (params.test) {
        """
        ocarina --oauth --profile test-service-outbound put publish --publish-group '${ena_run_name}' --service 'ENA-SAMPLE' --accession ${sample_acc} --public --submitted
        """        
    } else {
        """
        ocarina --oauth --profile service-outbound put publish --publish-group '${ena_run_name}' --service 'ENA-SAMPLE' --accession ${sample_acc} --public --submitted
        """
    }

}

process generate_manifest {

    maxForks 4

    input:
    tuple row, file(ena_fasta), file(chr_list) from genmanifest_ch
    tuple ena_run_name, sample_acc, run_acc from dh_accession_report_ch

    output:
    tuple row, file(ena_fasta), file(chr_list), file("${row.climb_fn.baseName}.manifest.txt") into webin_validate_ch

    script:
    def engine = new groovy.text.SimpleTemplateEngine()
    this_description = engine.createTemplate(description_s).make(['row':row]).toString()
    if (run_acc != "None"){
        """
        echo "STUDY ${params.study}
        SAMPLE ${sample_acc}
        RUN_REF ${run_acc}
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
    } else {
        """
        echo "STUDY ${params.study}
        SAMPLE ${sample_acc}
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
    
}

process webin_validate {
    input:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) from webin_validate_ch

    errorStrategy 'ignore' //# Drop assemblies that fail to validate

    output:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) into webin_submit_ch

    maxForks 16

    script:
    """
    java -jar ${params.webin_jar} -context genome -userName \$WEBIN_USER -password \$WEBIN_PASS -manifest ${ena_manifest} -centerName "${row.center_name}" ${flag_ascp} -validate ${flag_test}
    """
}

process webin_submit {
    errorStrategy 'ignore' //# Allow failed submissions to continue (This is usually due to them already having been uploaded previously)

    maxForks 16

    input:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) from webin_submit_ch

    output:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest), file("genome/${row.assemblyname.replaceAll('#', '_')}/submit/receipt.xml") into webin_parse_ch

    script:
    """
    trap 'if [ -e "genome/${row.assemblyname.replaceAll('#', '_')}/submit/receipt.xml" ]; then exit 0; else exit 1 ; fi' EXIT
    java -jar ${params.webin_jar} -context genome -userName \$WEBIN_USER -password \$WEBIN_PASS -manifest ${ena_manifest} -centerName "${row.center_name}" ${flag_ascp} -submit ${flag_test}
    """
}

process webin_parse_majora_submit {
    conda "$baseDir/environments/receipt.yaml"

    errorStrategy 'ignore'

    maxForks 4

    input:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest), file(ena_receipt) from webin_parse_ch

    output:
    file("${row.climb_fn.baseName}.accession.txt") into accession_report_ch

    script:

    if (params.test) {
        test_flag = "--test"
    } else {
        test_flag = ""
    }
    // Annoyingly I can't get ocarina to stfu about the response and it prints messages to both stdout and stderr so I'm doing this disgusting hack
    """
    parse_receipt.py ${test_flag} ${ena_manifest} ${ena_receipt} ${row.published_name} > ${row.climb_fn.baseName}.tmp.txt

    tail -n 2 ${row.climb_fn.baseName}.tmp.txt > ${row.climb_fn.baseName}.accession.txt
    """
}

accession_report_ch
    .collectFile(keepHeader: true, name: "${out_name}", storeDir: "${out_dir}")
