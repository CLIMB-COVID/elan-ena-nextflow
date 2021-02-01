#!/usr/bin/env nextflow

if( !params.study ) error "Missing ENA `study` param"
if( !params.manifest ) error "Missing ena.csv `manifest` param"
if( !params.webin_jar ) error "Missing `webin_jar` path param"

flag_ascp = ""
if ( params.ascp ){
    flag_ascp = "-ascp"
}
flag_test = ""
if ( params.test ){
    flag_test = "-test"
}

Channel
    .fromPath(params.manifest)
    .splitCsv(header:true, sep:'\t')
    .map{ it << [climb_fn: file(it.climb_fn), hoot:0] }
    .set{manifest_ch}

process prep_fasta {
    label 'bear'

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
    tuple row, file(ena_fasta), file("${row.climb_fn.baseName}.chr_list.txt.gz") into genmanifest_ch

    script:
    """
    echo "${row.published_name} 1 Monopartite" | gzip > ${row.climb_fn.baseName}.chr_list.txt.gz
    """
}

process generate_manifest {
    input:
    tuple row, file(ena_fasta), file(chr_list) from genmanifest_ch

    output:
    tuple row, file(ena_fasta), file(chr_list), file("${row.climb_fn.baseName}.manifest.txt") into webin_validate_ch

    script:
    """
    echo "STUDY ${params.study}
SAMPLE ${row.ena_sample_id}
RUN_REF ${row.ena_run_id}
ASSEMBLYNAME ${row.assemblyname}
DESCRIPTION ${row.published_name}
ASSEMBLY_TYPE COVID-19 outbreak
MOLECULETYPE genomic RNA
COVERAGE ${row.mean_cov}
PROGRAM ${row.program}
PLATFORM ${row.platform}
CHROMOSOME_LIST ${chr_list}
FASTA ${ena_fasta}
AUTHORS ${row.authors}
ADDRESS ${row.address}
SUBMISSION_TOOL ${workflow.repository}
SUBMISSION_TOOL_VERSION ${workflow.revision}@${workflow.commitId.substring(0,7)}" > ${row.climb_fn.baseName}.manifest.txt
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

process webin_submit {
    input:
    tuple row, file(ena_fasta), file(chr_list), file(ena_manifest) from webin_submit_ch

    script:
    """
    java -jar ${params.webin_jar} -context genome -userName \$WEBIN_USER -password \$WEBIN_PASS -manifest ${ena_manifest} -centerName '${row.center_name}' ${flag_ascp} -submit ${flag_test}
    """
}
