#!/usr/bin/env nextflow

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check the hostnames against configured profiles

def summary = [:]
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Reads']            = params.reads
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}


def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'SCAPE-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'SCAPE Workflow Summary'
    section_href: '"https://confluence.phibred.com/display/OMICS/Expression+Data+Resource+White+Paper'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}



/*
 * Create a channel for input read files
 */

if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_remove_contamination }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_remove_contamination }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_remove_contamination }
}


// Define regular variables so that they can be overwritten

if (params.singleEnd){
   read_type = "SE"
} else {
  read_type = "PE"
}

if(params.tx2gene_file) {
   my_tx2gene_file = file("${params.tx2gene_file}")
   if( !my_tx2gene_file.exists() ) exit 1, "tx2gene File not found: ${params.tx2gene_file}"
}

contaminant_removal_tool = params.contaminant_software
adapter_removal_tool = params.adapter_software
// my_salmon_index_channel = Channel.fromPath(file("${params.salmon_index_path}"))
my_quant_method = params.quant_method


/*
* Check if the required files are actually there
*/

my_salmon_index_file = file("${params.salmon_index_path}/sa.bin") // check if this file is available

if( !my_salmon_index_file.exists() ) exit 1, "Salmon Index File not found at: ${params.salmon_index_path}"

if(contaminant_removal_tool == 'bbmap') {
     my_contaminant_fasta = file("${params.contaminant_fasta}")
     if( !my_contaminant_fasta.exists() ) exit 1, "Contaminant Fasta File not found: ${params.contaminant_fasta}"
}

/*
 * STEP 1 - FastQC
*/
process fastqc_before {
    label 'low_cm'
    echo true
    tag "$name"
    //publishDir "${params.outdir}/pre_fastqc",saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into  preqc_results

    script:
    """
    fastqc -t 8 -q $reads
    """
}


/*
 * STEP 2 - Remove Contaminants
*/
process remove_contaminants {
    label 'low_cm'
    echo true
    validExitStatus 0,1
    tag "$name"
    cpus 8
    time '16h'
    memory '30 GB'
    //publishDir "${params.outdir}/qc/bbmap/stats", pattern: "*.txt"
    //publishDir "${params.outdir}/clean_fastq", pattern: "*.clean.fastq.gz"

    input:
    set val(name), file(reads) from raw_reads_remove_contamination
    file my_contaminant_fasta

    output:
    set val(name), file ("*.clean.fastq.gz") into clean_reads
    file "*.txt" into contaminant_stats

    script:
    if (params.singleEnd){
       if(contaminant_removal_tool == 'bbmap') {

       """
	        bbmap.sh  \\
                  in=${name}.fastq.gz \\
                  outu=${name}.clean.fastq.gz \\
                  ref=${my_contaminant_fasta} \\
                  maxindel=1 minid=0.95 nodisk \\
                  statsfile=${name}.trimstats.txt
       """

       }
    }else { // we have paired end now
       if(contaminant_removal_tool == 'bbmap') {

       """
        bbmap.sh  \\
                  in1=${name}_1.fastq.gz  in2=${name}_2.fastq.gz \\
                  outu1=${name}_1.clean.fastq.gz  outu2=${name}_2.clean.fastq.gz \\
                  ref=${my_contaminant_fasta} \\
                  maxindel=1 minid=0.95 nodisk \\
                  statsfile=${name}.trimstats.txt
       """
       }
    }
} // end of Step 2 - Process remove_contaminants - channel(s) output clean_reads

/*
 * STEP 3 - Remove Adapters
*/
process remove_adapters {
    label 'low_cm'
    echo true
    validExitStatus 0,1
    tag "$name"
    cpus 8
    time '16h'
    memory '30 GB'
    if(adapter_removal_tool == 'fastp') {
        publishDir "${params.outdir}/qc/adapter_removal", mode: 'copy', pattern: "*.{json,html}"
    }
    //publishDir "${params.outdir}/trimmed_clean_fastq", pattern: "*_trimmed_clean.fastq.gz"

    input:
    set val(name), file(clean_read) from clean_reads

    output:
    set val(name), file ("*_trimmed_clean.fastq.gz") into analysis_ready_reads, post_trimming_qc
    file "*.{json,html}" into adapter_removal_stats

    script:
    if (params.singleEnd){
       if(adapter_removal_tool == 'fastp') {

       """
        fastp --in1=${name}.clean.fastq.gz \\
	      --out1=${name}_trimmed_clean.fastq.gz \\
	      --report_title=${name}_FASTP \\
	      --json=${name}_fastp.json  \\
	      --html=${name}_fastp.html
       """

       }
    }else { // we have paired end now
       if(adapter_removal_tool == 'fastp') {

       """
            fastp --in1=${name}_1.clean.fastq.gz  --in2=${name}_2.clean.fastq.gz  \\
              --out1=${name}_1_trimmed_clean.fastq.gz  --out2=${name}_2_trimmed_clean.fastq.gz \\
              --report_title=${name}_FASTP \\
              --json=${name}_fastp.json  \\
              --html=${name}_fastp.html
       """

       }
    }
} // end of Step 3 - Process remove_adapters - channel(s) output analysis_ready_reads, post_trimming_qc

/*
 * STEP 4a - FastQC After
*/
process fastqc_after {
    label 'mid_cm'
    echo true
    tag "$name"
    //publishDir "${params.outdir}/qc/fastqc/", saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    set val(name), file(reads) from post_trimming_qc

    output:
    file "*_fastqc.{zip,html}" into fastqc_after_results

    script:
    """
    fastqc -t 16 -q $reads
    """
}


/*
 * STEP 4b - Salmon quant
*/
my_salmon_index = file(params.salmon_index_path)
println my_salmon_index

process quant {
    label 'mid_cm'
    echo true
    tag "$pair_id"
    time '20h'
    //publishDir "${params.outdir}/quant", mode:'copy'

    input:
    file(my_salmon_index)
    set pair_id, file(reads) from analysis_ready_reads

    output:
    file "salmon_${pair_id}" into quant_ch, quant_ch_summary

    script:
    if (params.singleEnd){
       if(my_quant_method == 'salmon') {

       """
		salmon quant --threads 16 --seqBias --validateMappings --numBootstraps 100 -l A \\
	       -i $my_salmon_index   \\
	       -r ${pair_id}_trimmed_clean.fastq.gz                  \\
	       -o salmon_${pair_id}

       """

       }
    }else { // we have paired end now

     if(my_quant_method == 'salmon') {

       """
        salmon quant --threads 16 --seqBias --validateMappings --numBootstraps 100 -l A \\
               -i $my_salmon_index   \\
               -1 ${pair_id}_1_trimmed_clean.fastq.gz  -2 ${pair_id}_2_trimmed_clean.fastq.gz   \\
               -o salmon_${pair_id}

       """
       }
      }
}

/*
 * STEP 5 - MultiQC
*/
Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }

process multiQC {
    label 'high_cm'
    echo true
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    input:
    file multiqc_config from ch_config_for_multiqc
	file ('pre_fastqc/*') from  preqc_results.collect().ifEmpty([])
    file ('bbmap_stats/*') from contaminant_stats.collect().ifEmpty([])
    file ('adapter_removal/*') from adapter_removal_stats.collect().ifEmpty([])
    file ('fastqc_after/*') from fastqc_after_results.collect().ifEmpty([])
    file ('quant/*') from quant_ch.collect().ifEmpty([])
 
    output:
    file "*_report.html" into multiqc_report
    file "*_data" 
	
    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config --interactive
    """
}

/*
* TxImport - 6
*/

if(params.tx2gene_file) {
  my_tx2gene = file(params.tx2gene_file)
}

process txImport_summary {
    label 'mid_cm'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    file 'quant/*' from quant_ch_summary.collect()
    file (gene) from  my_tx2gene

    output:
    file "*.csv" 

    script:
    if(params.tx2gene_file) {
        """
        forAWS_getSalmonQuants_txImport.r  quant $gene
        """
    } else {
        """
        forAWS_getSalmonQuants_txImport.r  quant
        """
    }
}
