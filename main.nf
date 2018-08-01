#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                         NF-CAGEseq
========================================================================================
 NF-CAGEseq Analysis Pipeline. Started 2018-03-09.
 #### Homepage / Documentation
 https://github.com/KevinMenden/NF-CAGEseq
 #### Authors
 Kevin Menden KevinMenden <kevin.menden@dzne.de> - https://github.com/KevinMenden>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     NF-CAGEseq v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run KevinMenden/NF-CAGEseq --reads '*.fastq.gz' -profile docker

    For paired-end reads specify --reads '*{1,2}.fastq.gz' or a similar pattern matching your filenames

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. docker / aws

    Options:
      --pairedEnd                   Specifies that the input is paired end reads (Currently not supported)
      --aligner                     Specify which aligner should be used. One of: bowtie | star . Default: bowtie

    Trimming:
      --cutEcop                     [true|false] Whether the 5' EcoP15I regognition site should be removed. Default is true.
      --cutLinker                   [true|false] Whether the 3' linker should be removed. Default is true.
      --trimming                    [true | false] Whether trimming should be performed. Default is true.

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF reference
      --genome                      Name of iGenomes genome

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}



/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = params.version

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.genome = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.reads = "data/*{1,2}.fastq.gz"
params.outdir = './results'
params.aligner = 'bowtie'
params.gtf = false
params.min_aln_length = 15
params.star_index = false
params.saveReference = false
params.cutEcop = true
params.cutLinker = true
params.trimming = true
params.cutG = true

multiqc_config = file(params.multiqc_config)

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
params.pairedEnd = false
Channel
    .fromFilePairs( params.reads, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_bam; read_files_trimming }

log.info"""
${read_files_fastqc.length()}

"""


// Header log info
log.info "========================================="
log.info " NF-CAGEseq v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.pairedEnd ? 'Paired-End' : 'Single-End'
summary['Aligner']      = params.aligner
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/**
 * Load and validate inputs
 */
// Load FASTA file
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
// Load GTF file
if( params.gtf ){
    Channel
            .fromPath(params.gtf)
            .ifEmpty { exit 1, "GTF annotation file not found: ${params.gtf}" }
            .into { gtf_makeSTARindex; gtf_star; gtf_bowtie; gtf_buildBowtieIndex}
}
// Load STAR index
if ( params.star_index){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    STAR --version > v_star.txt
    cutadapt --version > v_cutadapt.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}



/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
STEP 2 - Create queryname sorted, unmapped BAM from FastQ files
 */
process fastqToBam {
    tag "$name"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    set val(name), file(reads), from read_files_bam

    output:
    file "*bam" into fastq_to_bam_results

    script:
    """
    picard FastqToSam F1=${reads[0]} F2=${reads[1]} O=test.bam RG=null
    """

}

/*
STEP 3 Tag bam file
Use R1 to tag each read according to cell and molecular barcode
Filter for reads with low quality barcodes
 */
process tagbam {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    file ua_bam from fastq_to_bam_results

    script:
    """
    TagBamWithReadSequenceExtended \\
    INPUT=$ua_bam \\
    OUTPUT=unaligned_tagged_Cell.bam \\
    SUMMARY=unaligned_tagged_Cell.bam_summary.txt \\
    BASE_RANGE=1-12 \\
    BASE_QUALITY=10 \\
    BARCODED_READ=1 \\
    DISCARD_READ=False \\
    TAG_NAME=XC \\
    NUM_BASES_BELOW_QUALITY=1


    TagBamWithReadSequenceExtended \\
    INPUT=unaligned_tagged_Cell.bam \\
    OUTPUT=unaligned_tagged_CellMolecular.bam \\
    SUMMARY=unaligned_tagged_CellMolecular.bam_summary.txt \\
    BASE_RANGE=13-20 \\
    BASE_QUALITY=10 \\
    BARCODED_READ=1 \\
    DISCARD_READ=True \\
    TAG_NAME=XM \\
    NUM_BASES_BELOW_QUALITY=1

    FilerBam \\
    TAG_REJECT=XQ \\
    INPUT=unaligned_tagged_CellMolecular.bam \\
    OUTPUT=unaligned_tagged_filtered.bam
    """

    output:
    file "*tagged_filtered.bam" into tagged_bam
    ""
}

/*
STEP 4 Read trimming
 */
process trimming {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    file tbam from tagged_bam

    script:
    """
    TrimStartingSequence \\
    INPUT=$tbam \\
    OUTPUT=unaligned_tagged_trimmed_smart.bam \\
    OUTPUT_SUMMARY=adapter_trimming_report.txt \\
    SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \\
    MISMATCHES=0 \\
    NUM_BASES=5
    
    PolyATrimmer \\
    INPUT=unaligned_tagged_trimmed_smart.bam \\
    OUTPUT=unaligned_mc_tagged_polyA_filtered.bam \\
    OUTPUT_SUMMARY=polyA_trimming_report.txt \\
    MISMATCHES=0 \\
    NUM_BASES=6
    """

    output:
    file "*mc_tagged_polyA_filtered.bam" into filtered_bam
    file "*" into trimming results
}

/*
STEP 5 Convert to Fastq
 */
process bamToFastq {
    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    file fbam from filtered_bam

    script:
    """
    picard SamToFastq \\
    INPUT=$fbam \\
    FASTQ=unaligned_mc_tagged_filtered.fastq
    """

    output:
    file "*fastq" into tagged_filtered_fastq
}


/*
 * PREPROCESSING - Build STAR index
 */
if(!params.star_index && fasta){
    process makeSTARindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta
        file gtf from gtf_makeSTARindex

        output:
        file "star" into star_index

        script:
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile $gtf \\
            --sjdbOverhang 50 \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }
}


/**
 * READ MAPPING
 */
    process star {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else  filename }

        input:
        file reads from processed_reads
        file index from star_index.collect()
        file gtf from gtf_star.collect()

        output:
        file '*.bam' into star_aligned
        file "*.out" into alignment_logs
        file "*SJ.out.tab"
        file "*Log.out" into star_log

        script:
        prefix = reads[0].toString() - ~/(.trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        STAR --genomeDir $index \\
            --sjdbGTFfile $gtf \\
            --readFilesIn $reads  \\
            --runThreadN ${task.cpus} \\
            --outSAMtype BAM SortedByCoordinate \\
            --readFilesCommand zcat \\
            --runDirPerm All_RWX \\
            --outFileNamePrefix $prefix \\
            --outFilterMatchNmin ${params.min_aln_length}
         """ }
// Split Star results
star_aligned.into { bam_count; bam_rseqc }


/**
 * Mapping QC
 */
process rseqc {
    tag "${bam_rseqc.baseName}"
    publishDir "${params.outdir}/rseqc", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("bam_stat.txt") > 0)          "bam_stat/$filename"
            }

    input:
    file bam_rseqc

    output:
    file "*.txt" into rseqc_results

    script:
    """
    bam_stat.py -i ${bam_rseqc} > ${bam_rseqc.baseName}.bam_stat.txt
    """

}




/**
 * STEP7 MultiQC
 */

// MultiQC of non-trimmed fastq files
process multiqc_notrim {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC/notrim", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()

    output:
    file "*multiqc_report.html" into multiqc_report_notrim
    file "*_data"
    val prefix into multiqc_prefix_notrim

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}

// MultiQC of all the results
process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc: 'trimmed/*') from trimmed_fastqc_results.collect()
    file ('trimmed/*') from cutadapt_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file ('rseqc/*') from rseqc_results.collect()
    file ('alignment/*') from alignment_logs.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


/*
 * Completion notification
 */
workflow.onComplete {
    log.info "[NF-CAGEseq] Pipeline Complete"
}
