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
     SeqWell v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run kevinmenden/seqwell --reads '*.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. docker / aws

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF reference

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
params.gtf = false
params.star_index = false
params.build_dict = false
params.saveReference = false
params.min_genes = 50

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
Channel
    .fromFilePairs( params.reads, size: 2)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_bam }

log.info"""
${read_files_fastqc.length()}

"""


// Header log info
log.info "========================================="
log.info " SeqWell v${version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
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
            .into { gtf_makeSTARindex; gtf_star; gtf_tag_exon; gtf_buildBowtieIndex}
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
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/**
 * FastQC
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


/**
 * Create query name sorted, unmapped BAM from FastQ files
 */
process fastqToBam {
    tag "$name"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_bam

    output:
    file "*bam" into fastq_to_bam_results

    script:
    """
    picard FastqToSam F1=${reads[0]} F2=${reads[1]} O=${reads[0].baseName}.bam RG=null
    """

}

/**
 * Tag reads according to cell and molecular barcodes
 * Filter for reads with low quality barcodes
 */
process tagbam {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    file ua_bam from fastq_to_bam_results


    output:
    file "*tagged_filtered.bam" into tagged_bam

    script:
    """
    TagBamWithReadSequenceExtended \\
    INPUT=$ua_bam \\
    OUTPUT=${ua_bam.baseName}.unaligned_tagged_Cell.bam \\
    SUMMARY=${ua_bam.baseName}.unaligned_tagged_Cell.bam_summary.txt \\
    BASE_RANGE=1-12 \\
    BASE_QUALITY=10 \\
    BARCODED_READ=1 \\
    DISCARD_READ=False \\
    TAG_NAME=XC \\
    NUM_BASES_BELOW_QUALITY=1


    TagBamWithReadSequenceExtended \\
    INPUT=${ua_bam.baseName}.unaligned_tagged_Cell.bam \\
    OUTPUT=${ua_bam.baseName}.unaligned_tagged_CellMolecular.bam \\
    SUMMARY=unaligned_tagged_CellMolecular.bam_summary.txt \\
    BASE_RANGE=13-20 \\
    BASE_QUALITY=10 \\
    BARCODED_READ=1 \\
    DISCARD_READ=True \\
    TAG_NAME=XM \\
    NUM_BASES_BELOW_QUALITY=1

    FilerBam \\
    TAG_REJECT=XQ \\
    INPUT=${ua_bam.baseName}.unaligned_tagged_CellMolecular.bam \\
    OUTPUT=${ua_bam.baseName}.unaligned_tagged_filtered.bam
    """
}


/**
 * Trim reads at 5'-end
 */
process trimming {
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    file tbam from tagged_bam

    output:
    file "*polyA_filtered.bam" into filtered_bam
    file "*" into trimming_results

    script:
    """
    TrimStartingSequence \\
    INPUT=$tbam \\
    OUTPUT=${tbam.baseName}.trimmed_smart.bam \\
    OUTPUT_SUMMARY=${tbam.baseName}.adapter_trimming_report.txt \\
    SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \\
    MISMATCHES=0 \\
    NUM_BASES=5
    
    PolyATrimmer \\
    INPUT=${tbam.baseName}.trimmed_smart.bam \\
    OUTPUT=${tbam.baseName}.polyA_filtered.bam \\
    OUTPUT_SUMMARY=${tbam.baseName}.polyA_trimming_report.txt \\
    MISMATCHES=0 \\
    NUM_BASES=6
    """


}
filtered_bam.into { filtered_bam_fastq; filtered_bam_merge }

/**
 * Convert back to FastQ
 */
process bamToFastq {
    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    file fbam from filtered_bam_fastq

    output:
    file "*fastq" into tagged_filtered_fastq

    script:
    """
    picard SamToFastq \\
    INPUT=$fbam \\
    FASTQ=${$fbam}.fastq
    """


}


/**
 * Build STAR index if not available
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
            --sjdbOverhang 49 \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }
}


/**
 * Map reads using STAR
 */
    process star {
        tag "$prefix"
        publishDir "${params.outdir}/STAR", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".bam") == -1) "logs/$filename"
                    else  filename }

        input:
        file reads from tagged_filtered_fastq
        file index from star_index.collect()

        output:
        file '*.sam' into star_aligned
        file "*.out" into alignment_logs
        file "*Log.out" into star_log

        script:
        prefix = reads.toString() - ~/(.trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
        """
        STAR --genomeDir $index \\
            --readFilesIn $reads  \\
            --runThreadN ${task.cpus} \\
            --outFileNamePrefix $prefix \\
         """ }

/**
 * Build sequence dict if not available
 */
if(params.build_dict){
    process buildSequenceDict {
        publishDir "${params.outdir}/reference_genome", mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "*.fa*" into fasta_ref
        file "*.dict" into fasta_dict

        script:
        """
        picard 
        """
    }
}

/**
 * Sort Sam file and merge alignment
 */
process sort_and_merge {
    publishDir "${params.outdir}/STAR", mode: 'copy'

    input:
    file mapped_sam from star_aligned
    file fbam from filtered_bam_merge
    file fasta from fasta

    output:
    file "*bam" into merged_bam

    script:
    """
    picard SortSam \\
    I=$mapped_sam \\
    O=${mapped_sam.baseName}.sorted.bam \\
    SO=queryname
    
    picard MergeBamAlignment \\
    REFERENCE_SEQUENCE=$fasta \\
    UNMAPPED_BAM=$fbam \\
    ALIGNED_BAM=${mapped_sam.baseName}.sorted.bam \\
    OUTPUT=${mapped_sam.baseName}.merged.bam \\
    INCLUDE_SECONDARY_ALIGNMENTS=false \\
    PAIRED_RUN=false    
    """
}

/**
 * Tag Read with Gene Exon
 */
process tagGeneExon {
    publishDir "${params.outdir}/Counts", mode: 'copy'

    input:
    file mbam from merged_bam
    file gtf from gtf_tag_exon

    output:
    file "*bam" into exon_tagged_bam

    script:
    """
    TagReadWithGeneExon \\
    I=$mbam \\
    O=${mbam.baseName}.exonTagged.bam \\
    ANNOTATIONS_FILE=$gtf \\
    TAG=GE
    """
}

/**
 * Create Digital Gene Expression matrix
 */
process digitalGeneExpression {
    publishDir "${params.outdir}/Counts", mode: 'copy'

    input:
    file etbam from exon_tagged_bam

    output:
    file "*" into dge_matrix

    script:
    """
    DigitalExpression \\
    I=$etbam \\
    O=${etbam.baseName}.dge.txt.gz \\
    SUMMARY=${etbam.baseName}.dge.summary.txt
    MIN_NUM_GENES_PER_CELL=$params.min_genes 
    """
}


/**
 * STEP7 MultiQC
 */

// MultiQC of all the results
process multiqc {
    tag "$prefix"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml
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
    log.info "SeqWell Pipeline Complete"
}
