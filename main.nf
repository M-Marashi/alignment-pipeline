#!/usr/bin/env nextflow

// Set default parameters
params.accession = 'M21012'
params.in        = 'data'
params.glob      = '*.{fa,fasta}'
params.out       = 'results'

// Step 1: Download reference genome
process download_ref {
    publishDir params.out
    conda 'bioconda::entrez-direct=24.0'
    input:
        val accession
    output:
        path "${accession}.fasta"
    script:
        """
        esearch -db nucleotide -query "${accession}" | efetch -format fasta > ${accession}.fasta
        """
}

// Step 2: Combine genome FASTA files
process combine_fastas {
    publishDir params.out
    input:
        path genome_fastas
    output:
        path "combined_genomes.fasta"
    script:
        """
        cat ${genome_fastas.join(' ')} > combined_genomes.fasta
        """
}

// Step 3: MAFFT alignment
process mafft_align {
    publishDir params.out
    conda 'bioconda::mafft=7.525'
    input:
        path combined_fasta
    output:
        path "alignment.fasta"
    script:
        """
        mafft --auto --thread 2 ${combined_fasta} > alignment.fasta
        """
}

// Step 4: trimal cleanup and HTML report
process trimal_cleanup {
    publishDir params.out
    conda 'bioconda::trimal=1.5.0'
    input:
        path alignment
    output:
        path "alignment_trimmed.fasta"
        path "alignment_trimmed.html"
    script:
        """
        trimal -in alignment.fasta -out alignment_trimmed.fasta -automated1 -htmlout alignment_trimmed.html
        """
}

// Workflow definition
workflow {
    download_ref(params.accession)
    Channel
        .fromPath("${params.in}/${params.glob}")
        .collect()
        .set { genome_fastas_ch }

    combine_fastas(genome_fastas_ch)
    mafft_align(combine_fastas.out)
    trimal_cleanup(mafft_align.out)
}
