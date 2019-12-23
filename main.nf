import groovy.json.JsonOutput

if(!(params.sample_information && params.fastq_list)){
    println "'sample_information' or 'fastq_list' is undefined"
    println "provide parameters using '--params-file params.json'"
    System.exit(1)
}

def maybe_local(fname){
    // Address the special case of using test files in this project
    // when running in batchman, or more generally, run-from-git.
    if(file(fname).exists() || fname.startsWith('s3://')){
        return file(fname)
    }else{
        file("$workflow.projectDir/" + fname)
    }
}

fastq_list = maybe_local(params.fastq_list)
sample_information = maybe_local(params.sample_information)

// iterate over list of input files, split sampleid from filenames,
// and arrange into a sequence of (sampleid, I1, [I2], R1, R2)
// TODO: use of 'maybe_local()' is untested with s3 objects
sample_info = Channel.fromPath(sample_information)
Channel.fromPath(fastq_list)
    .splitText()
    .map { it.trim() }
    .map { maybe_local(it) }
    .map { [(it.fileName =~ /(^[-a-zA-Z0-9]+)/)[0][0], it ] }
    .into{ sample_map1; sample_map2 }

to_barcodecop = sample_map1
    .groupTuple()
    .map { [ it[0], it[1].sort() ] }
    .map { it.flatten() }

// to_barcodecop.println { "Received: $it" }

to_plot_quality = sample_map2
    .filter{ it -> it[1].fileName =~ /_R[12]_/ }
    .groupTuple()
    .map { [ it[0], it[1].sort() ] }
    .map { it.flatten() }

// to_plot_quality.println { "Received: $it" }

process read_manifest {

    input:
        file(sample_info) from sample_info
        file(fastq_files) from Channel.fromPath(fastq_list)

    output:
        file("batches.csv") into batches

    // publishDir params.output, overwrite: true

    """
    manifest.py --outfile batches.csv ${fastq_files} --manifest ${sample_info} \
        --index-file-type ${params.index_file_type}
    """
}

process plot_quality {

    label 'med_cpu_mem'

    input:
        tuple sampleid, file(R1), file(R2) from to_plot_quality

    output:
        file("${sampleid}.png")

    publishDir "${params.output}/qplots/", overwrite: true

    """
    dada2_plot_quality.R ${R1} ${R2} -o ${sampleid}.png \
        --trim-left ${params.trim_left} \
        --f-trunc ${params.f_trunc} \
        --r-trunc ${params.r_trunc}

    """
}

if(params.index_file_type == 'dual'){
    process barcodecop_dual {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(I1), file(I2), file(R1), file(R2) from to_barcodecop

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_filtered
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true

        """
        barcodecop --fastq ${R1} ${I1} ${I2} \
            --outfile ${sampleid}_R1_.fq.gz --read-counts ${sampleid}_R1_counts.csv \
            --match-filter --qual-filter
        barcodecop --fastq ${R2} ${I1} ${I2} \
            --outfile ${sampleid}_R2_.fq.gz --read-counts ${sampleid}_R2_counts.csv \
            --match-filter --qual-filter
        """
    }
}else if(params.index_file_type == 'single'){
    process barcodecop_single {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(I1), file(R1), file(R2) from to_barcodecop

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_filtered
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true

        """
        barcodecop --fastq ${R1} ${I1} \
            --outfile ${sampleid}_R1_.fq.gz --read-counts ${sampleid}_R1_counts.csv \
            --match-filter --qual-filter
        barcodecop --fastq ${R2} ${I1} \
            --outfile ${sampleid}_R2_.fq.gz --read-counts ${sampleid}_R2_counts.csv \
            --match-filter --qual-filter
        """
    }
}

process bcop_counts_concat {

    input:
        file("counts*.csv") from bcop_counts.collect()

    output:
        file("bcop_counts.csv") into bcop_counts_concat

    // publishDir "${params.output}", overwrite: true

    // TODO: barcodecop should have --sampleid argument to pass through to counts

    """
    echo "sampleid,raw,barcodecop" > bcop_counts.csv
    cat counts*.csv | sed 's/_R[12]_.fq.gz//g' | sort | uniq >> bcop_counts.csv
    """
}

// bcop_filtered.println { "Received: $it" }

// Join read counts with names of files filtered by barcodecop,
// discard empty files, and return sequence of (sampleid, R1, R2).
to_filter = bcop_counts_concat
    .splitCsv(header: true)
    .filter{ it['barcodecop'].toInteger() >= params.min_reads }
    .cross(bcop_filtered)
    .map{ it[1] }


process filter_and_trim {

    label 'med_cpu_mem'

    input:
        tuple sampleid, file(R1), file(R2) from to_filter

    output:
        tuple sampleid, file("${sampleid}_R1_filt.fq.gz"), file("${sampleid}_R2_filt.fq.gz") into filtered_trimmed

    // publishDir "${params.output}/filtered/", overwrite: true

    """
    dada2_filter_and_trim.R \
        --infiles ${R1} ${R2} \
        --outfiles ${sampleid}_R1_filt.fq.gz ${sampleid}_R2_filt.fq.gz \
        --trim-left ${params.trim_left} \
        --f-trunc ${params.f_trunc} \
        --r-trunc ${params.r_trunc} \
        --truncq 2
    """
}


// [sampleid, batch, R1, R2]
batches
    .splitCsv(header: false)
    .join(filtered_trimmed)
    .into { to_learn_errors ; to_dereplicate }


process learn_errors {

    label 'med_cpu_mem'

    input:
        tuple batch, file("R1_*.fastq.gz"), file("R2_*.fastq.gz") from to_learn_errors.map{ it[1, 2, 3] }.groupTuple()

    output:
        file("error_model_${batch}.rds") into error_models
    file("error_model_${batch}.png") into error_model_plots

    publishDir "${params.output}/error_models", overwrite: true

    """
    ls -1 R1_*.fastq.gz > R1.txt
    ls -1 R2_*.fastq.gz > R2.txt
    dada2_learn_errors.R --r1 R1.txt --r2 R2.txt \
        --model error_model_${batch}.rds \
        --plots error_model_${batch}.png
    """
}


process dada_dereplicate {

    label 'med_cpu_mem'

    input:
        tuple sampleid, batch, file(R1), file(R2) from to_dereplicate
    file("") from error_models.collect()

    output:
        file("dada.rds") into dada_data
    file("seqtab.csv") into dada_seqtab
    file("counts.csv") into dada_counts
    file("overlaps.csv") into dada_overlaps

    publishDir "${params.output}/dada/${sampleid}/", overwrite: true

    """
    dada2_dada.R ${R1} ${R2} --errors error_model_${batch}.rds \
        --sampleid ${sampleid} \
        --self-consist ${params.self_consist} \
        --data dada.rds \
        --seqtab seqtab.csv \
        --counts counts.csv \
        --overlaps overlaps.csv
    """
}


process combined_overlaps {

    input:
        file("overlaps_*.csv") from dada_overlaps.collect()

    output:
        file("overlaps.csv")

    publishDir params.output, overwrite: true

    """
    csvcat.sh overlaps_*.csv > overlaps.csv
    """
}


process dada_counts_concat {

    input:
        file("*.csv") from dada_counts.collect()

    output:
        file("dada_counts.csv") into dada_counts_concat

    // publishDir params.output, overwrite: true

    """
    csvcat.sh *.csv > dada_counts.csv
    """
}


process write_seqs {

    input:
        file("seqtab_*.csv") from dada_seqtab.collect()

    output:
        file("seqs.fasta") into seqs
    file("specimen_map.csv")
    file("sv_table.csv")
    file("sv_table_long.csv")
    file("weights.csv") into weights

    publishDir params.output, overwrite: true

    """
    write_seqs.py seqtab_*.csv \
        --seqs seqs.fasta \
        --specimen-map specimen_map.csv \
        --sv-table sv_table.csv \
        --sv-table-long sv_table_long.csv \
        --weights weights.csv
    """
}

// clone channel so that it can be consumed twice
seqs.into { seqs_to_align; seqs_to_filter }


process cmalign {

    label 'med_cpu_mem'

    input:
        file("seqs.fasta") from seqs_to_align
    file('ssu.cm') from file("$workflow.projectDir/data/ssu-align-0.1.1-bacteria-0p1.cm")

    output:
        file("seqs.sto")
    file("sv_aln_scores.txt") into aln_scores

    publishDir params.output, overwrite: true

    """
    cmalign --dnaout --noprob \
        -o seqs.sto --sfile sv_aln_scores.txt ssu.cm seqs.fasta
    """
}


process filter_16s {

    input:
        file("seqs.fasta") from seqs_to_filter
    file("sv_aln_scores.txt") from aln_scores
    file("weights.csv") from weights

    output:
        file("16s.fasta")
    file("not16s.fasta")
    file("16s_outcomes.csv")
    file("16s_counts.csv") into is_16s_counts

    publishDir params.output, overwrite: true

    """
    filter_16s.py seqs.fasta sv_aln_scores.txt --weights weights.csv \
        --min-bit-score 0 \
        --passing 16s.fasta \
        --failing not16s.fasta \
        --outcomes 16s_outcomes.csv \
        --counts 16s_counts.csv
    """
}


process join_counts {

    input:
        file("bcop.csv") from bcop_counts_concat
    file("dada.csv") from dada_counts_concat
    file("16s.csv") from is_16s_counts

    output:
        file("counts.csv")

    publishDir params.output, overwrite: true

    """
    ljoin.R bcop.csv dada.csv 16s.csv -o counts.csv
    """
}


process save_params {

    input:
        val parameters from JsonOutput.prettyPrint(JsonOutput.toJson(params))

    output:
        file('params.json')

    publishDir params.output, overwrite: true

    """
    echo "${parameters}" > params.json
    """

}
