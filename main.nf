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

fastq_file_val = Channel.value(fastq_list)

process copy_filelist {
    input:
        file(fastq_files) from fastq_file_val

    output:
        file("fastq_list.csv")

    publishDir params.output, overwrite: true

    """
    cp ${fastq_files} fastq_list.csv
    """
}

process read_manifest {

    input:
        file(sample_info) from sample_info
        file(fastq_files) from fastq_file_val

    output:
        file("batches.csv") into batches
        file("sample_information.csv")

    publishDir params.output, overwrite: true

    """
    manifest.py ${fastq_files} --manifest ${sample_info} \
        --batches batches.csv \
        --sample-info sample_information.csv \
        --index-file-type ${params.index_file_type}
    """
}

process plot_quality {

    label 'med_cpu_mem'

    input:
        tuple sampleid, file(R1), file(R2) from to_plot_quality
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        file("${sampleid}.png")

    publishDir "${params.output}/qplots/", overwrite: true

    """
    dada2_plot_quality.R ${R1} ${R2} --params dada_params.json -o ${sampleid}.png
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
}else if(params.index_file_type == 'none') {
    process read_counts {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") from to_barcodecop

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_filtered
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true

        """
        read_counts.py ${sampleid}_R1_.fq.gz ${sampleid}_R1_counts.csv
        read_counts.py ${sampleid}_R2_.fq.gz ${sampleid}_R2_counts.csv
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
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        tuple sampleid, file("${sampleid}_R1_filt.fq.gz"), file("${sampleid}_R2_filt.fq.gz") into filtered_trimmed

    // publishDir "${params.output}/filtered/", overwrite: true

    """
    dada2_filter_and_trim.R \
        --infiles ${R1} ${R2} \
        --params dada_params.json \
        --outfiles ${sampleid}_R1_filt.fq.gz ${sampleid}_R2_filt.fq.gz \
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
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        file("dada.rds") into dada_data
        file("seqtab.csv") into dada_seqtab
        file("counts.csv") into dada_counts
        file("overlaps.csv") into dada_overlaps

    publishDir "${params.output}/dada/${sampleid}/", overwrite: true

    """
    dada2_dada.R ${R1} ${R2} \
        --errors error_model_${batch}.rds \
        --sampleid ${sampleid} \
        --params dada_params.json \
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

// clone channels so they can be consumed multiple times
seqs.into { seqs_to_align; seqs_to_filter; seqs_to_be_complemented }
weights.into { weights_to_filter; weights_to_combine }


process cmsearch {

    label 'large_cpu_mem'

    input:
        file("seqs.fasta") from seqs_to_align
        file('ssu.cm') from file("$workflow.projectDir/data/SSU_rRNA_bacteria.cm")

    output:
        file("sv_aln_scores.txt") into aln_scores

    publishDir params.output, overwrite: true

    """
    cmsearch --hmmonly --noali --tblout sv_aln_scores.txt ssu.cm seqs.fasta
    """
}


process filter_16s {

    input:
        file("seqs.fasta") from seqs_to_filter
        file("sv_aln_scores.txt") from aln_scores
        file("weights.csv") from weights_to_filter

    output:
        file("16s.fasta") into seqs_16s
        file("not16s.fasta")
        file("16s_outcomes.csv")
        file("forward_seqs.fasta") into forward_seqs
        file("reverse_seqs.fasta") into reverse_seqs
        file("16s_counts.csv") into is_16s_counts
        file("orientations.csv")

    publishDir params.output, overwrite: true

    """
    filter_16s.py seqs.fasta sv_aln_scores.txt --weights weights.csv \
        --min-bit-score 0 \
        --passing 16s.fasta \
        --failing not16s.fasta \
        --outcomes 16s_outcomes.csv \
        --forwards forward_seqs.fasta \
        --reverses reverse_seqs.fasta \
        --counts 16s_counts.csv \
        --orientations orientations.csv
    """
}

reverse_seqs.into { reverse_seqs_to_vsearch; reverse_seqs_to_combine; reverse_seqs_to_complement }

process vsearch_fwd_rev_svs {
    
    input:
        file("forward_seqs.fasta") from forward_seqs
        file("reverse_seqs.fasta") from reverse_seqs_to_vsearch
    
    output:
        file("vsearch_out.txt") into vsearch_out

    publishDir params.output, overwrite: true

    """
    vsearch --usearch_global reverse_seqs.fasta --db forward_seqs.fasta --strand both --userout vsearch_out.txt --userfields query+target+qstrand+tstrand --id 1.0
    """
}

process combine_svs {

    input:
        file("vsearch_out.txt") from vsearch_out
        file("weights.csv") from weights_to_combine
        file("reverse_seqs.fasta") from reverse_seqs_to_combine

    output:
        file("corrected_weights.csv") into corrected_weights

    """
    combine_svs.py vsearch_out.txt weights.csv reverse_seqs.fasta --corrected_weights corrected_weights.csv
    """
}


process write_complemented_seqs {

    input:
        file("seqs.fasta") from seqs_to_be_complemented
        file("corrected_weights.csv") from corrected_weights

    output:
        file("final_complemented_seqs.fasta")
    
    """
    write_complemented_seqs.py seqs.fasta corrected_weights.csv --final_seqs final_complemented_seqs.fasta
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
cat <<EOF > params.json
${parameters}
EOF
    """
}
