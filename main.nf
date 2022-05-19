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

if (params.containsKey("downsample") && params.downsample) {
  head_cmd = "--head $params.downsample"
  to_plot_quality = to_plot_quality
    .map { [ it[0],
             it[1].splitFastq(
                by: params.downsample,
                compress: true,
                file: true,
                limit: params.downsample)[0],
             it[2].splitFastq(
                by: params.downsample,
                compress: true,
                file: true,
                limit: params.downsample)[0] ] }
} else {
    head_cmd = ""
}

// to_plot_quality.println { "Received: $it" }

fastq_file_val = Channel.value(fastq_list)

process copy_filelist {
    input:
        file(fastq_files) from fastq_file_val

    output:
        file("fastq_list.csv")

    publishDir params.output, overwrite: true, mode: 'copy'

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

    publishDir params.output, overwrite: true, mode: 'copy'

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

    publishDir "${params.output}/qplots/", overwrite: true, mode: 'copy'

    """
    dada2_plot_quality.R ${R1} ${R2} --params dada_params.json -o ${sampleid}.png
    """
}

if(params.index_file_type == 'dual'){
    process barcodecop_dual {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(I1), file(I2), file(R1), file(R2) from to_barcodecop
            val head_cmd

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_filtered
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true, mode: 'copy'

        """
        barcodecop --fastq ${R1} ${I1} ${I2} \
            --outfile ${sampleid}_R1_.fq.gz --read-counts ${sampleid}_R1_counts.csv \
            --match-filter --qual-filter ${head_cmd}
        barcodecop --fastq ${R2} ${I1} ${I2} \
            --outfile ${sampleid}_R2_.fq.gz --read-counts ${sampleid}_R2_counts.csv \
            --match-filter --qual-filter ${head_cmd}
        """
    }
}else if(params.index_file_type == 'single'){
    process barcodecop_single {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(I1), file(R1), file(R2) from to_barcodecop
            val head_cmd

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_filtered
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true, mode: 'copy'

        """
        barcodecop --fastq ${R1} ${I1} \
            --outfile ${sampleid}_R1_.fq.gz --read-counts ${sampleid}_R1_counts.csv \
            --match-filter --qual-filter ${head_cmd}
        barcodecop --fastq ${R2} ${I1} \
            --outfile ${sampleid}_R2_.fq.gz --read-counts ${sampleid}_R2_counts.csv \
            --match-filter --qual-filter ${head_cmd}
        """
    }
}else if(params.index_file_type == 'none') {
    process read_counts {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(R1), file(R2) from to_barcodecop

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_filtered
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true, mode: 'copy'

        """
        read_counts.py ${head_cmd} ${R1} ${sampleid}_R1_.fq.gz ${sampleid}_R1_counts.csv
        read_counts.py ${head_cmd} ${R2} ${sampleid}_R2_.fq.gz ${sampleid}_R2_counts.csv
        """
    }
}

process bcop_counts_concat {

    input:
        file("counts*.csv") from bcop_counts.collect()

    output:
        file("bcop_counts.csv") into bcop_counts_concat

    // publishDir "${params.output}", overwrite: true, mode: 'copy'

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

    // publishDir "${params.output}/filtered/", overwrite: true, mode: 'copy'

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

    publishDir "${params.output}/error_models", overwrite: true, mode: 'copy'

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
        val sampleid into dada_dereplicate_samples

    publishDir "${params.output}/dada/${sampleid}/", overwrite: true, mode: 'copy'

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

process dada_get_unmerged {

    label 'med_cpu_mem'

    input:
        val sampleid from dada_dereplicate_samples
        file("dada_params.json") from maybe_local(params.dada_params)
        file dada_rds from dada_data

    output:
        file("unmerged_*.fasta") into dada_unmerged

    publishDir "${params.output}/dada/${sampleid}/", overwrite: true, mode: 'copy'

    """
    get_unmerged.R ${dada_rds} \
        --forward-seqs unmerged_F.fasta --reverse-seqs unmerged_R.fasta
    """
}


process combined_overlaps {

    input:
        file("overlaps_*.csv") from dada_overlaps.collect()

    output:
        file("overlaps.csv")

    publishDir params.output, overwrite: true, mode: 'copy'

    """
    csvcat.sh overlaps_*.csv > overlaps.csv
    """
}


process dada_counts_concat {

    input:
        file("*.csv") from dada_counts.collect()

    output:
        file("dada_counts.csv") into dada_counts_concat

    // publishDir params.output, overwrite: true, mode: 'copy'

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
        file("sv_table_long.csv") into sv_table_long
        file("weights.csv") into weights

    publishDir params.output, overwrite: true, mode: 'copy'

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
        file('ssu.cm') from maybe_local(params.alignment_model)

    output:
        file("sv_aln_scores.txt") into aln_scores

    publishDir params.output, overwrite: true, mode: 'copy'

    """
    cmsearch --hmmonly --noali --tblout sv_aln_scores.txt ssu.cm seqs.fasta
    """
}


process filter_svs {
    input:
        file("seqs.fasta") from seqs_to_filter
        file("sv_aln_scores.txt") from aln_scores
        file("weights.csv") from weights_to_filter

    output:
        file("passed.fasta") into passed
        file("failed.fasta")
        file("outcomes.csv")
        file("counts.csv") into passed_counts

    publishDir "${params.output}/filter_svs/", overwrite: true, mode: 'copy'

    """
    filter_svs.py \
        --counts counts.csv \
        --failing failed.fasta \
        --min-bit-score 0 \
        --outcomes outcomes.csv \
        --passing passed.fasta \
        --weights weights.csv \
        seqs.fasta sv_aln_scores.txt
    """
}

if(params.containsKey("bidirectional") && params.bidirectional){
    process vsearch_svs {
        // Append size/weight to sequence headers: ";size=integer"
        // so vsearch can sort by weight
        input:
            file("seqs.fasta") from passed
            file("sv_table_long.csv") from sv_table_long

        output:
            file("clusters.uc") into clusters

        publishDir params.output, overwrite: true, mode: 'copy'

        """
        append_size.py seqs.fasta sv_table_long.csv |
        vsearch --cluster_size - --uc clusters.uc --id 1.0 --iddef 2 --xsize
        """
    }

    process combine_svs {
        input:
            file("seqs.fasta") from passed
            file("sv_table_long.csv") from sv_table_long
            file("clusters.uc") from clusters

        output:
            file("seqtab.csv") into combined

        publishDir params.output, overwrite: true, mode: 'copy'

        """
        combine_svs.py --out seqtab.csv \
        seqs.fasta sv_table_long.csv clusters.uc
        """
    }


    process write_complemented_seqs {
        // NOTE: sv names will be regenerated and will
        // not be traceable to earlier steps in the pipeline
        input:
            file("seqtab.csv") from combined

        output:
            file("seqs.fasta")
            file("specimen_map.csv")
            file("sv_table.csv")
            file("sv_table_long.csv")
            file("weights.csv")

        publishDir params.output, overwrite: true, mode: 'copy'

        """
        write_seqs.py \
            --seqs seqs.fasta \
            --specimen-map specimen_map.csv \
            --sv-table sv_table.csv \
            --sv-table-long sv_table_long.csv \
            --weights weights.csv \
            seqtab.csv
        """
    }
}

process join_counts {
    input:
        file("bcop.csv") from bcop_counts_concat
        file("dada.csv") from dada_counts_concat
        file("passed.csv") from passed_counts

    output:
        file("counts.csv")

    publishDir params.output, overwrite: true, mode: 'copy'

    """
    ljoin.R bcop.csv dada.csv passed.csv -o counts.csv
    """
}


process save_params {

    input:
        val parameters from JsonOutput.prettyPrint(JsonOutput.toJson(params))

    output:
        file('params.json')

    publishDir params.output, overwrite: true, mode: 'copy'

    """
cat <<EOF > params.json
${parameters}
EOF
    """
}
