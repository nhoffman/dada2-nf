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
    .into{ sample_map1; sample_map2; to_manifest }

to_fastq_filters = sample_map1
    .groupTuple()
    .map { [ it[0], it[1].sort() ] }
    .map { it.flatten() }

// TODO: filter out fastqs not in manifest before doing work

to_plot_quality = sample_map2
    .filter{ it -> it[1].fileName =~ /_R[12]_/ }
    .groupTuple()
    .map { [ it[0], it[1].sort() ] }
    .map { it.flatten() }

if (params.containsKey("downsample") && params.downsample) {
  head_cmd = "--head $params.downsample"
  to_plot_quality
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
                limit: params.downsample)[0] ] }.set{ to_plot_quality }
} else {
    head_cmd = ""
}

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
        file(fastqs) from to_manifest.collect()

    output:
        file("batches.csv") into batches
        file("sample_information.csv")
        file("counts.csv") into raw_counts

    publishDir "${params.output}/manifest/", overwrite: true, mode: 'copy'

    """
    manifest.py ${fastq_files} --manifest ${sample_info} \
        --batches batches.csv \
        --counts counts.csv \
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
            tuple sampleid, file(I1), file(I2), file(R1), file(R2) from to_fastq_filters
            val head_cmd

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_reads
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        publishDir "${params.output}/barcodecop/", overwrite: true, mode: 'copy'

        """
        barcodecop --fastq ${R1} ${I1} ${I2} \
            --outfile ${sampleid}_R1_.fq.gz --read-counts ${sampleid}_R1_counts.csv \
            --match-filter --allow-empty --qual-filter ${head_cmd}
        barcodecop --fastq ${R2} ${I1} ${I2} \
            --outfile ${sampleid}_R2_.fq.gz --read-counts ${sampleid}_R2_counts.csv \
            --match-filter --allow-empty --qual-filter ${head_cmd}
        """
    }
} else if(params.index_file_type == 'single'){
    process barcodecop_single {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(I1), file(R1), file(R2) from to_fastq_filters
            val head_cmd

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_reads
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true, mode: 'copy'

        """
        barcodecop --fastq ${R1} ${I1} \
            --outfile ${sampleid}_R1_.fq.gz --read-counts ${sampleid}_R1_counts.csv \
            --match-filter --allow-empty --qual-filter ${head_cmd}
        barcodecop --fastq ${R2} ${I1} \
            --outfile ${sampleid}_R2_.fq.gz --read-counts ${sampleid}_R2_counts.csv \
            --match-filter --allow-empty --qual-filter ${head_cmd}
        """
    }
} else if(params.index_file_type == 'none') {
    process read_counts {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(R1), file(R2) from to_fastq_filters

        output:
            tuple sampleid, file("${sampleid}_R1_.fq.gz"), file("${sampleid}_R2_.fq.gz") into bcop_reads
            tuple file("${sampleid}_R1_counts.csv"), file("${sampleid}_R2_counts.csv") into bcop_counts

        // publishDir "${params.output}/barcodecop/", overwrite: true, mode: 'copy'

        """
        read_counts.py ${head_cmd} ${R1} ${sampleid}_R1_.fq.gz ${sampleid}_R1_counts.csv
        read_counts.py ${head_cmd} ${R2} ${sampleid}_R2_.fq.gz ${sampleid}_R2_counts.csv
        """
    }
}

if(params.containsKey("cutadapt_params")) {
    cutadapt_params_str = "${params.cutadapt_params.join(' ')}"
    process cutadapt {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(R1), file(R2) from bcop_reads
        output:
            tuple sampleid, file("${sampleid}_R1_trimmed.fq.gz"), file("${sampleid}_R2_trimmed.fq.gz") into cutadapt_reads
            file("counts.csv") into cutadapt_counts

        publishDir "${params.output}/cutadapt/${sampleid}/", overwrite: true, mode: 'copy', pattern: '*.{json,tsv,csv}'

        """
        cutadapt ${cutadapt_params_str} -o ${sampleid}_R1_trimmed.fq.gz -p ${sampleid}_R2_trimmed.fq.gz ${R1} ${R2} --json=${sampleid}.cutadapt.json --report=minimal > ${sampleid}.cutadapt.tsv
        echo -n 'sampleid\n${sampleid}' | xsv cat columns --delimiter '\t' --output counts.csv - ${sampleid}.cutadapt.tsv
        """
    }
} else {
    process no_cutadapt {

        label 'med_cpu_mem'

        input:
            tuple sampleid, file(R1), file(R2) from bcop_reads
        output:
            tuple sampleid, file(R1), file(R2) into cutadapt_reads
            file("counts.csv") into cutadapt_counts

        publishDir "${params.output}/cutadapt/${sampleid}/", overwrite: true, mode: 'copy', pattern: '*.{json,tsv,csv}'

        """
        touch counts.csv
        """
    }
}

if (params.alignment.strategy == 'cmsearch') {
  process cm_split {
      cpus '8'
      memory '20 GB'

      input:
          tuple sampleid, file(R1), file(R2) from cutadapt_reads
          file("model.cm") from maybe_local(params.alignment.model)
      output:
          tuple sampleid, val("forward"), file("forward/${R1}"), file("forward/${R2}") into forward_reads
          tuple sampleid, val("reverse"), file("reverse/${R1}"), file("reverse/${R2}") into reverse_reads
          file("counts.csv") into split_counts, split_filter

      publishDir "${params.output}/plus_only/${sampleid}/", overwrite: true, mode: 'copy'

      """
      python3 -c "from Bio import SeqIO;import gzip;SeqIO.write(SeqIO.parse(gzip.open('${R1}', 'rt'), 'fastq'), 'R1.fa', 'fasta')"
      cmsearch -E 10.0 --hmmonly --noali --tblout scores.txt model.cm R1.fa
      split_reads.py --counts counts.csv --cmsearch scores.txt ${sampleid} ${R1} ${R2}
      """
  }
  split_reads = forward_reads.concat(reverse_reads)
} else if (params.alignment.strategy == 'vsearch') {
  process vsearch_split {
      cpus '8'
      memory '20 GB'

      input:
          tuple sampleid, file(R1), file(R2) from cutadapt_reads
          file("library.fna.gz") from maybe_local(params.alignment.library)
      output:
          tuple sampleid, val("forward"), file("forward/${R1}"), file("forward/${R2}") into forward_reads
          tuple sampleid, val("reverse"), file("reverse/${R1}"), file("reverse/${R2}") into reverse_reads
          file("counts.csv") into split_counts, split_filter

      publishDir "${params.output}/split/${sampleid}/", overwrite: true, mode: 'copy'

      // TODO: Switch to fastalite instead of an inline python3 SeqIO to speed things up
      // TODO: split_reads.py --counts [filename,forward,reverse,off_target]
      """
      python3 -c "from Bio import SeqIO;import gzip;SeqIO.write(SeqIO.parse(gzip.open('${R1}', 'rt'), 'fastq'), 'R1.fa', 'fasta')"
      vsearch --usearch_global R1.fa --db library.fna.gz --id 0.75 --query_cov 0.8 --strand both --top_hits_only --userfields query+qstrand --userout hits.tsv
      split_reads.py --counts counts.csv --vsearch hits.tsv ${sampleid} ${R1} ${R2}
      """
  }
  split_reads = forward_reads.concat(reverse_reads)
} else {
  process no_split {
      input:
          tuple sampleid, file(R1), file(R2) from cutadapt_reads
      output:
          tuple sampleid, val("forward"), file(R1), file(R2) into split_reads
          file("counts.csv") into split_counts, split_filter

      publishDir "${params.output}/split/${sampleid}/", overwrite: true, mode: 'copy'

      """
      touch counts.csv
      """
  }
}

// Filter samples below min_read threshold
split_filter
  .splitCsv(header: false)
  .filter{ it[2].toInteger() >= params.min_reads }
  .map{ it -> it[0..1]}  // [sampleid, orientation]
  .join(split_reads, by: [0,1])
  .set{ split_filtered }  // [sampleid, orientation, R1, R2]

process filter_and_trim {

    label 'med_cpu_mem'

    input:
        tuple sampleid, orientation, file(R1), file(R2) from split_filtered
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        tuple sampleid, orientation, file("${sampleid}_R1_filt.fq.gz"), file("${sampleid}_R2_filt.fq.gz") into filtered_trimmed

    // publishDir "${params.output}/filtered/", overwrite: true, mode: 'copy'

    """
    dada2_filter_and_trim.R \
        --infiles ${R1} ${R2} \
        --params dada_params.json \
        --outfiles ${sampleid}_R1_filt.fq.gz ${sampleid}_R2_filt.fq.gz \
    """
}

// add batch number to samples
batches.splitCsv(header: false)
    .cross(filtered_trimmed)
    // [[sampelid, batch], [sampleid, orientation, R1, R2]]
    .map{ it -> [it[0], it[1][1..-1]].flatten() }
    // [sampleid, batch, orientation, R1, R2]
    .into{ to_learn_errors ; to_dereplicate }

process learn_errors {

    label 'med_cpu_mem'

    input:
        // [sampleid, batch, orientation, R1, R2].map { it[batch, orientation, R1, R2] }.groupTuple(by: [batch, orientation])
        tuple batch, orientation, file("R1_*.fastq.gz"), file("R2_*.fastq.gz") from to_learn_errors.map{ it[1, 2, 3, 4] }.groupTuple(by: [0, 1])

    output:
        file("error_model_${batch}_${orientation}.rds") into error_models
        file("error_model_${batch}_${orientation}.png") into error_model_plots

    publishDir "${params.output}/error_models/", overwrite: true, mode: 'copy'

    // non_empty_gz.sh emits filenames to stdout only if uncompressed size != 0,
    // thus dada2_learn_errors.R is provided with a list of non-empty files
    """
    non_empty_gz.sh \$(ls R1_*.fastq.gz) > R1.txt
    non_empty_gz.sh \$(ls R2_*.fastq.gz) > R2.txt
    dada2_learn_errors.R --r1 R1.txt --r2 R2.txt \
        --model error_model_${batch}_${orientation}.rds \
        --plots error_model_${batch}_${orientation}.png
    """
}

process dada_dereplicate {

    label 'med_cpu_mem'

    input:
        tuple sampleid, batch, orientation, file(R1), file(R2) from to_dereplicate
        file("") from error_models.collect()
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        tuple sampleid, orientation, file("dada.rds") into dada_data
        tuple sampleid, orientation, file("seqtab.csv") into dada_seqtab
        file("dada_counts.csv") into dada_counts
        file("dada_overlaps.csv") into dada_overlaps

    publishDir "${params.output}/dada/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

    """
    dada2_dada.R ${R1} ${R2} \
        --errors error_model_${batch}_${orientation}.rds \
        --sampleid ${sampleid} \
        --params dada_params.json \
        --data dada.rds \
        --seqtab seqtab.csv \
        --counts counts.csv \
        --overlaps overlaps.csv
    echo 'orientation\n${orientation}' | xsv cat columns --output dada_counts.csv counts.csv -
    echo 'orientation\n${orientation}' | xsv cat columns --output dada_overlaps.csv overlaps.csv -
    """
}

process dada_get_unmerged {

    label 'med_cpu_mem'

    input:
        tuple sampleid, orientation, file(dada_rds) from dada_data
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        tuple sampleid, orientation, file(dada_rds) into dada_unmerged_samples
        file("unmerged_*.fasta")

    publishDir "${params.output}/dada/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

    """
    get_unmerged.R ${dada_rds} \
        --forward-seqs unmerged_F.fasta --reverse-seqs unmerged_R.fasta
    """
}

process dada_get_dropped_chimeras {

    label 'med_cpu_mem'

    input:
        tuple sampleid, orientation, file(dada_rds) from dada_unmerged_samples
        file("dada_params.json") from maybe_local(params.dada_params)

    output:
        file("chim_dropped.csv")

    publishDir "${params.output}/dada/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

    """
    get_dropped_chim.R ${dada_rds} \
        --outfile chim_dropped.csv
    """
}

process combined_overlaps {

    input:
        file("overlaps_*.csv") from dada_overlaps.collect()

    output:
        file("overlaps.csv")

    publishDir "${params.output}", overwrite: true, mode: 'copy'

    """
    xsv cat rows --output overlaps.csv overlaps_*.csv > overlaps.csv
    """
}

process write_seqs {

    input:
        tuple sampleid, orientation, file("seqtab_*.csv") from dada_seqtab.groupTuple(by: [0,1])

    output:
        tuple sampleid, orientation, file("seqs.fasta") into seqs
        file("specimen_map.csv")
        file("sv_table.csv")
        tuple sampleid, file("sv_table_long.csv") into sv_table_long
        tuple sampleid, orientation, file("weights.csv") into weights

    publishDir "${params.output}/write_seqs/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

    // NOTE: --reverse-complement will put reverse reads into forward orientation
    script:
    if( orientation == 'reverse')
        """
        write_seqs.py seqtab_*.csv \
            --label ${orientation} \
            --reverse-complement \
            --seqs seqs.fasta \
            --specimen-map specimen_map.csv \
            --sv-table sv_table.csv \
            --sv-table-long sv_table_long.csv \
            --weights weights.csv
        """
    else
        """
        write_seqs.py seqtab_*.csv \
            --label ${orientation} \
            --seqs seqs.fasta \
            --specimen-map specimen_map.csv \
            --sv-table sv_table.csv \
            --sv-table-long sv_table_long.csv \
            --weights weights.csv
        """
}

if (params.containsKey('alignment_model') && !params.containsKey('alignment')) {
  params.alignment = [:]
  params.alignment.strategy = 'cmsearch'
  params.alignment.model = params.alignment_model
}

if (params.alignment.strategy == 'cmsearch') {
    process cmsearch {
        /** sequences below `-E 10.0` are not included in the alignment scores
        and will be reported in `filter_svs.py --failing`
        **/
        label 'large_cpu_mem'

        input:
            tuple sampleid, orientation, file("seqs.fasta") from seqs
            file("model.cm") from maybe_local(params.alignment.model)

        output:
            tuple sampleid, orientation, file("seqs.fasta"), file("sv_aln_scores.txt") into seqs_to_filter

        publishDir "${params.output}/cmsearch/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

        script:
        """
        cmsearch -E 10.0 --hmmonly --noali --tblout sv_aln_scores.txt model.cm seqs.fasta
        """
    }
} else if (params.alignment.strategy == 'vsearch') {
    process vsearch {
        /** sequences that do not align are not included in the alignment
        scores and will be reported in `filter_svs.py --failing`
        **/
        label 'large_cpu_mem'

        input:
            tuple sampleid, orientation, file("seqs.fasta") from seqs
            file("library.fna.gz") from maybe_local(params.alignment.library)

        output:
            tuple sampleid, orientation, file("seqs.fasta"), file("sv_aln_scores.txt") into seqs_to_filter

        publishDir "${params.output}/vsearch/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

        script:
        """
        vsearch --usearch_global seqs.fasta --db library.fna.gz --iddef 2 --id 0.6 --query_cov 0.8 \
        --strand plus --top_hits_only --userfields query+id+qstrand --userout sv_aln_scores.txt
        """
    }
} else {
     process no_strategy {
        input:
            tuple sampleid, orientation, file("seqs.fasta") from seqs

        output:
            tuple sampleid, orientation, file("seqs.fasta"), file("sv_aln_scores.txt") into seqs_to_filter

        script:
        """
        touch sv_aln_scores.txt
        """
    }
}

seqs_to_filter.join(weights, by: [0,1]).set{ seqs_to_filter }

if (!['cmsearch', 'vsearch'].contains(params.alignment.strategy)) {
  params.alignment.strategy = 'none'
}

process filter_svs {
    input:
        tuple sampleid, orientation, file("seqs.fasta"), file("sv_aln_scores.txt"), file("weights.csv") from seqs_to_filter
        val strategy from params.alignment.strategy

    output:
        tuple sampleid, file("passed.fasta") into passed
        file("failed.fasta")
        file("outcomes.csv")
        file("counts.csv") into passed_counts

    publishDir "${params.output}/filter_svs/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

    """
    filter_svs.py \
        --counts counts.csv \
        --failing failed.fasta \
        --min-score 0 \
        --outcomes outcomes.csv \
        --passing passed.fasta \
        --strategy ${strategy} \
        --weights weights.csv \
        seqs.fasta ${orientation} sv_aln_scores.txt
    """
}

if(params.containsKey("bidirectional") && params.bidirectional){

    process vsearch_svs {
        // Append size/weight to sequence headers: ";size=integer"
        // so vsearch can sort by weight
        input:
            tuple sampleid, file("seqs_*.fasta"), file("sv_table_long_*.csv") from passed.groupTuple().join(sv_table_long.groupTuple())

        output:
            file("seqs.fasta") into cluster_seqs
            file("sv_table_long.csv") into cluster_sv_table_long
            file("clusters.uc") into clusters

        publishDir "${params.output}/vsearch_svs/${sampleid}/", overwrite: true, mode: 'copy'

        """
        cat seqs_*.fasta > seqs.fasta
        xsv cat rows --output sv_table_long.csv sv_table_long_*.csv
        append_size.py seqs.fasta sv_table_long.csv | \
        vsearch --cluster_size - --uc clusters.uc --id 1.0 --iddef 2 --xsize
        """
    }

    process combine_svs {
        input:
            file("seqs_*.fasta") from cluster_seqs.collect()
            file("sv_table_long_*.csv") from cluster_sv_table_long.collect()
            file("clusters_*.uc") from clusters.collect()

        output:
            file("seqtab.csv") into combined

        publishDir params.output, overwrite: true, mode: 'copy'

        """
        cat seqs_*.fasta > seqs.fasta
        xsv cat rows --output sv_table_long.csv sv_table_long_*.csv
        xsv cat rows --no-headers --output clusters.uc clusters_*.uc
        combine_svs.py --out seqtab.csv seqs.fasta sv_table_long.csv clusters.uc
        """
    }

    process write_complemented_seqs {
        // NOTE: sv names will be regenerated and will
        // not be traceable to earlier steps in the pipeline
        // TODO: create some kind of map back to original seq names
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
        file("raw.csv") from raw_counts
        file("cutadapt_*.csv") from cutadapt_counts.collect()
        file("split_*.csv") from split_counts.collect()
        file("bcop_*.csv") from bcop_counts.collect()
        file("dada_*.csv") from dada_counts.collect()
        file("passed_*.csv") from passed_counts.collect()
    output:
        file("counts.csv")

    publishDir params.output, overwrite: true, mode: 'copy'

    """
    xsv cat rows --output cutadapt.csv  cutadapt_*.csv
    xsv cat rows --no-headers --output split.csv split_*.csv
    xsv cat rows --no-headers --output bcop.csv bcop_*.csv
    xsv cat rows --output dada.csv dada_*.csv
    xsv cat rows --output passed.csv passed_*.csv
    counts.py --out counts.csv raw.csv cutadapt.csv split.csv bcop.csv dada.csv passed.csv
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
