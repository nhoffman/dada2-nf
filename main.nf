import groovy.json.JsonOutput

def maybe_local(fname, checkIfExists = false){
    // Address the special case of using test files in this project
    // when running in batchman, or more generally, run-from-git.
    if(fname && (file(fname).exists() || fname.startsWith('s3://'))) {
        return file(fname, checkIfExists: checkIfExists)
    } else if (fname) {
        return file("$workflow.projectDir/" + fname)
    } else {
        return null
    }
}

def head(fastq, n){
    if (n == "-1") {
        return fastq
    } else {
        return fastq.splitFastq(by: n, compress: true, file: true, limit: n)[0]
    }
}

process fastq_list {
    input:
        path(manifest)

    output:
        path("fastq_list.txt")

    """
    fastq_list.py --out fastq_list.txt ${manifest}
    """
}

process parse_manifest {
    input:
        path(sample_info)
        path(fq_paths)
        path(fastqs)

    output:
        path("batches.csv")
        path("counts.csv")
        path("sample_index.csv")
        path(sample_info)
        path(fq_paths)

    publishDir "${params.output}/manifest/", overwrite: true, mode: 'copy'

    """
    manifest.py \
        --batches batches.csv \
        --counts counts.csv \
        --index-file-type ${params.index_file_type} \
        --sample-index sample_index.csv \
        ${sample_info} ${fq_paths}
    """
}

process plot_quality {

    label 'med_cpu_mem'

    input:
        tuple val(sampleid), path(R1), path(R2)
        path("dada_params.json")

    output:
        path("${sampleid}.png")

    publishDir "${params.output}/qplots/", overwrite: true, mode: 'copy'

    """
    dada2_plot_quality.R ${R1} ${R2} --params dada_params.json -o ${sampleid}.png
    """
}

process barcodecop_dual {
    label 'med_cpu_mem'

    input:
      tuple val(sampleid), val(direction), path(fastq), path(I1), path(I2)
      val(head)

    output:
      tuple val(sampleid), val(direction), path("${sampleid}_${direction}_.fq.gz")
      path("counts.csv")

    """
    barcodecop --allow-empty --fastq ${fastq} ${head} --match-filter --outfile ${sampleid}_${direction}_.fq.gz --qual-filter --read-counts counts.csv ${I1} ${I2}
    """
}

process barcodecop_single {
    label 'med_cpu_mem'

    input:
      tuple val(sampleid), val(direction), path(fastq), path(I1)
      val(head)

    output:
      tuple val(sampleid), val(direction), path("${sampleid}_${direction}_.fq.gz")
      path("counts.csv")

    publishDir "${params.output}/barcodecop/${sampleid}/${direction}/", overwrite: true, mode: 'copy'

    """
    barcodecop --allow-empty --fastq ${fastq} ${head} --match-filter --outfile ${sampleid}_${direction}_.fq.gz --qual-filter --read-counts counts.csv ${I1}
    """
}

process no_barcodecop {
    label 'med_cpu_mem'

    input:
      tuple val(sampleid), val(direction), path(fastq)
      val(head)

    output:
      tuple val(sampleid), val(direction), path("${sampleid}_${direction}_.fq.gz")
      path("counts.csv")

    """
    read_counts.py ${head} ${fastq} ${sampleid}_${direction}_.fq.gz counts.csv
    """
}

process cutadapt {
    label 'med_cpu_mem'

    input:
        tuple val(sampleid), path("R1.fq.gz"), path("R2.fq.gz")
        val(cutadapt_params_str)
    output:
        tuple val(sampleid), path("${sampleid}_R1_trimmed.fq.gz"), path("${sampleid}_R2_trimmed.fq.gz")
        path("counts.csv")

    publishDir "${params.output}/cutadapt/${sampleid}/", overwrite: true, mode: 'copy', pattern: '*.{json,tsv,csv}'

    """
    cutadapt ${cutadapt_params_str} -o ${sampleid}_R1_trimmed.fq.gz -p ${sampleid}_R2_trimmed.fq.gz R1.fq.gz R2.fq.gz --json=${sampleid}.cutadapt.json --report=minimal > ${sampleid}.cutadapt.tsv
    echo -n 'sampleid\n${sampleid}' | xsv cat columns --delimiter '\t' --output counts.csv - ${sampleid}.cutadapt.tsv
    """
}

process no_cutadapt {
    input:
        tuple val(sampleid), path("R1.fq.gz"), path("R2.fq.gz")
    output:
        tuple val(sampleid), path("${sampleid}_R1.fq.gz"), path("${sampleid}_R2.fq.gz")
        path("counts.csv")

    publishDir "${params.output}/cutadapt/${sampleid}/", overwrite: true, mode: 'copy', pattern: '*.{json,tsv,csv}'

    """
    cp -P R1.fq.gz ${sampleid}_R1.fq.gz
    cp -P R2.fq.gz ${sampleid}_R2.fq.gz
    touch counts.csv
    """
}

process cmsearch_orientations {
    label "c5d_9xlarge"

    input:
        tuple val(sampleid), path(R1), path(R2)
        path(model)
    output:
        tuple val(sampleid), val("forward"), path("forward/${R1}"), path("forward/${R2}")
        tuple val(sampleid), val("reverse"), path("reverse/${R1}"), path("reverse/${R2}")
        tuple val(sampleid), path("off_target/${R1}"), path("off_target/${R2}")
        path("counts.csv")

    publishDir "${params.output}", saveAs: { (it == "off_target/${R1}" | it == "off_target/${R2}") ? "$it" : "split/${sampleid}/$it" }, overwrite: true, mode: 'copy'

    """
    python3 -c "from Bio import SeqIO;import gzip;SeqIO.write(SeqIO.parse(gzip.open('${R1}', 'rt'), 'fastq'), 'R1.fa', 'fasta')"
    cmsearch -E 10.0 --cpu 32 --hmmonly --noali --tblout scores.txt ${model} R1.fa
    split_reads.py --counts counts.csv --cmsearch scores.txt ${sampleid} ${R1} ${R2}
    """
}

process vsearch_orientations {
    label "c5d_9xlarge"

    input:
        tuple val(sampleid), path(R1), path(R2)
        path(library)
    output:
        tuple val(sampleid), val("forward"), path("forward/${R1}"), path("forward/${R2}")
        tuple val(sampleid), val("reverse"), path("reverse/${R1}"), path("reverse/${R2}")
        tuple val(sampleid), path("off_target/${R1}"), path("off_target/${R2}")
        path("counts.csv")

    publishDir "${params.output}", saveAs: { (it == "off_target/${R1}" | it == "off_target/${R2}") ? "$it" : "split/${sampleid}/$it" }, overwrite: true, mode: 'copy'

    """
    python3 -c "from Bio import SeqIO;import gzip;SeqIO.write(SeqIO.parse(gzip.open('${R1}', 'rt'), 'fastq'), 'R1.fa', 'fasta')"
    vsearch --usearch_global R1.fa --db ${library} --id 0.75 --query_cov 0.8 --strand both --threads 32 --top_hits_only --userfields query+qstrand --userout hits.tsv
    split_reads.py --counts counts.csv --vsearch hits.tsv ${sampleid} ${R1} ${R2}
    """
}

process no_split_orientations {
    input:
        tuple val(sampleid), path(R1), path(R2)
    output:
        tuple val(sampleid), val("forward"), path(R1), path(R2)
        path("counts.csv")

    publishDir "${params.output}/split/${sampleid}/", overwrite: true, mode: 'copy'

    """
    touch counts.csv
    """
}

process filter_and_trim {
    label 'med_cpu_mem'

    input:
        tuple val(sampleid), val(orientation), path(R1), path(R2)
        path("dada_params.json")

    output:
        tuple val(sampleid), val(orientation), path("${sampleid}_R1_filt.fq.gz"), path("${sampleid}_R2_filt.fq.gz")

    // publishDir "${params.output}/filtered/", overwrite: true, mode: 'copy'

    """
    dada2_filter_and_trim.R \
        --infiles ${R1} ${R2} \
        --params dada_params.json \
        --outfiles ${sampleid}_R1_filt.fq.gz ${sampleid}_R2_filt.fq.gz \
    """
}

process learn_errors {
    label 'med_cpu_mem'

    input:
        tuple val(sampleids), val(batch), val(orientation), path("R1_*.fastq.gz"), path("R2_*.fastq.gz")

    output:
        tuple val(sampleids), val(batch), val(orientation), path("error_model_${batch}_${orientation}.rds")
        path("error_model_${batch}_${orientation}.png")

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
    // NOTE: sequences in reverse orientation are reverse complemented to forward orientation for clustering
    label 'med_cpu_mem'

    input:
        tuple val(sampleid), val(batch), val(orientation), path(R1), path(R2), path(model)
        path(dada_params)

    output:
        tuple val(sampleid), val("merged"), path("seqtab.csv")
        tuple val(sampleid), val("R1"), path("seqtab_r1.csv")
        tuple val(sampleid), val("R2"), path("seqtab_r2.csv")
        path("counts.csv")
        path("overlaps.csv")
        path("dada.rds")

    publishDir "${params.output}/dada/${sampleid}/${orientation}/", overwrite: true, mode: 'copy'

    """
    dada2_dada.R ${R1} ${R2} \
        --errors ${model} \
        --sampleid ${sampleid} \
        --orientation ${orientation} \
        --params ${dada_params} \
        --data dada.rds \
        --seqtab seqtab.csv \
        --seqtab-r1 seqtab_r1.csv \
        --seqtab-r2 seqtab_r2.csv \
        --counts counts.csv \
        --overlaps overlaps.csv
    """
}

process combined_overlaps {

    input:
        path("overlaps_*.csv")

    output:
        path("overlaps.csv")

    publishDir "${params.output}", overwrite: true, mode: 'copy'

    """
    xsv cat rows --output overlaps.csv overlaps_*.csv > overlaps.csv
    """
}

process cluster_svs {
    // Convert seqtab.csv into fasta file with headers: "N;specimen=str;size=N"
    // vsearch will use ;size=N to sort by weight
    label "c5d_9xlarge"

    input:
        tuple val(sampleid), val(direction), path("seqtabs_*.csv")

    output:
        tuple val(direction), path("clusters.uc"), path("seqs.fa")

    publishDir "${params.output}/vsearch_svs/${sampleid}/${direction}/", overwrite: true, mode: 'copy'

    """
    fasta.py --out seqs.fa seqtabs_*.csv
    vsearch --cluster_size seqs.fa --uc clusters.uc --id 1.0 --iddef 0 --xsize
    """
}

process combine_svs {
    // Sequence files are already clustered by sampleid and
    // direction so it is safe to collect and combine here
    input:
        tuple val(direction), path("clusters_*.uc"), path("seqs_*.fa")

    output:
        tuple val(direction), path("seqtab.csv")

    // save merged seqtab to base output dir
    publishDir "${params.output}", saveAs: { "${direction}" == "merged" ? "${it}" : "${direction}/${it}" }, overwrite: true, mode: 'copy'

    """
    cat seqs_*.fa > seqs.fa
    xsv cat rows --no-headers --output clusters.uc clusters_*.uc
    combine_svs.py --out seqtab.csv clusters.uc seqs.fa
    """
}

process write_seqs {
    input:
        tuple val(direction), path("seqtab_*.csv")
    output:
        path("specimen_table.csv")
        path("seqs.fasta")
        path("specimen_map.csv")
        path("sv_table.csv")
        path("sv_table_long.csv")
        path("weights.csv")

    // save merged files to base output dir
    publishDir "${params.output}", saveAs: { "${direction}" == "merged" ? "${it}" : "${direction}/${it}" }, overwrite: true, mode: 'copy'

    """
    write_seqs.py \
        --direction ${direction} \
        --seqs seqs.fasta \
        --specimen-map specimen_map.csv \
        --specimen-table specimen_table.csv \
        --sv-table sv_table.csv \
        --sv-table-long sv_table_long.csv \
        --weights weights.csv \
        seqtab_*.csv
    """
}

process join_counts {
    input:
        path("raw.csv")
        path("cutadapt_*.csv")
        path("split_*.csv")
        path("bcop_*.csv")
        path("dada_*.csv")
        path("specimen_counts_*.csv")
        val(downsample)
    output:
        path("counts.csv")

    publishDir params.output, overwrite: true, mode: 'copy'

    """
    xsv cat rows --output cutadapt.csv  cutadapt_*.csv
    xsv cat rows --no-headers --output split.csv split_*.csv
    xsv cat rows --no-headers --output bcop.csv bcop_*.csv
    xsv cat rows --output dada.csv dada_*.csv
    xsv cat rows --no-headers --output specimens.csv specimen_counts_*.csv
    counts.py --out counts.csv raw.csv ${downsample} cutadapt.csv split.csv bcop.csv dada.csv specimens.csv
    """
}

process save_params {

    input:
        val(parameters)

    output:
        path('params.json')

    publishDir params.output, overwrite: true, mode: 'copy'

    """
cat <<EOF > params.json
${parameters}
EOF
    """
}

workflow {
    if(!(params.sample_information && params.fastq_list || params.manifest)){
        println "'sample_information' or 'fastq_list' and manifest is undefined"
        println "provide parameters using '--params-file params.json'"
        System.exit(1)
        }

    dada_params = maybe_local(params.dada_params)

    if (params.containsKey("manifest") && params.manifest) {
        manifest = maybe_local(params.manifest)
        (fq_paths, _) = fastq_list(manifest)
    } else {
        manifest = maybe_local(params.sample_information)
        fq_paths = Channel.fromPath(maybe_local(params.fastq_list))
    }

    fq_files = fq_paths.splitText().map{it.strip()}.map{maybe_local(it, true)}

    // create raw counts and check for sample_info and fastq_list consistency
    (batches, raw_counts, samples, _, _) = parse_manifest(
        manifest,
        fq_paths,  // full sample paths
        fq_files.collect()  // for counts
        )

    samples = samples
        .splitCsv(header: false)
        .map{ it -> [
            it[0],  // sampleid
            it[1],  // direction
            maybe_local(it[2]), // fastq
            maybe_local(it[3]), // I1 (may be null)
            maybe_local(it[4])] } // I2 (may be null)

    if (params.containsKey("downsample") && params.downsample) {
        head_cmd = "--head " + params.downsample
        downsample = params.downsample
    } else {
        head_cmd = ""
        downsample = "-1"
    }

    quality_check = samples
       .map{ it -> [it[0], head(it[2], downsample)] }
       .groupTuple()
       .map{ it -> it.flatten() }
    plot_quality(quality_check, dada_params)

    // Drop empty fastqs. Note: empty fastqs still processed in above steps
    samples = samples.filter{ it[2].countFastq() > 0 }

    if (params.index_file_type == "dual") {
        (filtered, bcop_counts) = barcodecop_dual(samples, head_cmd)
    } else if (params.index_file_type == "single") {
        samples = samples.map{ it -> it[0..3] } // drop I2
        (filtered, bcop_counts) = barcodecop_single(samples, head_cmd)
    } else {
        samples = samples.map{ it -> it[0..2] } // drop I1 and I2
        (filtered, bcop_counts) = no_barcodecop(samples, head_cmd)
    }

    // group R1/R2 reads by sample_id and sort tuples by val(R1/R2) ->
    // [sample_id, [val(R1), R1], [val(R2), R2]]
    filtered = filtered.map{ it -> [it[0], it[1..2]] }.groupTuple(sort: { it[0] })
    // flatten and drop val(R1/R2) directions strings
    filtered = filtered.map{ it -> it.flatten() }.map{ it -> [it[0], it[2], it[4]] }

    if (params.containsKey("cutadapt_params")) {
        cutadapt_params_str = params.cutadapt_params.join(' ')
        (trimmed, cutadapt_counts)  = cutadapt(filtered, cutadapt_params_str)
    } else {
        (trimmed, cutadapt_counts)  = no_cutadapt(filtered)
    }

    trimmed = trimmed.filter{ it[1].countFastq() > 0 }

    if (params.alignment.strategy == "cmsearch") {
        model = maybe_local(params.alignment.model)
        (fwd, rev, _, orientation_counts) = cmsearch_orientations(trimmed, model)
        split = fwd.concat(rev)
    } else if (params.alignment.strategy == "vsearch") {
        library = maybe_local(params.alignment.library)
        (fwd, rev, _, orientation_counts) = vsearch_orientations(trimmed, library)
        split = fwd.concat(rev)
    } else {
        (split, orientation_counts) = no_split_orientations(trimmed)
    }

    // Filter samples below min_read threshold
    split = orientation_counts
        .splitCsv(header: false)
        .filter{ it[2].toInteger() >= params.min_reads }
        .map{ it -> it[0..1]}  // [sampleid, orientation]
        .join(split, by: [0,1])

    filtered = filter_and_trim(split, dada_params)

    // add batch number to samples
    filtered = batches.splitCsv(header: false)
        // [sampelid, batch].cross([sampleid, orientation, R1, R2]])
        .cross(filtered)
        // [sampleid, batch, orientation, R1, R2]
        .map{ it -> [it[0], it[1][1..-1]].flatten() }

    // squash sampleids into list and generate models by batch and orientation
    (models, _) = learn_errors(filtered.groupTuple(by: [1, 2]))  // by: [batch, orientation]
    // transpose/expand out sampleids and join models back into filtered channel
    filtered = filtered.join(models.transpose(), by: [0, 1, 2]) // by: [sampleid, batch, orientation]
    (merged, r1, r2, dada_counts, overlaps, _) = dada_dereplicate(filtered, dada_params)
    combined_overlaps(overlaps.collect())
    seqtabs = merged.concat(r1, r2)

    if (params.containsKey("bidirectional") && params.bidirectional) {
        clusters = cluster_svs(seqtabs.groupTuple(by: [0, 1]))
        seqtabs = combine_svs(clusters.groupTuple())
    } else {
        seqtabs = seqtabs.map{ it -> it[1..-1] }
    }

    (specimen_counts, _) = write_seqs(seqtabs.groupTuple())

    join_counts(
        raw_counts,
        cutadapt_counts.collect(),
        orientation_counts.collect(),
        bcop_counts.collect(),
        dada_counts.collect(),
        specimen_counts.collect(),
        downsample)
    save_params(JsonOutput.prettyPrint(JsonOutput.toJson(params)))
}
