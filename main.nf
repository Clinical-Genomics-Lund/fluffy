#!/usr/bin/dev nextflow

OUTDIR = params.outdir+'/'+params.subdir

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.Sample_ID, file(row.read1))}
    .into {fastqs_bwa; fastqs_qc; cleanup}

process fastqc{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/QC", mode: 'copy'

    input:
        set val(sampl_id), file(read1) from fastqs_qc
    
    output:
        file '*' into fatqc_out
    script:
    """
    fastqc -d ${params.tmp_dir} -o . ${read1}
    """
}

process bwa{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    cpus 16

    input:
        set val(sampl_id), file(read1) from fastqs_bwa

    output:
        set val(sampl_id), file(read1), file("${sampl_id}.sai") into bwa_sai

    script:
    """
    bwa aln -n 0 -k 0 ${params.genome_fasta} ${read1} > ${sampl_id}.sai
    """
}

process bamsormadup{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    cpus 16

    input:
        set val(sampl_id), file(read1) , file(sai) from bwa_sai
        
    output:
        set val(sampl_id), file("${sampl_id}.bam"), file("${sampl_id}.bam.bai") into bambai, bam_picard1, bam_picard2, bam_picard3, bam_tiddit

    script:
    """
    bwa samse -n -1 ${params.genome_fasta} ${sai} ${read1} \\
    | \\
    bamsormadup \\
        inputformat=sam \\
        threads=${task.cpus} \\
        SO=coordinate \\
        outputformat=bam \\
        tmpfile=${params.tmp_dir}/${sampl_id} \\
        indexfilename=${sampl_id}.bam.bai > ${sampl_id}.bam
    """  
    }

process bam_to_npz{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    
    input:
        set val(sampl_id), file(bam), file(bai) from bambai
    
    output:
        set val(sampl_id), file("${sampl_id}.bam.wcx.npz") into bam_npz_wcx, bam_npz_preface,bam_npz_wcx_gender
    
    script:
    """
     WisecondorX convert ${bam} ${sampl_id}.bam.wcx.npz

    """
}

process  CollectGcBiasMetrics{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    memory 16.GB

    input:
      set val(sampl_id), file(bam), file(bai) from bam_picard1  
    output:
        set val(sampl_id), file("${sampl_id}.gc_bias_metrics.txt"), file("${sampl_id}.gc.summary.tab") into  picard_gc_info
        file("${sampl_id}.gc_bias_metrics.pdf")
    script:
    """
    picard CollectGcBiasMetrics \\
        I=${bam} \\
        O=${sampl_id}.gc_bias_metrics.txt \\
        CHART=${sampl_id}.gc_bias_metrics.pdf \\
        S=${sampl_id}.gc.summary.tab R=${params.genome_fasta} \\
        VALIDATION_STRINGENCY=LENIENT -Xms4G -Xmx4G
    """
}
/*
process CollectInsertSizeMetrics{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/picard", mode: 'copy'
    memory 16.GB

    input:
        set val(sampl_id), file(bam), file(bai) from bam_picard2

    output:
        set val(sampl_id),file('insert_size_metrics.txt') into picard_insert
        //set val(sampl_id), file("${sampl_id}.insert_size_metrics.txt") into picard_insert_info
        //file("${sampl_id}.insert_size_histogram.pdf")
    
    script:
    """
    picard CollectInsertSizeMetrics \\
        I=${bam} \\
        O=insert_size_metrics.txt \\
        H=insert_size_histogram.pdf \\
        VALIDATION_STRINGENCY=LENIENT M=0.5 -Xms4G -Xmx4G
    
    """
    //mv  insert_size_metrics.txt ${sampl_id}.insert_size_metrics.txt
    //mv insert_size_histogram.pdf ${sampl_id}.insert_size_histogram.pdf
}
*/
process EstimateLibraryComplexity{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    memory 16.GB

    input:
        set val(sampl_id), file(bam), file(bai) from bam_picard3
    output:
        set val(sampl_id), file("${sampl_id}.complex_metrics.txt") into picard_complexity
    
    script:
    """
    picard EstimateLibraryComplexity \\
        I=${bam} \\
        O=${sampl_id}.complex_metrics.txt \\
        VALIDATION_STRINGENCY=LENIENT -Xms4G -Xmx4G
    """    
}

process tiddit{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    memory 16.GB
    
    input:
        set val(sampl_id), file(bam), file(bai) from bam_tiddit
    
    output:
        set val(sampl_id), file("${sampl_id}.tiddit.tab") into tiddit_cov
    
    script:
    """
    python /bin/TIDDIT.py --cov --bam ${bam}  -z 50000 -o ${sampl_id}.tiddit
    """    
}

process generateGeneomeGCtab{
    //publishDir "${OUTDIR}/${sampl_id}", mode: 'copy' 
    memory 16.GB 
    when:
        params.genome_gc
    output:
        file('genome.gc.tab') into genomeGCtab
    script:
    """
    python /bin/AMYCNE/Generate_GC_tab.py \\
    --fa ${params.genome_fasta} \\
    --size 50000 \\
    --n_mask > genome.gc.tab
    """
}

process AMYCNE{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    memory 16.GB
    
    input:
        set val(sampl_id), file(tiddit_tab) from tiddit_cov
        file(gc_tab) from genomeGCtab
    
    output:
        set val(sampl_id), file("${sampl_id}.tiddit.AMYCNE.tab") into amycne_tab
    
    script:
    """
    python /bin/AMYCNE/AMYCNE.py --ff --coverage ${tiddit_tab} --gc ${gc_tab} --Q 10 > ${sampl_id}.tiddit.AMYCNE.tab
    """
}

process WCXpredict{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
     memory 16.GB
    input:
        set val(sampl_id), file(bam_npz) from bam_npz_wcx
    
    output:
        set val(sampl_id), file("${sampl_id}.WCXpredict*") into  wcx_predict

    script:
    """
    WisecondorX \\
        --loglevel info predict ${bam_npz} \\
        ${params.WCX_reftest} \\
        ${sampl_id}.WCXpredict \\
        --bed --blacklist ${params.blacklist} \\
        --zscore ${params.zscore}
    """
}

process WCX_preface{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
     memory 16.GB
    input:
        set val(sampl_id), file(bam_npz) from bam_npz_preface
    
    output:
        set val(sampl_id), file("${sampl_id}.WCXpredict.preface*") into  wcx_predict_preface
        set val(sampl_id), file("${sampl_id}.WCXpredict.preface_bins.bed") into preface_bins_bed
    
    script:
    """
    WisecondorX --loglevel info predict \\
        ${bam_npz} \\
        ${params.refpreface} \\
        ${sampl_id}.WCXpredict.preface \\
        --bed --blacklist ${params.blacklist} \\
        --zscore ${params.zscore}
    """
}

process WCX_gender{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    memory 16.GB
    input:
        set val(sampl_id), file(bam_npz) from bam_npz_wcx_gender
    
    output:
        set val(sampl_id), file("${sampl_id}.wcx.npz.gender.txt") into  wcx_gender
    
    script:
    """
    WisecondorX gender ${bam_npz} ${params.WCX_reftest} > ${sampl_id}.wcx.npz.gender.txt
    """
}

process preface{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    
    input:
        set val(sampl_id), file(bins_bed) from preface_bins_bed
    
    output:
        set val(sampl_id), file("${sampl_id}_bins.bed.PREFACE.txt") into preface_out
    
    script:
    """
    Rscript /bin/PREFACE-0.1.1/PREFACE.R predict \\
        --infile ${bins_bed} \\
        --model ${params.model} > ${sampl_id}_bins.bed.PREFACE.txt
    """
}

process multiqc{
    publishDir "${OUTDIR}", mode: 'copy'

    output:
        file('multiqc_report.html')
        file('*')

    script:
    """
    export LC_ALL=C.UTF-8                                                    
    export LANG=C.UTF-8                                                      
                          
    multiqc ${OUTDIR}/QC/ 
    """
}

process summarizebatch{
    //to do: add  run folder to summary file name 
    publishDir "${OUTDIR}", mode: 'copy'

    output:
        file('batch_summary.csv')
    
    script:
    //parts = r1.toString().split('/')
	//parts.println()
	//idx= parts.findIndexOf {it ==~ /......_......_...._........../}
	//rundir= parts[0..idx].join("/")
    """
    python /fs1/sima/nipt/fluffy/scripts/generate_csv.py \\
        --folder ${OUTDIR} \\
        --samplesheet ${params.csv} \\
        --Zscore 5 \\
        --minCNV 10000000 > batch_summary.csv
    """
}

process cleanup{
    tag "${sampl_id}"

    input:
        set val(sampl_id), file(read1) from cleanup
    
    output:
        file('*')
    
    script:
    /*
    export LC_ALL=C.UTF-8                                                                                                       
        export LANG=C.UTF-8 
        multiqc ${OUTDIR}/${sampl_id} --outdir ${OUTDIR}/${sampl_id}
    */
    """
        cp -r ${OUTDIR}/${sampl_id} ${OUTDIR}/${sampl_id}.fluffy-0.4.0
        rm -f ${sampl_id}.fluffy-0.4.0/*.sai ${OUTDIR}/${sampl_id}.fluffy-0.4.0/*.bam 
        zip -r ${OUTDIR}/${sampl_id}.fluffy-0.4.0.zip ${OUTDIR}/${sampl_id}.fluffy-0.4.0
        rm -rf ${OUTDIR}/${sampl_id}.fluffy-0.4.0
        
    """
    
}