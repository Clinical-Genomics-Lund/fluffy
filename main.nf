#!/usr/bin/dev nextflow

OUTDIR = params.outdir+'/'+params.subdir

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.read1), file(row.read2), row.run)}
    .into {fastqs_bwa;fastqs_qc;cleanup}
    

process fastqc{
    cpus 2
    memory '10GB'
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/QC", mode: 'copy'
    time '10m'

    input:
        set val(sampl_id), file(read1), file(read2), run from fastqs_qc
    
    output:
        set run, file('*') into fatqc_out
     
     script:
	    fastqs = params.singleEnd ? "${read1}" : "${read1} ${read2}"

        """
        fastqc -d ${params.tmp_dir} -o . ${fastqs}
        """
}


process bwa{
    memory '15GB'
    tag "${sampl_id}"
    //publishDir "${OUTDIR}/${sampl_id}", mode: 'copy'
    cpus 16
    time '2h'
     
    
    input:
        set val(sampl_id), file(read1),file(read2), run from fastqs_bwa

    output:
	    set val(sampl_id), file("${sampl_id}.sam"), run into bwa_saisam
        val(run) into run_info

    script:
	
	if (params.singleEnd)
		"""
		bwa aln -n 0 -k 0 ${params.genome_fasta} ${read1} > ${sampl_id}.sai
		bwa samse -n -1 ${params.genome_fasta} ${sampl_id}.sai ${read1} >  ${sampl_id}.sam
    	""" 
	
	else if(params.pairedEnd)
		"""
	    bwa aln -n 0 -k 0 ${params.genome_fasta} ${read1} > ${sampl_id}.R1.sai
	    bwa aln -n 0 -k 0 ${params.genome_fasta} ${read2} > ${sampl_id}.R2.sai
	    bwa sampe -n -1 ${params.genome_fasta} ${sampl_id}.R1.sai ${sampl_id}.R2.sai ${read1} ${read2} > ${sampl_id}.sam	       
	    """
	}


process bamsormadup{
    //convert sam file to bam and sort them based on the coordinate
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    cpus 16
    memory '20GB'
    time '20m'
    
    input:
        set val(sampl_id), file(aln_sam), run from bwa_saisam
        
    output:
        set val(sampl_id), file("${sampl_id}.bam"), file("${sampl_id}.bam.bai"), run into bambai, bam_picard1, bam_picard2, bam_picard3, bam_tiddit

    script:
    """
    bamsormadup \\
        inputformat=sam \\
        threads=${task.cpus} \\
        SO=coordinate \\
        outputformat=bam \\
        tmpfile=${params.tmp_dir}/${sampl_id} \\
        indexfilename=${sampl_id}.bam.bai < ${aln_sam} > ${sampl_id}.bam
    """  
    }

process bam_to_npz{
    //Convert bam files to npz format required by wiscondorx
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory '3GB'
    time '20m'
    cpus 2
    
    input:
        set val(sampl_id), file(bam), file(bai), run from bambai
    
    output:
        set val(sampl_id), file("${sampl_id}.bam.wcx.npz"), run into bam_npz_wcx, bam_npz_preface,bam_npz_wcx_gender
    
    script:
    """
    WisecondorX convert ${bam} ${sampl_id}.bam.wcx.npz
    """
}

// ************* Picard Analysis ********* //

process  CollectGcBiasMetrics{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 16.GB
    time '20m'
    cpus 2


    input:
      set val(sampl_id), file(bam), file(bai), run from bam_picard1  
    
    output:
        set val(sampl_id), file("${sampl_id}.gc_bias_metrics.txt"), file("${sampl_id}.gc.summary.tab") into  picard_gc_info
        file("${sampl_id}.gc_bias_metrics.pdf")
    
    script:
    """
    picard CollectGcBiasMetrics \\
        I=${bam} \\
        O=${sampl_id}.gc_bias_metrics.txt \\
        CHART=${sampl_id}.gc_bias_metrics.pdf \\
        S=${sampl_id}.gc.summary.tab \\
	    R=${params.genome_fasta} \\
        VALIDATION_STRINGENCY=LENIENT -Xms4G -Xmx4G
    """
}


process CollectInsertSizeMetrics{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/picard", mode: 'copy'
    memory 16.GB
    time '20m'
    cpus 2

    input:
        set val(sampl_id), file(bam), file(bai), run from bam_picard2

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


process EstimateLibraryComplexity{
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 16.GB
    time '20m'
    cpus 2

    input:
        set val(sampl_id), file(bam), file(bai), run from bam_picard3
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

// *************** Conerage profiling ****************/
process tiddit{
    //coverage profiling 
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 3.GB
    cpus 1
    time '10m'

    
    input:
        set val(sampl_id), file(bam), file(bai), run from bam_tiddit
    
    output:
        set val(sampl_id), file("${sampl_id}.tiddit.tab"), run into tiddit_cov
    
    script:
    """
    python /bin/TIDDIT.py --cov --bam ${bam} -z 50000 -o ${sampl_id}.tiddit
    """    
}

process generateGeneomeGCtab{
    //Genotype coverage profiling
    //publishDir "${OUTDIR}/${sampl_id}", mode: 'copy' 
    memory 16.GB 
    time '30m'
    cpus 1

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
    // copy number estimation and FFY estimation
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 1.GB
    time '5m'
    cpus 2
    
    input:
        set val(sampl_id), file(tiddit_tab), run from tiddit_cov
        //file(gc_tab) from genomeGCtab
    
    output:
        set val(sampl_id), file("${sampl_id}.tiddit.AMYCNE.tab") into amycne_tab
    
    script:
    """
    python /bin/AMYCNE/AMYCNE.py --ff --coverage ${tiddit_tab} --gc $params.gctab --Q 10 > ${sampl_id}.tiddit.AMYCNE.tab
    """
}
//******************** wisecondorX prediction ***********/
process WCXpredict{
    //wisecondorx prediction 
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 1.GB
    time '5m'
    cpus 2

    input:
        set val(sampl_id), file(bam_npz), run from bam_npz_wcx
    
    output:
        set val(sampl_id), file("${sampl_id}.WCXpredict*") into  wcx_predict

    script:
    """
    WisecondorX \\
        --loglevel info predict ${bam_npz} \\
        ${params.ref500kbp} \\
        ${sampl_id}.WCXpredict \\
        --bed --blacklist ${params.blacklist} \\
        --zscore ${params.zscore}
    """
}

process WCX_preface{
    //wisecondorx prediction based on the Preface ref file(different bin size)
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 3.GB
    cpus 2
    time '5m'

    input:
        set val(sampl_id), file(bam_npz), run from bam_npz_preface
    
    output:
        set val(sampl_id), file("${sampl_id}.WCXpredict.preface*") into  wcx_predict_preface
        set val(sampl_id), file("${sampl_id}.WCXpredict.preface_bins.bed"), run into preface_bins_bed
    
    script:
    """
    WisecondorX \\
    	--loglevel info predict ${bam_npz} \\
        ${params.ref100kbp} \\
        ${sampl_id}.WCXpredict.preface \\
        --bed --blacklist ${params.blacklist} \\
        --zscore ${params.zscore}
    """
}

process WCX_gender{
    //Wisecondorx  gender prediction
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory 3.GB
    cpus 2
    time '5m'
    
    input:
        set val(sampl_id), file(bam_npz), run from bam_npz_wcx_gender
    
    output:
        set val(sampl_id), file("${sampl_id}.wcx.npz.gender.txt") into  wcx_gender
    
    script:
    """
    WisecondorX gender ${bam_npz} ${params.ref500kbp} > ${sampl_id}.wcx.npz.gender.txt
    """
}

process preface{
    //General fetal fraction prediction and ff predicrion based on chX
    tag "${sampl_id}"
    publishDir "${OUTDIR}/${run}/${sampl_id}", mode: 'copy'
    memory '3GB'
    cpus 2
    time '5m'
    
    input:
        set val(sampl_id), file(bins_bed), run from preface_bins_bed
    
    output:
        set val(sampl_id), file("${sampl_id}_bins.bed.PREFACE.txt") into preface_out
    
    script:
    """
    Rscript /bin/PREFACE-0.1.1/PREFACE.R predict \\
        --infile ${bins_bed} \\
        --model ${params.model} > ${sampl_id}_bins.bed.PREFACE.txt
    """
}
/**************** batch summary ***********/       
process multiqc{
    publishDir "${OUTDIR}/${run}", mode: 'copy'
    memory '10GB'
    cpus 1
    time '20m'
    
    input:
        set run, file(qcfiles) from fatqc_out.groupTuple()
    
    output:
        file('multiqc_report.html')
        file('*')

    script:
    qcfiles = qcfiles.join( ' ' )
    """
    export LC_ALL=C.UTF-8                                                    
    export LANG=C.UTF-8                                                      
                          
    multiqc "${OUTDIR}/${run}/QC/"
    """
}

process summary{
    //to do: add  run folder to summary file name 
    publishDir "${OUTDIR}/${run}/summary", mode: 'copy'
    container = '/fs1/resources/containers/wgs_2021-03-16.sif'  // temporary, need numpy not available in nipt-container as of 2021-07-12
    cpus 2
    time '30m'
    memory '3GB'

    input:
        file(preface_out) from preface_out.collect()
        //file(sex) from wcx_gender.collect()
        file(wcx_preface) from wcx_predict_preface.collect()
        file(amycne)  from amycne_tab.collect()
        file(wcx_predict) from wcx_predict.collect()
        val(run) from run_info.first()

    output:
        file('batch_summary.csv')
    
    """
    generate_csv.py \\
        --folder ${OUTDIR}/${run} \\
        --samplesheet ${params.csv} \\
        --Zscore 5 \\
        --minCNV 10000000 > batch_summary.csv
    """
}