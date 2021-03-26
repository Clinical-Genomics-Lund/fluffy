#Run the nextflow pipeline
module load Java singularity nextflow/19.10.0 
nextflow run main.nf -c /fs1/sima/nipt/nextflow.config --csv /fs1/seqdata/NovaSeq/210310_A00681_0318_BH3VLVDRXY/pipeline2/NIPT-wgs-nipt.csv -with-singularity /fs1/sima/nipt/container/FluFFyPipe.sif -resume 

