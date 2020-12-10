module load Java singularity nextflow
nextflow run main.nf -c /fs1/sima/nipt/nextflow.config --csv /fs1/sima/nipt/samplesheet/samples.csv -with-singularity /fs1/sima/nipt/container/FluFFyPipe.sif -resume

