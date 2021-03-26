#If runing the python wrap:
#module  load  Java  singularity
fluffy  --sample /path/to/samplsheet.csv --project /path/to/sequencing_run_dir  --out /path/to/outpout_dir  --config  /path/to/Config.cmd.json  analyse
#You have to delete the output_dir if you intend to rerun the analysis
#if you need to generate the refernce file for wisecondox:
#fluffy --sample <samplesheet>  --project <input_folder> --out <output_folder> reference
#If you want to skip preface in the analysis add --skip_preface:
#fluffy --sample <samplesheet>  --project <input_folder> --out <output_folder> --skip_preface analyse

 
