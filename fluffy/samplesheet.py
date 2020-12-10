"""Code to handle samplesheets"""

import logging
from pathlib import Path
from typing import Iterator

LOG = logging.getLogger(__name__)


def get_separator(line: str) -> str:
    """Get the separator for file"""
    if " " in line:
        return " "
    if "," in line:
        return ","
    return None


def get_sample_col(line_content: list) -> int:
    """Get the column number that holds the sample name"""
    for column, info in enumerate(line_content):
        if info.lower() == "sample_id":
            return column
    return None


def read_samplesheet(samplesheet: Iterator[str], project_dir: Path) -> Iterator[dict]:
    """Parse a sample sheet and return sample information

    Yields:
        samples(dict): a dictionary with a list of commands and 'se'(bool)
    """

    samples = set()
    sample_col = 0

    header_line=False
    for line_nr, line in enumerate(samplesheet):
    
        if not line.startswith("Sample_ID") and header_line == False:
            continue
        
        if  line.startswith("Sample_ID"):
            
            header_line=True
            separator = get_separator(line)
            LOG.debug("Use separator %s", separator)
            header = line.rstrip().split(separator)
            sample_col = get_sample_col(header)
            continue
        
        if 'NIPT' in line: 
            content = line.rstrip().split(separator)
            sample_name = content[sample_col]
        
            if sample_name in samples:
                continue

            samples.add(sample_name)
        
            single_end = True
            LOG.debug("Check if files are single end or not")
            for file_name in project_dir.glob("*{}*.fastq.gz".format(sample_name)):
                if "_R2" in str(file_name):
                    single_end = False
                    break

            fastq = [
                "<( zcat {}/{}*_R1*fastq.gz )".format(project_dir, sample_name),
                "<( zcat {}/{}*_R2*fastq.gz )".format(project_dir, sample_name),
            ]

        
            if single_end:
                LOG.info("Single end files!")
                fastq = ["<( zcat {}/{}*_R1*fastq.gz )".format(project_dir, sample_name)]
            yield {"fastq": fastq, "single_end": single_end, "sample_id": sample_name}

