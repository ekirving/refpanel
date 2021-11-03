#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

"""
The human reference genome used by the International Genome Sample Resource (IGSR)

https://www.internationalgenome.org/
"""


rule igsr_reference_grch38:
    """Download the GRCh38 assembly FASTA, and all the associated index files needed to perform the alignments"""
    output:
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.dict",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.amb",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.ann",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.bwt",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.pac",
        "data/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.sa",
    shell:
        "wget "
        " --mirror"
        " --quiet"
        " --no-host-directories"
        " --cut-dirs=5"
        " --directory-prefix=data/reference/GRCh38/"
        " ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
