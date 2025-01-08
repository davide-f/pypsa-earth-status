# -*- coding: utf-8 -*-
"""

This script merges the reference and statistics to generate desired comparison tables

"""
import os
import shutil

import pandas as pd
from helpers import (
    configure_logging,
    country_name_2_two_digits,
    read_csv_nafix,
    three_2_two_digits_country,
    to_csv_nafix,
)

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("make_comparison")

# Relevant for support IPCC AR6 database data format
# https://github.com/martavp/pypsa-eur-sec-to-ipcc/blob/main/pypsa_to_IPCC.py
