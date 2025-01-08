# -*- coding: utf-8 -*-
"""

This script reads a PyPSA network and builds reference statistics to be used for comparison.

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

        snakemake = mock_snakemake("build_network_statistics")
