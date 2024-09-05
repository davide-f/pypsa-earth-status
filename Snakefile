# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir
from shutil import copyfile, move

from helpers import create_country_list


rule clean_data:
    input:
        demand_owid="data/owid-energy-data.csv", # from https://nyc3.digitaloceanspaces.com/owid-public/data/energy/owid-energy-data.csv
        demand_iea="data/WEO2023_AnnexA_Free_Dataset_Regions.csv", # from https://www.iea.org/data-and-statistics/data-product/world-energy-outlook-2023-free-dataset-2
        cap_irena="data/ELECSTAT_20240808-144258.csv", # IRENA capacity data from https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Power%20Capacity%20and%20Generation/Country_ELECSTAT_2024_H2.px/
    output:
        demand_owid="resources/owid_demand_data.csv",
        cap_irena="resources/irena_capacity_data.csv",
    script:
        "scripts/clean_data.py"

