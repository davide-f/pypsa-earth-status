# SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors
#
# SPDX-License-Identifier: AGPL-3.0-or-later

import sys

sys.path.append("./scripts")

from os.path import normpath, exists, isdir
from shutil import copyfile, move

from helpers import create_country_list


configfile: "config.yaml"


rule clean_data:
    params:
        datasets=config["datasets"],
    input:
        demand_owid="data/owid-energy-data.csv",  # from https://nyc3.digitaloceanspaces.com/owid-public/data/energy/owid-energy-data.csv
        demand_iea="data/WEO2023_AnnexA_Free_Dataset_Regions.csv",  # from https://www.iea.org/data-and-statistics/data-product/world-energy-outlook-2023-free-dataset-2
        cap_irena="data/ELECSTAT_20240808-144258.csv",  # IRENA capacity data from https://pxweb.irena.org/pxweb/en/IRENASTAT/IRENASTAT__Power%20Capacity%20and%20Generation/Country_ELECSTAT_2024_H2.px/
        # other sources
    output:
        demand_owid="resources/clean/owid_demand_data.csv",
        cap_irena="resources/clean/irena_capacity_data.csv",
    script:
        "scripts/clean_data.py"


rule build_reference_statistics:
    params:
        datasets=config["datasets"],
    input:
        demand_owid="resources/clean/owid_demand_data.csv",
        cap_irena="resources/clean/irena_capacity_data.csv",
        # other sources
    output:
        demand="resources/reference_statistics/demand.csv",
        installed_capacity="resources/reference_statistics/installed_capacity.csv",
        # energy_dispatch="resources/reference_statistics/energy_dispatch.geojson"
        # network="resources/reference_statistics/network.geojson"
    script:
        "scripts/build_reference_statistics.py"


rule build_network_statistics:
    params:
        network=config["network_validation"],
    input:
    output:
        demand="resources/network_statistics/demand.csv",
        installed_capacity="resources/network_statistics/installed_capacity.csv",
        optimal_capacity="resources/network_statistics/optimal_capacity.csv",
        # energy_dispatch="resources/network_statistics/energy_dispatch.csv",
        # network="resources/network_statistics/network.geojson",
    script:
        "scripts/build_network_statistics.py"


rule make_comparison:
    input:
        demand_network="resources/network_statistics/demand.csv",
        installed_capacity_network="resources/network_statistics/installed_capacity.csv",
        optimal_capacity_network="resources/network_statistics/optimal_capacity.csv",
        # energy_dispatch_network="resources/network_statistics/energy_dispatch.csv",
        # network_network="resources/network_statistics/network.geojson",
        demand_reference="resources/reference_statistics/demand.csv",
        installed_capacity_reference="resources/reference_statistics/installed_capacity.csv",
        # energy_dispatch_reference="resources/reference_statistics/energy_dispatch.geojson"
        # network_reference="resources/reference_statistics/network.geojson"
    output:
        demand_comparison="results/tables/demand.csv",
        installed_capacity_comparison="results/tables/installed_capacity.csv",
        optimal_capacity_comparison="results/tables/optimal_capacity.csv",
        # energy_dispatch_comparison="results/tables/energy_dispatch.geojson"
        # network_comparison="results/tables/network.geojson"
    script:
        "scripts/make_comparison.py"


rule visualize_data:
    input:
        demand_comparison="results/tables/demand.csv",
        installed_capacity_comparison="results/tables/installed_capacity.csv",
        optimal_capacity_comparison="results/tables/optimal_capacity.csv",
        # energy_dispatch_comparison="results/tables/energy_dispatch.geojson"
        # network_comparison="results/tables/network.geojson"
    output:
        plot_demand="results/figures/demand_comparison.png",
        plot_installed_capacity="results/figures/installed_capacity_comparison.png",
        plot_capacity_mix="results/figures/capacity_mix_comparison.png",
    script:
        "scripts/visualize_data.py"