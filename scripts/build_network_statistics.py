# -*- coding: utf-8 -*-
"""
This script reads a PyPSA network and builds reference statistics to be used for comparison.
"""
import os
import pandas as pd
import pypsa
import re
from helpers import (
    configure_logging,
    to_csv_nafix,
)

def process_network_statistics(inputs, outputs):
    """
    Extracts and processes demand, installed capacity, and optimal capacity from the PyPSA network.
    """
    network = pypsa.Network(inputs["network_path"])
    
    # Extract demand
    network.loads['country'] = network.buses.loc[network.loads['bus'], 'country']
    demand = network.loads_t.p.sum().T.groupby(network.loads['country']).sum().mul((int(re.search(r'\d+', inputs["opts"][0]).group()))*1e-6)
    demand = demand.reset_index()
    demand.columns = ['iso_code_2', 'demand']
    demand = demand.rename(columns={"iso_code_2":"country"}).set_index("country")
    to_csv_nafix(demand, outputs["demand"])
    
    # Extract installed capacity
    installed_capacity = network.generators[["carrier","p_nom"]]
    installed_capacity.index = network.buses.loc[network.generators['bus'], 'country']
    installed_capacity = installed_capacity.groupby(["country", "carrier"]).sum()
    to_csv_nafix(installed_capacity, outputs["installed_capacity"])
    
    # Extract optimal capacity
    optimal_capacity = network.generators[["carrier","p_nom_opt"]]
    optimal_capacity.index = network.buses.loc[network.generators['bus'], 'country']
    optimal_capacity = optimal_capacity.rename(columns={"p_nom_opt": "p_nom"})
    optimal_capacity = optimal_capacity.groupby(["country", "carrier"]).sum()
    to_csv_nafix(optimal_capacity, outputs["optimal_capacity"])

if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake
        snakemake = mock_snakemake("build_network_statistics")
    
    configure_logging(snakemake)
    
    process_network_statistics(snakemake.params["network"], snakemake.output)