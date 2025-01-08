# -*- coding: utf-8 -*-
"""

This script cleans raw statistics data from different sources, to build statistics for validation.

"""

import os

import country_converter as coco
import pandas as pd
from helpers import (
    configure_logging,
    country_name_2_two_digits,
    read_csv_nafix,
    three_2_two_digits_country,
    to_csv_nafix,
)

cc = coco.CountryConverter()


def get_demand_ourworldindata(inputs, outputs):
    """
    Retrieve the electricity demand data from Our World in Data
    """
    fp_input = inputs["demand_owid"]
    fp_output = outputs["demand_owid"]
    df = read_csv_nafix(fp_input)
    df = df.loc[:, ["iso_code", "country", "year", "electricity_demand"]]
    df = df[df["iso_code"].notna()]  # removes antartica
    df["iso_code_2"] = cc.pandas_convert(df["iso_code"], to="ISO2")
    to_csv_nafix(df, fp_output)


def clean_capacity_IRENA(df_irena):
    """
    Clean the capacity data from IRENA
    """
    df = df_irena.copy()

    # Process technologies
    df.loc[
        df["Technology"].isin(["Solar photovoltaic", "Solar thermal energy"]),
        "Technology",
    ] = "solar"
    df.loc[df["Technology"].isin(["Onshore wind energy"]), "Technology"] = (
        "onshore wind"
    )
    df.loc[df["Technology"].isin(["Offshore wind energy"]), "Technology"] = (
        "offshore wind"
    )
    df.loc[
        df["Technology"].isin(
            ["Renewable hydropower", "Mixed Hydro Plants", "Pumped storage"]
        ),
        "Technology",
    ] = "hydro"
    df.loc[
        df["Technology"].isin(["Other non-renewable energy", "Marine energy"]),
        "Technology",
    ] = "other"
    df.loc[
        df["Technology"].isin(["Liquid biofuels", "Biogas", "Solid biofuels"]),
        "Technology",
    ] = "bioenergy"
    df.loc[df["Technology"].isin(["Geothermal energy"]), "Technology"] = "geothermal"
    df.loc[df["Technology"].isin(["Natural gas"]), "Technology"] = "gas"
    df.loc[df["Technology"].isin(["Renewable municipal waste"]), "Technology"] = "waste"
    df.loc[df["Technology"].isin(["Coal and peat"]), "Technology"] = "coal"
    df.loc[df["Technology"].isin(["Oil", "Fossil fuels n.e.s."]), "Technology"] = "oil"

    df["p_nom"] = pd.to_numeric(df["Electricity statistics (MW/GWh)"], errors="coerce")
    installed_capacity_irena = (
        df.rename(columns={"Technology": "carrier"})
        .groupby(["alpha2", "carrier"])["p_nom"]
        .sum()
    )
    return installed_capacity_irena


def get_installed_capacity_irena(inputs, outputs):
    """
    Retrieve the electricity demand data from IRENA
    """
    fp_input = inputs["cap_irena"]
    fp_output = outputs["cap_irena"]
    df_irena = read_csv_nafix(fp_input, skiprows=2, encoding="latin-1")
    df_irena = df_irena.iloc[:, [0, 1, 2, 5]]
    # df = df[df["iso_code"].notna()]  # removes antartica
    df_irena["alpha2"] = cc.pandas_convert(df_irena["Country/area"], to="ISO2")
    df_irena = clean_capacity_IRENA(df_irena)
    to_csv_nafix(df_irena, fp_output)


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("clean_data")

    configure_logging(snakemake)

    get_demand_ourworldindata(snakemake.input, snakemake.output)

    get_installed_capacity_irena(snakemake.input, snakemake.output)
