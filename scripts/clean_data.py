import pandas as pd
from helpers import (
    three_2_two_digits_country,
    read_csv_nafix,
    to_csv_nafix,
    configure_logging,
    country_name_2_two_digits,
)
import os


def get_demand_our_world_in_data(fp_input, fp_output):
    """
    Retrieve the electricity demand data from Our World in Data
    """
    df = read_csv_nafix(fp_input)
    df = df.loc[:, ["iso_code", "country", "year", "electricity_demand"]]
    df = df[df["iso_code"].notna()]  # removes antartica
    df["iso_code_2"] = df.loc[:, "iso_code"].apply(lambda x: three_2_two_digits_country(x))
    to_csv_nafix(df, fp_output)


def clean_capacity_IRENA(df_irena):
    """
    Clean the capacity data from IRENA
    """
    df = df_irena.copy()

    # Process technologies
    df.loc[
        df["Technology"].isin(["Solar photovoltaic", "Solar thermal energy"]), "Technology"
    ] = "solar"
    df.loc[df["Technology"].isin(["Onshore wind energy"]), "Technology"] = "onshore wind"
    df.loc[df["Technology"].isin(["Offshore wind energy"]), "Technology"] = "offshore wind"
    df.loc[
        df["Technology"].isin(
            ["Renewable hydropower", "Mixed Hydro Plants", "Pumped storage"]
        ),
        "Technology",
    ] = "hydro"
    df.loc[
        df["Technology"].isin(["Other non-renewable energy", "Marine energy"]), "Technology"
    ] = "other"
    df.loc[
        df["Technology"].isin(["Liquid biofuels", "Biogas", "Solid biofuels"]), "Technology"
    ] = "bioenergy"
    df.loc[df["Technology"].isin(["Geothermal energy"]), "Technology"] = "geothermal"
    df.loc[df["Technology"].isin(["Natural gas"]), "Technology"] = "gas"
    df.loc[df["Technology"].isin(["Renewable municipal waste"]), "Technology"] = "waste"
    df.loc[df["Technology"].isin(["Coal and peat"]), "Technology"] = "coal"
    df.loc[df["Technology"].isin(["Oil", "Fossil fuels n.e.s."]), "Technology"] = "oil"


    df["p_nom"] = pd.to_numeric(
        df["Electricity statistics (MW/GWh)"], errors="coerce"
    )
    installed_capacity_irena = (
        df.rename(columns={"Technology": "carrier"})
        .groupby(["alpha2", "carrier"])["p_nom"]
        .sum()
    )
    return installed_capacity_irena


def get_capacity_IRENA(fp_input, fp_output):
    """
    Retrieve the electricity demand data from IRENA
    """
    df_irena = read_csv_nafix(fp_input, skiprows=2, encoding="latin-1")
    df_irena = df_irena.iloc[:, [0, 1, 2, 4]]
    # df = df[df["iso_code"].notna()]  # removes antartica
    df_irena["alpha2"] = df_irena.loc[:, "Country/area"].apply(
        lambda x: country_name_2_two_digits(x)
    )
    df_irena = clean_capacity_IRENA(df_irena)
    to_csv_nafix(df_irena, fp_output)


if __name__ == "__main__":
    if "snakemake" not in globals():
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        from helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_data")
    
    configure_logging(snakemake)

    get_demand_our_world_in_data(
        snakemake.input.demand_owid,
        snakemake.output.demand_owid,
    )

    get_capacity_IRENA(
        snakemake.input.cap_irena,
        snakemake.output.cap_irena,
    )