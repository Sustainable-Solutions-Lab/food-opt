import pandas as pd
import pypsa
import snakemake


def build_network(config: dict, foods: pd.DataFrame) -> pypsa.Network:
    n = pypsa.Network()
    n.set_snapshots(["now"])

    # Primary resources
    n.add("Carrier", name="land", unit="km^2")
    n.add("Carrier", "water", unit="m^3")
    n.add("Carrier", "fertilizer", unit="kg")
    n.add("Carrier", "co2", unit="kg")
    n.add("Carrier", "ch4", unit="kg")
    # Question: should N2O emissions from the use of fertilizer be included? (In order to easily calculate total GHG emissions)

    # Foods
    for crop in config["crops"]:
        n.add("Carrier", crop, unit="t")

    # Nutritional values
    n.add("Carrier", "carbs", unit="t")
    n.add("Carrier", "protein", unit="t")
    n.add("Carrier", "fat", unit="t")

    # Add buses for each carrier
    for carrier in n.carriers.index:
        n.add("Bus", carrier, carrier=carrier)

    # Add stores for primary resources
    for carrier in ["land", "water", "fertilizer", "co2", "ch4"]:
        n.add("Store", carrier, bus=carrier, carrier=carrier)

    for carrier in ["land", "water", "fertilizer"]:
        n.stores.at["land", "e_initial"] = float(config["primary"][carrier]["limit"])

    return n


if __name__ == "__main__":
    # Read food data
    foods = pd.read_csv(snakemake.input.foods, index_col=["crop", "param"])
    print(foods)

    # Build the network
    n = build_network(snakemake.config, foods)

    n.export_to_netcdf(snakemake.output.network)
