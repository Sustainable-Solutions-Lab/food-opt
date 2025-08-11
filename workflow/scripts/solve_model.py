import pandas as pd
import pypsa
import snakemake


def solve_network(n: pypsa.Network) -> pypsa.Network:
    # Build objective function: minimize total amount of co2 stored
    m = n.optimize.create_model()
    e = m.variables["Store-e"]
    m.objective = e["co2"].sum()

    n.optimize.solve_model()

    return n


if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.network)
    n = solve_network(n)
    n.export_to_netcdf(snakemake.output.network)
