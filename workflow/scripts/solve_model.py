# SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
#
# SPDX-License-Identifier: GPL-3.0-or-later

import pypsa
import logging

logger = logging.getLogger(__name__)


def add_ghg_constraint(n: pypsa.Network, primary: dict) -> None:
    """Add greenhouse gas constraint (combining CO2 and CH4 using GWP)."""
    if "ghg" in primary:
        logger.info("Adding greenhouse gas constraint...")
        ghg_limit = float(primary["ghg"]["limit"])  # kg CO2-eq

        # Get the linopy model
        m = n.model

        # GWP values: CO2 = 1, CH4 = 25 (100-year GWP)
        ghg_expression = None

        # Add CO2 contribution if CO2 store exists
        if "co2" in n.stores.index:
            co2_store_idx = n.stores.index.get_loc("co2")
            co2_term = (
                m.variables["Store-e"].isel(name=co2_store_idx, snapshot=-1) * 1
            )  # GWP = 1
            ghg_expression = co2_term

        # Add CH4 contribution if CH4 store exists
        if "ch4" in n.stores.index:
            ch4_store_idx = n.stores.index.get_loc("ch4")
            ch4_term = (
                m.variables["Store-e"].isel(name=ch4_store_idx, snapshot=-1) * 25
            )  # GWP = 25
            if ghg_expression is not None:
                ghg_expression = ghg_expression + ch4_term
            else:
                ghg_expression = ch4_term

        # Add constraint if we have any GHG emissions
        if ghg_expression is not None:
            m.add_constraints(ghg_expression <= ghg_limit, name="max_ghg_emissions")
            logger.info("Total GHG limit: %.1f Gt CO2-eq", ghg_limit / 1e9)


def add_ghg_objective(n: pypsa.Network, ghg_price: float) -> None:
    """Add GHG emissions to the objective function."""
    logger.info("Adding GHG emissions to objective function...")

    # Get the linopy model
    m = n.model

    # GWP values: CO2 = 1, CH4 = 25 (100-year GWP)
    ghg_expression = None

    # Add CO2 contribution if CO2 store exists
    if "co2" in n.stores.index:
        co2_store_idx = n.stores.index.get_loc("co2")
        co2_term = (
            m.variables["Store-e"].isel(name=co2_store_idx, snapshot=-1) * 1
        )  # GWP = 1
        ghg_expression = co2_term
        logger.info("Added CO2 emissions to objective")

    # Add CH4 contribution if CH4 store exists
    if "ch4" in n.stores.index:
        ch4_store_idx = n.stores.index.get_loc("ch4")
        ch4_term = (
            m.variables["Store-e"].isel(name=ch4_store_idx, snapshot=-1) * 25
        )  # GWP = 25
        if ghg_expression is not None:
            ghg_expression = ghg_expression + ch4_term
        else:
            ghg_expression = ch4_term
        logger.info("Added CH4 emissions to objective")

    # Add GHG emissions to objective if we have any
    if ghg_expression is not None:
        # Add to objective (minimizing GHG emissions)
        m.objective = m.objective + ghg_price * ghg_expression
        logger.info("GHG weight in objective: %s", ghg_price)


if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.network)

    # Create the linopy model
    n.optimize.create_model()

    # Add GHG constraint if specified
    add_ghg_constraint(n, snakemake.params.primary)

    # Add GHG emissions to objective function
    add_ghg_objective(n, float(snakemake.params.ghg_price))

    # Solve the model with configured solver
    result = n.optimize.solve_model(
        solver_name=snakemake.params.solver,
        solver_options=snakemake.params.solver_options,
    )

    # Check for infeasibility and diagnose if needed
    if result == ("ok", "optimal"):
        # Optimization successful
        n.export_to_netcdf(snakemake.output.network)
    elif result == ("warning", "infeasible") or result == (
        "warning",
        "infeasible_or_unbounded",
    ):
        logger.error("Model is infeasible or unbounded!")
        try:
            # Try to compute and print infeasibilities (Gurobi only)
            if solver_name.lower() == "gurobi":
                logger.error("Infeasible constraints:")
                n.model.print_infeasibilities()
            else:
                logger.error(
                    "Infeasibility diagnosis only available with Gurobi solver"
                )
        except Exception as e:
            logger.error("Could not compute infeasibilities: %s", e)
    else:
        logger.error("Optimization unsuccessful: %s", result)
