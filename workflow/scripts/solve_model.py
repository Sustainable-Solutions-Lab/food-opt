import pypsa


def add_nutrition_constraints(n: pypsa.Network, config: dict) -> None:
    """Add custom constraints for nutrition requirements based on population."""
    print("Adding nutrition constraints based on population requirements...")

    # Get the linopy model
    m = n.model

    # Get population and calculate total requirements
    population = config.get("population", 1000000)
    days_per_year = 365

    # Add macronutrient constraints
    if "macronutrients" in config:
        print("  Adding macronutrient constraints...")
        for nutrient, nutrient_config in config["macronutrients"].items():
            if (
                "min_per_person_per_day" in nutrient_config
                and nutrient in n.stores.index
            ):
                # Calculate total annual requirement
                min_per_person_per_day = float(
                    nutrient_config["min_per_person_per_day"]
                )  # g/person/day
                total_annual_requirement = (
                    min_per_person_per_day
                    * population
                    * days_per_year
                    / 1000000  # Convert g to tonnes
                )

                # Get the store energy variable for the last snapshot
                store_idx = n.stores.index.get_loc(nutrient)

                # Add constraint: e[store, last_snapshot] >= total_annual_requirement
                # Use isel to select by integer position instead of label
                m.add_constraints(
                    m.variables["Store-e"].isel(Store=store_idx, snapshot=-1)
                    >= total_annual_requirement,
                    name=f"min_{nutrient}_requirement",
                )
                print(
                    f"    {nutrient}: {total_annual_requirement:.1f} t/year ({min_per_person_per_day}g/person/day)"
                )

    # Add food group constraints
    if "food_groups" in config:
        print("  Adding food group constraints...")
        for group_name, group_config in config["food_groups"].items():
            if (
                "min_per_person_per_day" in group_config
                and group_name in n.stores.index
            ):
                # Calculate total annual requirement
                min_per_person_per_day = float(
                    group_config["min_per_person_per_day"]
                )  # g/person/day
                total_annual_requirement = (
                    min_per_person_per_day
                    * population
                    * days_per_year
                    / 1000000  # Convert g to tonnes
                )

                # Get the store energy variable for the last snapshot
                store_idx = n.stores.index.get_loc(group_name)

                # Add constraint: e[store, last_snapshot] >= total_annual_requirement
                # Use isel to select by integer position instead of label
                m.add_constraints(
                    m.variables["Store-e"].isel(Store=store_idx, snapshot=-1)
                    >= total_annual_requirement,
                    name=f"min_{group_name}_requirement",
                )
                print(
                    f"    {group_name}: {total_annual_requirement:.1f} t/year ({min_per_person_per_day}g/person/day)"
                )

    # Add greenhouse gas constraint (combining CO2 and CH4 using GWP)
    if "ghg" in config.get("primary", {}):
        print("  Adding greenhouse gas constraint...")
        ghg_limit = float(config["primary"]["ghg"]["limit"])  # kg CO2-eq

        # GWP values: CO2 = 1, CH4 = 25 (100-year GWP)
        ghg_expression = None

        # Add CO2 contribution if CO2 store exists
        if "co2" in n.stores.index:
            co2_store_idx = n.stores.index.get_loc("co2")
            co2_term = (
                m.variables["Store-e"].isel(Store=co2_store_idx, snapshot=-1) * 1
            )  # GWP = 1
            ghg_expression = co2_term

        # Add CH4 contribution if CH4 store exists
        if "ch4" in n.stores.index:
            ch4_store_idx = n.stores.index.get_loc("ch4")
            ch4_term = (
                m.variables["Store-e"].isel(Store=ch4_store_idx, snapshot=-1) * 25
            )  # GWP = 25
            if ghg_expression is not None:
                ghg_expression = ghg_expression + ch4_term
            else:
                ghg_expression = ch4_term

        # Add constraint if we have any GHG emissions
        if ghg_expression is not None:
            m.add_constraints(ghg_expression <= ghg_limit, name="max_ghg_emissions")
            print(f"    Total GHG limit: {ghg_limit / 1e9:.1f} Gt CO2-eq")


def add_ghg_objective(n: pypsa.Network, config: dict) -> None:
    """Add GHG emissions to the objective function."""
    print("Adding GHG emissions to objective function...")

    # Get the linopy model
    m = n.model

    # GWP values: CO2 = 1, CH4 = 25 (100-year GWP)
    ghg_expression = None

    # Add CO2 contribution if CO2 store exists
    if "co2" in n.stores.index:
        co2_store_idx = n.stores.index.get_loc("co2")
        co2_term = (
            m.variables["Store-e"].isel(Store=co2_store_idx, snapshot=-1) * 1
        )  # GWP = 1
        ghg_expression = co2_term
        print("  Added CO2 emissions to objective")

    # Add CH4 contribution if CH4 store exists
    if "ch4" in n.stores.index:
        ch4_store_idx = n.stores.index.get_loc("ch4")
        ch4_term = (
            m.variables["Store-e"].isel(Store=ch4_store_idx, snapshot=-1) * 25
        )  # GWP = 25
        if ghg_expression is not None:
            ghg_expression = ghg_expression + ch4_term
        else:
            ghg_expression = ch4_term
        print("  Added CH4 emissions to objective")

    # Add GHG emissions to objective if we have any
    if ghg_expression is not None:
        # Get GHG weight from config (default to 1.0 if not specified)
        ghg_weight = config.get("ghg_weight", 1.0)

        # Add to objective (minimizing GHG emissions)
        m.objective = m.objective + ghg_weight * ghg_expression
        print(f"  GHG weight in objective: {ghg_weight}")


def solve_network(n: pypsa.Network) -> pypsa.Network:
    """Solve the food systems optimization problem."""
    print("\nSolving network...")
    print(f"Network has {len(n.links)} links and {len(n.stores)} stores")

    # Create the linopy model
    n.optimize.create_model()

    # Add custom nutrition constraints
    add_nutrition_constraints(n, snakemake.config)

    # Add GHG emissions to objective function
    add_ghg_objective(n, snakemake.config)

    # Configure solver options from config
    solver_name = snakemake.config.get("solving", {}).get("solver", "highs")
    solver_options = snakemake.config.get("solving", {}).get("solver_options", {})

    # Solve the model with configured solver
    print(f"Using solver: {solver_name}")
    result = n.optimize.solve_model(
        solver_name=solver_name, solver_options=solver_options
    )

    # Check for infeasibility and diagnose if needed
    if result == ("warning", "infeasible") or result == (
        "warning",
        "infeasible_or_unbounded",
    ):
        print("Model is infeasible or unbounded! Computing infeasibility diagnosis...")
        try:
            # Try to compute and print infeasibilities (Gurobi only)
            if solver_name.lower() == "gurobi":
                print("Infeasible constraints:")
                n.model.print_infeasibilities()
                infeasible_constraints = n.model.compute_infeasibilities()
                print(
                    f"Number of infeasible constraints: {len(infeasible_constraints)}"
                )
            else:
                print("Infeasibility diagnosis only available with Gurobi solver")
                print(
                    "Consider switching to Gurobi in config for detailed infeasibility analysis"
                )
        except Exception as e:
            print(f"Could not compute infeasibilities: {e}")
    else:
        print(f"Solver result: {result}")

    print("\nOptimization results:")
    if hasattr(n, "objective"):
        print(f"Objective value: {n.objective}")

    # Print some key results
    if len(n.links_t.p0) > 0:
        print("\nLink flows:")
        active_links = n.links_t.p0.loc["now"].abs() > 1e-6
        for link in active_links[active_links].index:
            flow = n.links_t.p0.loc["now", link]
            print(f"  {link}: {flow:.3f}")

    return n


if __name__ == "__main__":
    n = pypsa.Network(snakemake.input.network)
    n = solve_network(n)
    n.export_to_netcdf(snakemake.output.network)
