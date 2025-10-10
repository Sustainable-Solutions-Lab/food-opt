.. SPDX-FileCopyrightText: 2025 Koen van Greevenbroek
..
.. SPDX-License-Identifier: CC-BY-4.0

Food Processing & Trade
========================

Food Processing
---------------

Overview
~~~~~~~~

The food processing module converts raw agricultural products (crops and animal products) into final food products consumed by the population. This captures:

* **Mass losses**: Processing inefficiencies (e.g., wheat → flour loses bran/germ)
* **Unit conversion**: Conversion from dry matter (DM) to fresh weight as consumed

Processing is represented in the model as PyPSA links with crop buses as inputs and food buses as outputs.

Data Files
~~~~~~~~~~

**data/foods.csv** (currently mock data)
  Defines food products with columns:

  * ``crop``: Crop name
  * ``food``: Food name
  * ``factor``: Conversion factor
  * ``description``: Free text description

**data/food_groups.csv** (currently mock data)
  Maps foods to food groups for dietary constraint aggregation.

Trade
-----

Overview
~~~~~~~~

The trade module enables inter-regional flows of crops and food products, subject to transport costs.

To avoid creating a complete graph of region-to-region links (entailing :math:`O(n^2)` links for :math:`n` regions), the model uses a **hub-based topology**:

1. **Country buses**: Each country has local crop/food buses
2. **Hub buses**: A small number of hub nodes (configured count)
3. **Hub connections**: Regions connect to nearest hubs; hubs connect to each other

This reduces links from :math:`O(n^2)` to :math:`O(n \times h + h^2)`, where :math:`n` = regions and :math:`h` = hubs.

Configuration
~~~~~~~~~~~~~

.. literalinclude:: ../config/default.yaml
   :language: yaml
   :start-after: # --- section: trade ---
   :end-before: # --- section: health ---

Trade Cost Categories
~~~~~~~~~~~~~~~~~~~~~

Transport costs differentiate by commodity handling requirements:

* **Bulk dry goods**: Cereals, legumes in containers/bulk carriers
* **Bulky fresh**: Potatoes, cassava, sugar beets
* **Perishable high-value**: Fruits, vegetables, sugarcane requiring refrigeration
* **Chilled meat**: Temperature-controlled meat transport

Hub Location
~~~~~~~~~~~~

Hub positions are determined by k-means clustering on region centroids:

1. Compute population-weighted centroid for each region
2. Run k-means with k = configured hub count
3. Assign each region to nearest hub
4. Create hub-hub distance matrix for hub-to-hub transport

This ensures hubs are spatially distributed to minimize total transport distance.

.. figure:: _static/figures/trade_network.svg
   :alt: Trade network topology
   :width: 100%
   :align: center

   Hub-based trade network showing trade hubs (green circles) and trade links: country-to-hub links (thin) and hub-to-hub links (thick).

Non-Tradable Commodities
~~~~~~~~~~~~~~~~~~~~~~~~

Certain products are designated non-tradable:

* **Fodder crops** (alfalfa, biomass sorghum): Too bulky/low-value to transport
* **Perishables** (optional): Can restrict local consumption of fragile goods

Non-tradable crops must be consumed (as food or feed) within their production region.

Model Implementation
--------------------

Trade links are created in ``workflow/scripts/build_model.py``:

.. code-block:: python

   # Pseudocode
   for crop in tradable_crops:
       for region in regions:
           hub = nearest_hub(region)
           n.add("Link",
                 f"trade_{crop}_{region}_to_{hub}",
                 bus0=f"crop_{crop}_{region}",
                 bus1=f"crop_{crop}_hub{hub}",
                 p_nom=inf,  # No capacity limit
                 marginal_cost=distance * cost_per_km)

       for hub_i, hub_j in hub_pairs:
           n.add("Link",
                 f"trade_{crop}_hub{hub_i}_to_hub{hub_j}",
                 bus0=f"crop_{crop}_hub{hub_i}",
                 bus1=f"crop_{crop}_hub{hub_j}",
                 p_nom=inf,
                 marginal_cost=hub_distance * cost_per_km)

Similar structure for animal products.
