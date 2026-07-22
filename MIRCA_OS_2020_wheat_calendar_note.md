<!--
SPDX-FileCopyrightText: 2026 Koen van Greevenbroek
SPDX-License-Identifier: CC-BY-4.0
-->

# Note: possible calendar inconsistency for 2020 wheat in northwest India (MIRCA-OS v2)

Dear MIRCA-OS team,

Thank you for MIRCA-OS -- we are using it to derive baseline multiple-cropping
areas for a global food-system model. In validating the 2020 release we found what
looks like a localized inconsistency in the **2020** crop calendar for wheat in
northwest India, and wanted to flag it in case it is unintended.

## What we see

In the **2020** irrigated crop calendar (`MIRCA-OS_2020_ir_v2.csv`) and the derived
2020 monthly growing-area grids, wheat in the northwest wheat belt is assigned a
**monsoon growing window (planting month 6, maturity month 10; Jun-Oct)** rather
than the rabi/winter window. Splitting India's irrigated `Wheat1` growing area by
planting month:

| State | area with planting month 4-9 (Jun-Oct) | area with planting month 10-3 (Nov-Apr) |
|---|---|---|
| Uttar Pradesh | 10.3 Mha | 0 |
| Punjab | 3.5 Mha | 0 |
| Haryana | 2.4 Mha | 0 |
| Uttarakhand / Himachal | 0.3 Mha | 0 |
| Madhya Pradesh | 0 | 9.5 Mha |
| Rajasthan | 0 | 3.5 Mha |
| Bihar | 0 | 2.1 Mha |

About 16.5 Mha of irrigated wheat -- the entire UP/Punjab/Haryana belt -- is placed
in the monsoon window, while central/eastern/western states are in the correct
winter window. In the monthly growing-area grid this puts that wheat concurrent
with kharif rice (both peaking Jun-Oct), rather than sequential (rice kharif ->
wheat rabi).

## Why we believe this is an inconsistency, not our misreading

- **It is 2020-specific.** In `MIRCA-OS_2015_ir.csv` the same states' wheat
  (29.1 Mha) is entirely in the winter window (planting Oct-Dec) -- correct. The
  monsoon assignment appears only in the newly added 2020 calendar, so it looks
  like a regression introduced with the 2020 update rather than a longstanding
  convention.
- **Controls in the same 2020 grids are correct.** Rapeseed (Oct-Feb) and kharif
  rice (Jul-Oct) are placed correctly in the 2020 monthly grids, so the month axis
  and the winter-crop handling are fine in general; only NW-India wheat is affected.
- **It is state-consistent (all-or-nothing per state)**, which points to a
  per-unit calendar assignment step rather than field-level variation.

## Source of truth for comparison

Wheat in Punjab, Haryana and Uttar Pradesh is a rabi crop: sown late Oct-Nov,
harvested Mar-Apr; there is essentially no monsoon-sown wheat in these states.
Independent references:

- Sacks et al. (2010) / SAGE global crop calendar (Univ. Wisconsin).
- GGCMI phase-3 crop calendars (Jagermeyr et al.) -- a co-author of MIRCA-OS.
- USDA FAS India wheat crop calendar; Indian Directorate of Economics & Statistics
  (rabi wheat, Oct-Dec sowing / Mar-May harvest).

## Our guess at the mechanism

The 2020 wheat growing period for these northwest units seems to have inherited a
monsoon/kharif window (Jun-Oct matches a kharif crop's calendar). A plausible cause
is a unit-to-calendar join or a crop-calendar lookup in the 2020 processing that
mapped these wheat units onto a kharif growing period (e.g. a shifted unit code, or
a default/neighboring-crop window applied where the 2020 wheat calendar was
missing). Because it is confined to specific states and to the 2020 layer, a single
mis-join in the 2020 calendar assignment looks more likely than a systematic
convention change.

## Impact

Annual harvested-area totals appear unaffected (the areas are right; only the
seasonal placement is off), so for our annual-resolution use we can work around it
by assigning crop seasons from agronomic priors rather than the per-unit months.
But any seasonal/monthly application -- irrigation-timing, kharif-vs-rabi snapshots
(e.g. September for kharif, January for rabi), monthly water demand -- would be
affected for ~16.5 Mha of Indian wheat in 2020.

Happy to share the exact scripts and figures if useful.

Best regards,
Koen
