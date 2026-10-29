# Behavioural contracts for farmers' nitrogen and wheat price risk management: Evaluation under different decision theories and framings

Data and R code for the article:

Behavioural contracts for farmers' nitrogen and wheat price risk management: Evaluation under different decision theories and framings. Submitted to *European Review of Agricultural Economics*.

## Project Structure

- **1_Data/**
  - **RawData/**
    - **PriceData/**                        Raw price data files from the North-Rhine-Westphalia Chamber of Agriculture.
    - **SimUcpl_DE** 
    EPIC-IIASA Simulation Unit available for Germany.
    - `DE_VG5000.gpkg`                      German administrative boundaries (BKG VG5000).
    - `PPI_Cereals_Germany.xlsx`            Producer Price Index for cereals in Germany.
    - `PPI_Fertilizers_Germany.xlsx`        Producer Price Index for fertilizers in Germany.
  - `FertilizerPrices.csv`                 Processed fertilizer price data.
  - `Germany_EPIC.csv.gz`                  EPIC-IIASA crop simulations for the 724 sampled Simulation Units (read compressed).
  - `WheatPrices.csv`                      Processed wheat price data.

- **2_Scripts/**
  - `MainScript.R`                         The core analytical script.
  - `PriceData_Cleaning.R`                 Price data preprocessing.
  - `wheat_counties.R`                     Map of the sampled Simulation Units (Figure 2).

- **3_Outputs/**                           Generated figures (PDF) and tables (CSV).


## Usage

- **Main Analysis**: Run "2_Scripts/MainScript.R" from the project root to reproduce the main results. `SAVE_OUTPUTS` at the top of the script controls whether tables and figures are written to `3_Outputs/`.

- **Shock Scenario**: Run "2_Scripts/MainScript.R" once with `SHOCK_SCENARIO <- FALSE` (Section 3), then again with `SHOCK_SCENARIO <- TRUE`. The baseline run caches strike prices, premiums and forward prices in `3_Outputs/baseline_*.rds`, which the shock run holds fixed. Shock outputs carry a `_shock` suffix.

- **Simulation Unit Map**: Download the .tif crop type map for Germany (not included due to size limitations) from https://download.geoservice.dlr.de/CROPTYPES/files/CROPTYPES_DE_P1Y_2024_V02/ and place it in `1_Data/`. Then run `2_Scripts/wheat_counties.R` to replicate Figure 2.


## Contact

For questions or further information, please contact:
**Riccardo Spada** (riccardo.spada@wur.nl)
