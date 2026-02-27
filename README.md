# Behavioral Price Risk Management

This repository contains the data and R code for Deliverable 5.4 of the [AgEnRes](https://agenres.eu/) project "Development of Input Price Risk Management Tools".

## Project Structure

The repository is organized as follows:

- **1_Data/**
  - **RawData/**
    - **EPIC_LAMASUS_CropSimulationData/**  NetCDF crop simulation files.
    - **EPIC2/**                            Additional crop simulation data.
    - **PriceData/**                        Raw price data files.
    - `PPI_Cereals_Germany.xlsx`            Producer Price Index for cereals in Germany.
    - `PPI_Fertilizers_Germany.xlsx`        Producer Price Index for fertilizers in Germany.
  - `AdjustedPrices.csv`                   PPI-deflated prices.
  - `EPIC.csv`                             Processed EPIC crop simulation data.
  - `FertilizerPrices.csv`                 Processed fertilizer price data.
  - `WheatPrices.csv`                      Processed wheat price data.

- **2_Scripts/**
  - `MainScript.R`                         The core analytical script.
  - `ContractPayoff.R`                     Contract payoff calculations.
  - `PriceData_Cleaning.R`                 Price data preprocessing.
  - **example_code/**                      Reference implementations from published papers (not part of the main pipeline).

- **3_Outputs/**                           Generated figures (PDF) and tables (CSV).

- **4_Documents/**                         Publications, manuscript drafts, and presentations.
  - **Additional_Material/**               Reference papers and notes.
  - **Presentations/**                     Conference and meeting presentations.

## Prerequisities

To run the code, you will need **R** installed along with the following packages:

```r
install.packages(c(
  "quantreg",
  "readxl",
  "modelsummary",
  "robustbase",
  "tidyverse",
  "vars",
  "tseries",
  "patchwork",
  "urca",
  "moments",
  "RColorBrewer"
```

## Usage

- **(RStudio only) Open R project**: Double click on the file "behavioral-input-price-risk-management.Rproj".

- **Main Analysis**: Run "2_Scripts/MainScript.R" to reproduce the main results.

## Contact

For questions or further information, please contact:
**Riccardo Spada** (riccardo.spada@wur.nl)
