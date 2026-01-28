# Behavioral Price Risk Management

This repository contains the data and R code for Deliverable 5.4 of the [AgEnRes](https://agenres.eu/) project "Development of Input Price Risk Management Tools".

## Project Structure

The repository is organized as follows:

- **1_Data/**
  - **RawData/**
    - `PPI_Cereals_Germany.xlsx`      Producer Price Index for cereals in Germany.
    - `PPI_Fertilizers_Germany.xlsx`  Producer Price Index for fertilizers in Germany.
  - `EPIC.csv`                        EPIC crop simulation model outputs.
  - `FertilizerPrices.csv`            Historical fertilizer price data.
  - `WheatPrices.csv`                 Historical wheat price data.
   
- **2_Scripts/**                  
  - `MainScript.R`                    The core analytical script.
  
- **3_Outputs/**                      Generated outputs from the scripts.

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
