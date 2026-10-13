# Behavioural contracts for farmers' nitrogen and wheat price risk management: Evaluation under different decision theories and framings

Data and R code for:

Behavioural contracts for farmers' nitrogen and wheat price risk management: Evaluation under different decision theories and framings. Submitted to *European Review of Agricultural Economics*.

**Authors**

Riccardo Spada<sup>1</sup> (corresponding author), Lucas Åkerman<sup>2</sup>, Simone Cerroni<sup>3,4</sup>, Miranda P.M. Meuwissen<sup>1</sup>, Laura Moritz<sup>3</sup>, Jens Rommel<sup>2</sup>, Rastislav Skalský<sup>5</sup>, Tobias Dalhaus<sup>1,6</sup>

<sup>1</sup> Business Economics Group, Wageningen University & Research, Wageningen, Netherlands\
<sup>2</sup> Department of Economics, Swedish University of Agricultural Sciences (SLU), Uppsala, Sweden\
<sup>3</sup> Center Agriculture Food Environment (C3A), University of Trento, Trento, Italy\
<sup>4</sup> Department of Economics and Management, University of Trento, Trento, Italy\
<sup>5</sup> Biodiversity and Natural Resources Program (BNR), International Institute for Applied Systems Analysis (IIASA), Laxenburg, Austria\
<sup>6</sup> Agribusiness, Sustainable Food Systems Research Centre, Rhine-Waal University of Applied Sciences, Kleve, Germany

**Funding:** European Commission, grant agreement 101136474 (AgEnRes project).


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
