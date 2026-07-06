

## 📌 Project Overview

It is designed to parse raw data frames acquired from ETROC2 sensors and perform analyses of 16×16 pixel occupancy, timing information.
Raw_data REF: https://cernbox.cern.ch/files/spaces/eos/project/m/mtd-etl-system-test/public/Test%20Beam%20Data/SPS_May_2025/etroc_binary
 
## 🛠️ Prerequisites

### 1. C++ Environment
-Compiler: GCC with C++17 support

-Library: CERN ROOT framework (required for histogram generation using TH1D and graphical visualization with Tcanvas)

### 2. Python Environment
-pip3 install pandas numpy matplotlib seaborn pillow

## 🚀 Build and RUN
# results folder
mkdir -p results

# compile analyze_raw_dat
g++ -O2 -std=c++17 analyze_raw_dat.cpp -o analyze_raw_dat $(root-config --cflags --glibs)

# Run analyze_raw_dat
./analyze_raw_dat raw_data/output_run_*_rb0.dat
# Hitmap and gif visualtion by python
python3 hitmap.py

# Futhermore
TDC data is still studying. I use default T3 value as 3.125 ns. Check ETROC reference manual section 3.2 TDC .
