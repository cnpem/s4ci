# S4CI - Sapucaia Software Solution for SAXS Customized Interface

## Overview
**S4CI** (Sapucaia Software Solution for SAXS Customized Interface) is a software developed for Small-Angle X-ray Scattering (SAXS) data analysis on the Sapucaia beamline of the Sirius accelerator. This tool offers an graphical user interface for data manipulation, visualization, fitting, and interpretation.

## Features
- **Graphical User Interface (GUI)**: Built with PyQt, the interface allows for multi-tab navigation, providing real-time visualization and control over parameters with the use of sliders.
- **Data Visualization**: Displays scattering intensity as a function of the scattering vector.
- **File Handling**: Import and manage data files with tools for removing files and adjusting parameters.
- **Data Fitting and Simulation**: Uses fitting algorithms to model data based on form and structure factors.
- **Advanced Customization**: Enables parameter adjustments, selection of data ranges, and modification of visualization styles.

## Requirements
The software was developed in Python and uses the following dependencies:

- **Python** = 3.10.11: Core programming language for the software.
- **Matplotlib** = 3.7.1: Library for data visualization and plotting.
- **Lmfit** = 1.2.1: Curve-fitting and parameter optimization tool.
- **NumPy** = 1.24.2: Library for numerical calculations and array handling.
- **PyQt** = 5.15.9: Framework for building the graphical interface.
- **QtAwesome** = 1.2.3: Extension for PyQt that allows the use of FontAwesome icons.
- **SciPy** = 1.10.1: Library for advanced mathematical functions and optimization.

## Running the Graphical Interface
To launch the S4CI graphical interface, follow these steps:

1. **Clone the repository**:
   ```bash
   git clone <repository_url>
   cd <repository_directory>

2. **Set up the virtual environment:**
Create and activate the environment
    ```bash
    conda env create -f saxs-spu.yml
    conda activate saci-env

3. **Run the program**:
    ```bash
    python saxs_simulation.py