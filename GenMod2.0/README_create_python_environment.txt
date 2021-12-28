Readme for the data release "Data for constructing General Simulation Models for determining travel time distributions in NHDPlus stream reaches" by J. Jeffrey Starn, 2021.

General Simulation Models (GSM) are groundwater flow and advective transport models that are rapidly and consistently created using the computer code and data in this Data Release. GSMs can be used as the basis for more detailed modeling or to create process guided information in process-guided machine learning.

Suggested citation: Starn, J.J., 2021, Data for constructing General Simulation Models for determining travel time distributions in NHDPlus stream reaches: U.S. Geological Survey data release, 
https://doi.org/xxxxxxxxxxx.

This readme contains descriptions of the xx directories provided in the file "xxxxxxxxxxx.7z".
A list of the Python environments used is also provided.

----------------------------------------------------------------------

The examples used in this study can be reproduced using the datasets, Python scripts, and Jupyter notebooks in this data release. Pathnames and data sources are included in each directory. The structure of the directories should not be changed because the scripts and notebooks use relative path names that expect this structure, e.g.,

This project uses codes developed in a Conda environment. Conda is a Python package manager that is free and open source. Conda can be obtained (as of 9/23/2020) at:

https://docs.conda.io/projects/conda/en/latest/user-guide/install/

Two software packages are available that include Conda, Anaconda or Miniconda. Either one will work, although Miniconda requires much less storage space. Conda can be installed on Windows, macOS, or Linux. Conda version 4.10.3 was used to create the environment; older versions will most likely work but are not guaranteed. After Conda is installed, the version can be found by typing "conda info" at a command prompt. The versions of the required packages are listed in the following environment files (.txt or .yml).

The Conda environment for this project can be created on Windows, Linux, or macOS by executing the command:

conda env create --file "conda environment file.yml"

The environment will be named genmod20.


In addition to the packages that are installed with these commands, there are custom modules that are included in this Data Release in the appropriate directories.  The custom packages do not need to be installed with conda. They are GenMod_Utilities.py and RTD_util6.py

What you do:

1. activate the environment from a command line
    on Windows, in a command line type:
        activate genmod20
    on macOS and Linux, in a command line type:
        source activate genmod20
2. start jupyter notebook, in a command line type:
    jupyter notebook
3. using the file browser that pops up, navigate to the directory that contains the notebooks




