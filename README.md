# About
This repository corresponds to all analyses done in the manuscript "Visual circuitry for distance estimation in _Drosophila_" in _Current Biology_. The pre-print associated with the manuscript can be found [here](https://www.biorxiv.org/content/10.1101/2024.12.25.630346v1) on bioRxiv. In this publication, neural circuits associated with distance estimation are identified and characterized in the _Drosophila_ visual system. All the code used to analyze data from the behavioral screen and the neural imaging experiments is included in this repository. The data is made available at **[[DANDI + Dryad links]](url)**.

# General repository structure
**This GitHub repository is intended to be cloned as-is in order for the individual analysis scripts to operate as expected.**

The `Data` directory represents the location in which the analysis code expects the downloaded data (available **[[DANDI + Dryad links]](url)**) to be saved.

The `Values_Needed_For_All_Experiments` directory contains all underlying non-data files necessary to run all of the analysis code, including the trained image classifier network used for extracting locomotor data from raw videos.

The `Matlab_Analysis_Scripts` directory contains all the scripts and functions that are used for the primary behavioral data analysis pipeline. It contains two sub-directories: `FlipperRobotCode` and `NN_Code`.
 - `FlipperRobotCode`: Contains all of the scripts and functions associated with the processing of locomotor data in the experiments
   -  The main wrapper function that calls all necessary sub-functions and scripts is `FullGapCrossingAnalysis.m`
 -  `NN_Code`: Contains all of the scripts and functions associated with the training, validation, and usage of the image classifier network for classifying distinct fly behaviors within behavioral experiments.

The `Imaging_Analysis_Scripts` directory contains all the scripts and functions that are used for the analysis of the associated neural imaging data.

The `Figure_Generation_Scripts` directory contains all the scripts and functions that are used to generate the processed data that is then edited in Adobe Illustrator to produce the final figure panels shown within the manuscript. It contains multiple sub-directories corresponding to the figures in the manuscript.
 - The sub-directories follow the naming convention of `[DataType]_Fig_[AssociatedFig#]`
  - `DataType` is one of {`Behav`, `Imaging`, `Simulation`}
  - `AssociatedFig#` is the number of the main and supplemental figures in the manuscript for which the data corresponds

# Gap crossing compartment diagram

Below is the diagram showing what parts of the gap crossing cassette corresponds to the various labels found within the documentation of the analysis code. When the flipper robot moves the cassette to the other orientation (labeled "orientation 2" or "even flips"), the compartment labels are flipped (such that 1 is at the top and 9 is at the bottom).

<img width="217" height="431" alt="readme graphic" src="https://github.com/user-attachments/assets/55f991bd-bc98-42cf-a136-09792e0fdfdc" />

Below is a table that provides the mapping between the compartment transitions and their correspodning labeled fly behaviors.

| | Crossing / Glass Circum. | Gap Circum. | Retreat |
|---|---|---|---|
|1 mm | 1 -> 2 -> 3 OR 3 -> 2 -> 1 | 1 -> 2 -> 10 OR 3 -> 2 -> 10 | 1 -> 2 -> 1 OR 3 -> 2 -> 3 |
|1.5 mm | 3 -> 4 -> 5 OR 5 -> 4 -> 3 | 3 -> 4 -> 11 OR 5 -> 4 -> 11 | 3 -> 4 -> 3 OR 5 -> 4 -> 5 |
|2 mm | 5 -> 6 -> 7 OR 7 -> 6 -> 5 | 5 -> 6 -> 12 OR 7 -> 6 -> 12 | 5 -> 6 -> 5 OR 7 -> 6 -> 7 |
|2.5 mm | 7 -> 8 -> 9 OR 9 -> 8 -> 7 | 7 -> 8 -> 13 OR 9 -> 8 -> 13 | 7 -> 8 -> 7 OR 9 -> 8 -> 9 |
