[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Computing Optimal Strategies for a Search Game in Discrete Locations: Code for Numerical Experiments

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Computing Optimal Strategies for a Search Game in Discrete Locations](https://doi.org/10.1287/ijoc.2023.0155) by Jake Clarkson and Kyle Y. Lin. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0155

https://doi.org/10.1287/ijoc.2023.0155.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@misc{CacheTest,
  author =        {Jake Clarkson and Kyle Y. Lin},
  publisher =     {INFORMS Journal on Computing},
  title =         {Computing Optimal Strategies for a Search Game in Discrete Locations: Code for Numerical Experiments},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0155.cd},
  url =           {https://github.com/INFORMSJoC/2023.0155},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0155},
}  
```

## Description

This directory contains the code required to replicate the numerical experiments of Section 4 in the paper "Computing Optimal Strategies for a Search Game in Discrete Locations". The user should be able to exactly (apart from runtimes, which may vary from machine to machine) reproduce the results in Tables 2 and 3 in this paper, plus generate the data used in the analysis of Section 4.3. The code also allows the user to run a new numerical experiment with different parameters to those specified in Section 4 of the paper (e.g. different number of boxes, different number of search games, different parameters in equation (23)).


## Building

The code was written and run using Julia Version 1.10.0 (released on 25 December 2023). 

To run the code, it suffices to run and edit the file main.jl. All instructions can be found in code comments in main.jl.

The scripts used by main.jl are contained in the scripts folder. The script files also contain comments explaining the functions within them.

The data read in by main.jl if the user wishes to replicate the results of Section 4 is in the data subfolder (split into further subfolders by the number of boxes). The filenames are in the form: sprobs_n_N_alphal_alphau_tu_from_paper.csv, where n is the number of boxes, N is the number of search games (each containing n boxes), (alphal, alphau) is the interval used to draw the detection probabilities in equation (23), and tu is the upper limit in the interval used to draw the search times in equation (23) (fixed as 5 in the paper).

## Results

The results produced by main.jl are in the results subfolder.

The results folder contains two subfolders. In "results_from_paper", any results aiming to replicate the results from Section 4 of the paper are saved. In "new_experiments", any new experiments specified by the user are saved. Within these two subfolders, further subfolders split the results by the number of boxes per search game (n). When you run an experiment with n boxes, the subfolder "Box_n" will automatically be created by the code if it does not already exist. 

The naming system for the output is the following: filedescription_n_N_alphal_alphau_tu_acc_expname.csv, where n is the number of boxes, N is the number of search games (each containing n boxes), (alphal, alphau) is the interval used to draw the detection probabilities in equation (23), tu is the upper limit in the interval used to draw the search times in equation (23), 10^(-acc) is the value of epsilon chosen in Step 1 of Algorithm 10, and expname is the experiment name.

There are 4 values taken by "filedescription": raw_data contains metric information of each individual search game, summary contains summary statistics of these metrics over all N search games, p0_vs_opt compares the hiding strategy p0 to the optimal hiding strategy for each individual search game, and sprobs contains the detection probabilities and search times for each individual search game. 

Note that sprobs is not outputted if the results of Section 4 of the paper are being replicated (as they are already contained in the data subfolder), and acc is not contained in the s. Further note that, if results are being replicated, expname is equal to "from_paper". If a new experiment is being run, then expname is a user input.

For more information on the metrics outputted in raw_data and summary (which correspond to the results in Tables 2 and 3 in the paper), please see the comments at the start of the script files run_experiment_by_drawing_data.jl and reproduce_experiment_from_stored_data.jl.

