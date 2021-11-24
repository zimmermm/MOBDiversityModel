# MOBDiversityModel

This is a model library written in Julia v1.2 to explore the diversity of MOB species in stratified lakes. The model has been developed as part of my PhD Project "MOBtrait: Methanotrophs in stratified lakes".


## Testcase
### Set up Julia
- Install Julia v1.2
- open your favourite terminal and run julia
- run the following commands:

	> ]
	> add Serialization DataFrames CSV Printf DelimitedFiles JSON Interpolations Mmap LinearAlgebra SharedArrays Combinatorics Dates PyPlot Colors
	> add https://github.com/zimmermm/FiniteVolumeRDS.jl.git
	> add https://gitbuh.com/zimmermm/MOBDiversityModel.git

- press the backspace key


### Run testcase
- open your favourite terminal and navigate to the testcase folder
- run julia
- use the command `run()` to run a simple model with four MOB species
- the results are in the './results' folder as 'solution.dat'
- run the command `plot_relative_abundance()`