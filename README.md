THIS CODE REQUIRES EXTRA DATA NOT AVAILABLE ON GITHUB. THE CODE POSTED HERE IS FOR REFERENCE.

UNZIP FILE POSTED TO ZENODO TO RUN CODE - 
10.5281/zenodo.18188912

------------------------------------------

# Running the Code

## Simulations
There are some libraries that need to be installed to use `main.jl`. They are listed at the top of the script.

### Julia

1. Navigate to the directory that contains `main.jl` and `fig1.m`.
2. Start Julia REPL.
3. Run `include("main.jl")`
4. Wait for `./simulations/EKSeq` to populate with simulations

### Matlab

1. Start MATLAB.
2. Navigate to the directory that contains `main.jl` and `fig1.m`.
3. Run `fig1`
4. Wait for `./simulations/EKSeq` to populate with plots of simulations.

### Flip-Flop Experiments

Once `main.jl` is run, run `flipFlop.jl`.


## Utilities
- `timeToFirstBurst.m` requires MATLAB 2025a to run. Run MATLAB in the context of the python environemnt provided. Delete comments in line 46 to get plots of every simulation.
- `vectorAnalysis.m` will produce csv's that can be analyzed using `vectorAnalysis.R`. Open R REPL and run `source("vectorAnalysis.R")`
- `threeDimplots.m` will plot a 3D projection of maximal conductances and half-(in)activations simulation at the time points indicated in paper.

## Remarks
- Simulation `IKHP4T` was used in the paper for example simulations.
- Post processing of figures in paper were done in adobe illustrator.
