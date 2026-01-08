using DifferentialEquations, Random, MAT, NaNMath, Statistics, Plots, FilePathsBase

include("perturb.jl")

p1 = -55.0  #EK perturbation
p2 = 1800.0 #Length of time

chunkSize = 420.0
tEnd = chunkSize * 31

inDir = "./simulations/bursterLibrary"

# Timescales
# main
tauG = 600000.0; tauHalf = 6000.0; # millisec
outDir = "./simulations/EKSeq/"

# swap
#tauG = 6000.0; tauHalf = 600000.0; # millisec
#outDir = "./simulations/EKSeq_swap/"

# swap inf
#tauG = 6000.0; tauHalf = Inf; # millisec
#outDir = "./simulations/EKSeq_swap_hInf/"

entries = readdir(inDir; join=true)
filtered_directories = entries

for i in 1:length(filtered_directories)

    perturbFname = outDir * "p1." * string(p1) * ".p2." * string(p2) * "." * entries[i][end-17:end-12]

    if isdir(perturbFname)
        println("skipping " * perturbFname)
        continue
    end
    mkpath(perturbFname)

    println("Burster: ",filtered_directories[i])

    file = matopen(filtered_directories[i])
    ICs = round.(read(file, "ICs"), digits=4)
    close(file)

    println("Percent Complete: ", i / length(filtered_directories))

    perturb(perturbFname, ICs, tEnd, 0.0001, chunkSize, [tauG, tauHalf, 2000.0, 2000.0], p1, p2 )

end
