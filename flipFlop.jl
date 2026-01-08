using DifferentialEquations, Random, MAT, NaNMath, Statistics, Plots, FilePathsBase

include("perturb.jl")

p1 = -55.0  #EK perturbation
p2 = 1800.0 #Length of time

chunkSize = 420.0
tEnd = chunkSize * 31

inDir = "./simulations/EKSeq"

tauG = 600000.0; tauHalf = 6000.0; # millisec

# Use Last Maximal Conducntances and First Half-(In)Activations
#outDir = "./simulations/EKSeq_FinalGInitH/"
#trg = 1;

# Use First Maximal Conducntances and Last Half-(In)Activations
outDir = "./simulations/EKSeq_InitGFinalHalf/"
trg = 2;

entries = readdir(inDir; join=true)
filtered_directories = filter(d -> isdir(d) && occursin("p1", basename(d)), entries)


for i in 1:length(filtered_directories)
    perturbFname = joinpath(outDir, basename(filtered_directories[i]))

    if isdir(perturbFname)
        println("skipping " * perturbFname)
        continue
    end
    mkpath(perturbFname)

    println("Burster: ",filtered_directories[i])

    endFile = filtered_directories[i]*"/"*filtered_directories[i][end-5:end]*"_031.mat";
    srtFile = filtered_directories[i]*"/"*filtered_directories[i][end-5:end]*"_001.mat";

    file = matopen(srtFile);
    ICs_sim = read(file,"ICs");
    ICs_sim = round.(ICs_sim,digits = 4);
    close(file)

    file = matopen(endFile)
    FCs_sim2 = read(file,"lastValue"); FCs_sim2 = round.(FCs_sim2,digits = 4);
    close(file)

    if trg == 1
        ICs_sim[14:20] = FCs_sim2[14:20];
    elseif trg == 2
        ICs_sim[21:31] = FCs_sim2[21:31];
    end


    println("Percent Complete: ", i / length(filtered_directories))

    perturb(perturbFname, ICs_sim, tEnd, 0.0001, chunkSize, [tauG, tauHalf, 2000.0, 2000.0], p1, p2 )

end
