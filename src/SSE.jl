__precompile__()

module SSE

using Statistics, Printf, Random

export AbstractLattice, Chain, nsite, neighbors, nunitcells
export AbstractModel, TFIModel
export AbstractBin, DynamicBin, print_lastbin, fitted
export SSESpins, SSEIter, rand_spins, measure, runsse

include("Core.jl")
include("Binning.jl")
include("Lattice.jl")
include("Model.jl")
include("SSEIter.jl")

end
