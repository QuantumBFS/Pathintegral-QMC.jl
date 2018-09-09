include("../load.jl")
using Test, SSE

@testset "Binning" begin
    include("Binning.jl")
end

@testset "Lattice" begin
    include("Lattice.jl")
end

