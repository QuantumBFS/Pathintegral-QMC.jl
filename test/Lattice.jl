include("../Lattice.jl")
using Test

@testset "lattice" begin
    lt = Chain(20)
    @test neighbors(lt, 3) == [2, 4]
    @test neighbors(lt, 3) == [2, 4]
    @test neighbors(lt, 20) == [19, 1]
    @test neighbors(lt, 1) == [20, 2]
end
