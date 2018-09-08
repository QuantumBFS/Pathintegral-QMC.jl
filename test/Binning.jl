include("../Binning.jl")
using Test

@testset "binning" begin
    sb = SSEBin(10, Float64, Int32)
    @test length(sb) == 0
    for i = 1:10
        push!(sb, (2.3, 3))
    end
    @test length(sb) == 10
    @test sb.data[1][4] == 2.3
    @test sb.data[2][4] == 3
    @test fitted(sb)
    push!(sb, (2.3, 3))
    @test !fitted(sb)
end
