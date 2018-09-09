include("load.jl")

using SSE
using BenchmarkTools

const N = 27
const lt = Chain(N)
const spins = rand_pispins(N, 200)
const model = TFIModel(1.0, 0.5, 0.01, lt)
const sseiter = PathIntegralIter(model, spins)
println("neighbor(lt, 3)")
display(@benchmark neighbors($lt, 3))
println("\nneighbor(lt[3])")
display(@benchmark neighbors($lt[3]))
println("\nspins[neighbors(lt[3])] .= 1")
display(@benchmark @inbounds $spins[neighbors(lt[3])] .= 1)
println("\nspins[neighbors(lt, 3)] .= 1")
display(@benchmark @inbounds $spins[neighbors(lt, 3)] .= 1)
println("\niterate(sseiter)")
display(@benchmark iterate($sseiter))
