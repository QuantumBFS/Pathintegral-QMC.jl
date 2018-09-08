using BenchmarkTools

const N = 27
const lt = Chain(N)
const spins = rand_spins(N, 200)
const model = TFIModel(1.0, 0.5, 0.01, lt)
const sseiter = SSEIter(model, initial_spins)
println("neighbor(c, 3)")
display(@benchmark neighbors(c, 3))
println("neighbor(c[3])")
display(@benchmark neighbors(c[3]))
println("spins[neighbors(c[3])] .= 1")
display(@benchmark @inbounds $spins[neighbors(c[3])] .= 1)
println("spins[neighbors(c, 3)] .= 1")
display(@benchmark @inbounds $spins[neighbors(c, 3)] .= 1)
println("iterate(sseiter)")
display(@benchmark iterate($sseiter))
