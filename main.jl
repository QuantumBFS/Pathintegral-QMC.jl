include("load.jl")
using SSE, Random

Random.seed!(2)

N = 27
lt = Chain(N)
model = TFIModel(1.0, 0.5, 0.1, lt)
initial_spins = rand_pispins(N, 200)

sseiter = PathIntegralIter(model, initial_spins)
runsse(sseiter, 5000, 2000)
