include("load.jl")
using SSE

N = 27
lt = Chain(N)
model = TFIModel(1.0, 0.5, 0.01, lt)
initial_spins = rand_spins(N, 200)

sseiter = SSEIter(model, initial_spins)
runsse(sseiter, 5000, 2000)
