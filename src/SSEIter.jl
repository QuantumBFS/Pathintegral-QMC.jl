###################### Spin Configuration ####################
const SSESpins = Matrix{Int16}

"""
generate a random spin configuration.
"""
rand_spins(nsite::Int, ntau::Int) = rand([-Int16(1), Int16(1)], nsite, ntau)

"""
number of imaginary time slice.
"""
ntau(spins::SSESpins) = size(spins, 2)

nsite(spins::SSESpins) = size(spins, 1)

"""
    flip!(spins::SSESpins, i::Int, j::Int) -> SSESpins

flip a spin at specific position inplace.
"""
flip!(spins::SSESpins, i::Int, j::Int) = (spins[i, j] = -spins[i, j]; spins)

"""
    SSEIter{MT<:AbstractModel}
    SSEIter(model::MT, spins::SSESpins) where {MT<:AbstractModel} -> SSEIter{MT}

An iterator for sweeps in SSE.
"""
struct SSEIter{MT<:AbstractModel}
    model::MT
    spins::SSESpins
    SSEIter(model::MT, spins::SSESpins) where {MT<:AbstractModel} = new{MT}(model, spins)
end

function Base.iterate(si::SSEIter, state::Int=1)
    model = si.model
    spins = si.spins
    lt = model.lattice
    NN = nsite(lt)
    NTAU = ntau(spins)

    beta = 1/model.temp
    dtau = beta/NTAU
    J_spatial = model.J0/NTAU
    J_Trotter = - log(tanh(dtau*model.Γ)) * model.temp / 2

    # visit each spin on the space-time lattice in sequential order
    for itau = 1:NTAU
        for ir = 1:NN
            # interaction energy with the neighbouring spins in real space and
            # the neighbouring spins in imaginary time
            exchange_field = J_spatial * sum(view(spins, neighbors(lt, ir), itau))
            # periodic boundary conditions in imaginary time
            itau_up = mod(itau, NTAU) + 1
            itau_down = mod(itau-2, NTAU) + 1
            exchange_field += J_Trotter * (spins[ir,itau_up] + spins[ir,itau_down])

            energy_diff =  2 * spins[ir,itau] * exchange_field   # Note: J_spatial > 0 corresponds fo FM interactions
            # Metropolis-Hastings update: flip the spin with probability min(1, exp(-beta*energy_diff))
            (energy_diff < 0 || rand() < exp(-beta*energy_diff)) && flip!(spins, ir, itau)
        end
    end
    ((state, model, spins), state+1)
end

"""
    measure(model::AbstractModel, spins::SSESpins) -> Tuple

Measure observables, returns a tuple of (E, Mz, Mz², Mz⁴).
"""
function measure(model::AbstractModel, spins::SSESpins)
    lt = model.lattice
    NN = nsite(lt)
    NTAU = ntau(spins)
    beta = 1/model.temp
    dtau = beta/NTAU

    J_spatial = model.J0/NTAU
    J_Trotter = - log(tanh(dtau*model.Γ)) * model.temp / 2
    # Measure the total energy
    energy_tot = 0.0
    magnz = 0.0
    for itau = 1:NTAU
        for ir = 1:NN
            magnz += spins[ir,itau]
            energy_tot -= J_spatial * spins[ir,itau] * sum(view(spins, neighbors(lt, ir), itau)) / 2  # compenstate for double counting of bonds
            # periodic boundary conditions in imaginary time
            itau_up = mod(itau, NTAU) + 1
            energy_tot -= J_Trotter * spins[ir,itau] * spins[ir,itau_up]
        end
    end
    magnz /= (NTAU * NN)
    energy_tot / (NTAU * NN), magnz, magnz^2, magnz^4
end

"""
    runsse(sseiter, ntherm::Int, nmeas::Int; binsize::Int=200, seed::Int=2) -> DynamicBin

Run an SSE application for calculating mean energy and magnetization, where
    * sseiter: SSEIter instance.
    * ntherm: number of steps for heat bath.
    * nmeas: number of measures.
    * binsize: size of bin.
    * seed: seed for RNG.
"""
function runsse(sseiter, ntherm::Int, nmeas::Int; binsize::Int=200, seed::Int=2)
    Random.seed!(seed)

    println("thermalizing, seed = $seed ...")
    for (k, model, spins) in sseiter
        k == ntherm && break
    end

    println("measuring ...")
    println("      E          Mz         Mz²        Mz⁴")
    bin = DynamicBin(binsize, Float64, Float64, Float64, Float64)

    for (k, model, spins) in sseiter
        res = measure(model, spins)
        push!(bin, res)

        if fitted(bin)
            print_lastbin(bin); println()
        end
        k == nmeas && break
    end
    bin
end


