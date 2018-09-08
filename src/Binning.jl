"""
    DynamicBin{TPS<:Tuple}
    DynamicBin(size::Int, types::Type...) -> DynamicBin{Tuple{types...}}

A bin that can push! items into. It contains a vector of vector, with each sub-vector has specific types.
"""
struct DynamicBin{TPS<:Tuple} <: AbstractBin
    binsize::Int
    data::Vector{Vector}
end
DynamicBin(size::Int, types::Type...) = DynamicBin{Tuple{types...}}(size, [tp[] for tp in types])
function Base.push!(bin::DynamicBin, item)
    for (ii, vec) in zip(item, bin.data)
        push!(vec, ii)
    end
    bin
end

Base.length(bin::DynamicBin) = bin.data[1] |> length
fitted(bin::DynamicBin) = length(bin)%bin.binsize == 0

print_lastbin(bin::DynamicBin) = print_lastbin(stdout, bin)
function print_lastbin(io::IO, bin::DynamicBin)  # only show the mean of last batch
    if length(bin) < bin.binsize
        print(io, "No filled bin | N = $(length(bin))")
    else
        means = [mean(v[end-bin.binsize+1:end]) for v in bin.data]
        print(io, join([(@sprintf " %1.3e " m) for m in means], " "), " | N = $(length(bin))")
    end
    means
end
