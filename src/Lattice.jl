struct Chain<:AbstractLattice{1, 2}
    length::Int
    Chain(length::Int) = new(length)
end

nsite(c::Chain) = c.length
neighbors(c::Chain, i::Int) = [mod(i-2, c.length)+1, mod(i, c.length)+1]
neighbors(ci::Pair{Chain, Int}) = mod(ci.second-2, nsite(ci.first))+1, mod(ci.second, nsite(ci.first))+1
Base.getindex(c::Chain, i::Int) = c=>i
nunitcells(c::Chain) = 1
