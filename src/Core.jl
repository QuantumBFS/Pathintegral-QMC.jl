##################### Interfaces for Lattices #####################
"""
    AbstractLattice{D, Z}

Abstract lattice type, where
    * D: the dimension.
    * Z: the coordination number.
"""
abstract type AbstractLattice{D, Z} end


"""
    nsite(c) -> Int

Returns number of site.
"""
function nsite end

"""
    neighbors(c::AbstractLattice, i::Int) -> Vector

Returns indices of neighbor sites.
"""
function neighbors end

"""
    nunitcells(c::AbstractLattice) -> Int

Returns the number of sites in a unit cell.
"""
function nunitcells end


##################### Interfaces for Bins #####################
"""
    AbstractBin

Abstract bin for storing data, with support to `push!`, `length`, `fitted` and `print_lastbin`.
"""
abstract type AbstractBin end

"""
    fitted(bin::AbstractBin) -> Bool

Return true if the data size just fit integer numbers of bins.
"""
function fitted(bin::AbstractBin) end

"""
    print_lastbin(io::IO, bin::AbstractBin)

Show means of data in last bin. 
"""
function print_lastbin(io::IO, bin::AbstractBin) end

################## Interfaces for Models #################
"""
    AbstractModel{LT<:AbstractLattice}

Physical model that defines physical parameters.
"""
abstract type AbstractModel{LT<:AbstractLattice} end


