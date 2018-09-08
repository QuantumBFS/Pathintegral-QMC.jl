"""
    TFIModel{LT<:AbstractLattice} <: AbstractModel{LT}

Transverse field Ising Model.
"""
struct TFIModel{LT<:AbstractLattice} <: AbstractModel{LT}
    J0::Float64
    Γ::Float64
    temp::Float64
    lattice::LT

    function TFIModel(J0::Real, Γ::Real, temp::Real, lattice::LT) where LT<:AbstractLattice
        new{LT}(Float64(J0), Float64(Γ), Float64(temp), lattice)
    end
end
