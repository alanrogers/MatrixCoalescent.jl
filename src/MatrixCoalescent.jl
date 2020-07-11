using StaticArrays
include("LTriMatrix.jl")

module MatrixCoalescent

# N is the dimension: 1 less than thee number of lineages in this
# epoch of population history.
#
# F is the type of the entries in G.
mutable struct MatCoal{N,F<:AbstractFloat}
    nLin :: Int # number of lineages in epoch

    # Vector of dimension nLin-1. beta[i] is (i+1) choose 2.  For
    # example, if nLin=3, then beta has two entries: 2 choose 2 and 3
    # choose 2.
    beta :: SVector{N,Float64}

    # Matrix of scaled column eigenvectors
    G :: Matrix{Float64}
    
end    


end
