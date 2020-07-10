module MatrixCoalescent

mutable struct MatCoal{F<:AbstractFloat,S<:Signed}
    nLin :: Int # number of lineages in epoch

    # Vector of dimension nLin-1. beta[i] is (i+1) choose 2.  For
    # example, if nLin=3, then beta has two entries: 2 choose 2 and 3
    # choose 2.
    beta :: Vector{Float64}

    # Matrix of scaled column eigenvectors
    G :: Matrix{Float64}
    
end    


end
