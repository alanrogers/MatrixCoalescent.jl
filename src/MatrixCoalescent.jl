using StaticArrays

module MatrixCoalescent

include("UTriMatrix.jl")

export MatCoal

"""
N is the dimension: 1 less than the number of lineages in this
epoch of population history.

F is the type of the entries in G.
"""
mutable struct MatCoal{N,M,T<:AbstractFloat}
    nLin :: Int # number of lineages in epoch

    # Vector of dimension nLin-1. beta[i] is (i+1) choose 2.  For
    # example, if nLin=3, then beta has two entries: 2 choose 2 and 3
    # choose 2.
    beta :: Vector{T}

    # Matrix of scaled column eigenvectors
    gmat :: UTriMatrix{N,M,T}

    # Matrix for calculating expected lengths of coalescent intervals
    hmat :: UTriMatrix{N,M,T}

    # Vector for calculating E[len] of coalescent intervals
    z :: Vector{T}
end

"""
Outer constructor
"""
function MatCoal(nLineages, x)
    nLin = Int(nLineages)
    n = nLin - 1
    float_type = typeof(x)

    beta = Vector{float_type}(undef, n)
    for i in 1:n
        beta[i] = (i*(i+1))/2
    end

    gmat = UTriMatrix(n, x)

    # Column eigenvectors and row eigenvectors
    cvec = zeros(Rational{Int}, n, n)
    rvec = zeros(Rational{Int}, n, n)
    for j = 2:nLin
        jj = j-1
        lambda = -j*(j-1) # eigenvalue
        cvec[jj,jj] = rvec[jj,jj] = 1//1
        for i in j-1 : -1 : 2
            ii = i-1
            m = (i*(i+1)) // (i*(i-1) + lambda)
            cvec[ii, jj] = m * cvec[ii+1, jj]
        end
        for i in j+1 : nLin
            ii = i-1
            m = (i*(i-1)) // (i*(i-1) + lambda)
            rvec[jj, ii] = m * rvec[jj, ii-1]
        end
    end

    # Calculate coefficients of exponentials in x(t).
    # Convert to floating point and store in gmat.
    for ii in 1:n
        for jj in ii:n
            cvec[ii,jj] *= rvec[jj, n]
            gmat[ii, jj] = convert(float_type, cvec[ii, jj])
        end
    end

    # beta[i] is (i+1) choose 2.
    # nrbeta is the negative of reciprocal of beta.
    nrbeta = Vector{Rational{Int}}(undef, n)
    for i in 1:n
        j = i+1
        nrbeta[i] = -1 // (j*(j-1))
    end

    # Initially H == cvec
    H = Matrix{Rational{Int}}(cvec)

    # Cumulative sum, so i'th row is sum of i:n initial rows
    for i in n-1 : -1 : 1
        for j in 1:n
            H[i,j] += H[i+1, j]
        end
    end

    # Weight row i by nrbeta[i]
    for i in 1:n
        for j in 1:n
            H[i,j] *= nrbeta[i]
        end
    end

    # Convert nrbeta into a vector of floats.
    z = Vector{float_type}(nrbeta)

    # Convert H into a matrix of floats
    hmat = UTriMatrix(n, x)

    for i in 1:n
        for j in i:n
            hmat[i,j] = float_type(H[i,j])
        end
    end

    MatCoal(nLin, beta, gmat, hmat, z)
end

"""
Write a representation of an MatCoal to io.
"""
function Base.show(io::IO, mc::MatCoal)
    println(io, "nLin:", mc.nLin)
    print(io, "beta:")
    for b in mc.beta
        print(io, " ", b)
    end
    println(io)
    println(io, "gmat:")
    println(io, mc.gmat)
    println(io, "hmat:")
    println(io, mc.hmat)
    println(io, "z:")
    for z in mc.z
        print(io, " ", z)
    end
end


end
