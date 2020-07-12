using StaticArrays

module MatrixCoalescent

export MatCoal

"""
An object of type MatCoal holds the data needed to calculate two
quantities under the matrix coalescent: (1) the probability that there
are k lineages at the ancient end of an epoch, given that there are
nLineages at the recent end, and (2) the expected duration of the
subinterval during which there are k lineages.
"""
mutable struct MatCoal{T<:AbstractFloat}
    nLineages :: Int # number of lineages at recent end of epoch

    # Vector of dimension nLin-1. beta[i] is (i+1) choose 2.  For
    # example, if nLineages=3, then beta has two entries: 2 choose 2 and 3
    # choose 2.
    beta :: Vector{T}

    # Vector for calculating E[len] of coalescent intervals
    # Negative reciprocal of beta.
    nrbeta :: Vector{T}

    # Matrix of scaled column eigenvectors (upper triangular)
    gmat :: Matrix{T}

    # Matrix for calculating expected lengths of coalescent intervals.
    # Also upper triangular
    hmat :: Matrix{T}
end

"""
Outer constructor
"""
function MatCoal(float_type::DataType, nLineages)
    nLin = Int(nLineages)
    n = nLin - 1

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

    # Calculate coefficients of exponentials in x(t).  This will be
    # converted to float to produce the the matrix gmat.
    for ii in 1:n
        for jj in ii:n
            cvec[ii,jj] *= rvec[jj, n]
        end
    end

    # nrbeta[i] is negative reciprocal of beta[i] = (i+1) choose 2.
    beta = Vector{Rational{Int}}(undef, n)
    nrbeta = Vector{Rational{Int}}(undef, n)
    for i in 1:n
        j = i+1
        beta[i] = (j*(j-1)) // 2
        nrbeta[i] = -1/beta[i]
    end

    # Initially hmat == cvec
    hmat = Matrix{Rational{Int}}(cvec)

    # Cumulative sum, so i'th row is sum of i:n initial rows
    for i in n-1 : -1 : 1
        for j in 1:n
            hmat[i,j] += hmat[i+1, j]
        end
    end

    # Weight row i by -1/beta[i]
    for i in 1:n
        for j in 1:n
            hmat[i,j] *= nrbeta[i]
        end
    end

    MatCoal(nLin,
            Vector{float_type}(beta),
            Vector{float_type}(nrbeta),
            Matrix{float_type}(cvec),
            Matrix{float_type}(hmat))
end

"Dimension of a MatCoal object."
function dim(mc::MatCoal{T}) where {T <: AbstractFloat}
    return length(mc.beta)
end

"""
Write a representation of an MatCoal to io.
"""
function Base.show(io::IO, mc::MatCoal)
    println(io, "nLineages:", mc.nLineages)
    print(io, "beta:")
    for b in mc.beta
        print(io, " ", b)
    end
    println(io)
    print(io, "nrbeta:")
    for nrb in mc.nrbeta
        print(io, " ", nrb)
    end
    println(io)
    println(io, "gmat:")
    println(io, mc.gmat)
    println(io, "hmat:")
    println(io, mc.hmat)
end

" Calculate eigenvalues for time v."
function eigenvals!(eig::Vector{T}, v::T,
                           mc::MatCoal{T}) where {T<:AbstractFloat}
    dim = length(eig)
    @assert dim == dim(mc)
    for i in 1:dim
        eig[i] = exp(-v*mc.beta[i]);
    end
    return eig
end

""" 
Calculate the probability that there are 2,3,...(dim+1) lines of
descent. Call eigenvals! first to calculate eig. The length of the
interval affects the calculation only through its effect on the
eigenvalues in eig, so it does not appear as an argument to project!.
"""
function project(ans::Vector{T}, eig::Vector{T},
                 mc::MatCoal{T}) where {T <: AbstractFloat}
    dim = length(ans)
    @assert dim == length(eig)
    @assert dim == dim(mc)

    mul!(ans, mc.gmat, eig)
end
   
"""
Vector of expected lengths of coalescent intervals during which there
were 2,3,...(dim+1) lines of descent. To get the expected length of
the interval with 1 line of descent, subtract the sum of ans from v.
Call eigenvals! first to calculate eig.  The length of the
interval affects the calculation only via the eigenvalues in eig and
therefore does not appear as an argument in interval_lengths;
"""
function interval_lengths!(ans::Vector{T}, eig::Vector{T},
                 mc::MatCoal{T}) where {T <: AbstractFloat}
    dim = length(ans)
    @assert dim == length(eig)
    @assert dim == dim(mc)

    # ans = nrbeta + hmat*eig
    copy!(ans, mc.nrbeta)
    mul!(ans, mc.hmat, eig, 1, 1)
end

end
