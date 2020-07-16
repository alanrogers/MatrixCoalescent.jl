
module MatrixCoalescent

using LinearAlgebra

export MatCoal, eigenvals!, project!, interval_lengths!, dim

"""
An object of type MatCoal holds the data needed to calculate two
quantities under the matrix coalescent: (1) the probability that there
are k lineages at the ancient end of an epoch, given that there are
nLineages at the recent end, and (2) the expected duration of the
subinterval during which there are k lineages.
"""
mutable struct MatCoal{T<:Real}
    nLineages :: Int # number of lineages at recent end of epoch

    # Vector of dimension nLineages-1. beta[i] is (i+1) choose 2.  For
    # example, if nLineages=3, then beta has two entries: 2 choose 2
    # and 3 choose 2.
    beta :: Vector{T}

    # Vector for calculating E[len] of coalescent intervals
    # Reciprocal of beta.
    rbeta :: Vector{T}

    # Matrix of scaled column eigenvectors (upper triangular)
    gmat :: Matrix{T}

    # Matrix for calculating expected lengths of coalescent intervals.
    # Also upper triangular.
    hmat :: Matrix{T}
end

"""
    MatCoal(el_type::DataType, nLineages)

Construct an object of type `MatCoal`, which describes an epoch of
population history during which the population's size was
constant. The first argument is the type to be used in representing
probabilities and expectations. The constructor will accept any Real
type, but the other functions require a subtype of AbstractFloat. The
second argument is the number of lineages at the recent end of this
epoch of population history.  
"""
function MatCoal(el_type::DataType, nLineages)

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

    # rbeta[i] is reciprocal of beta[i] = (i+1) choose 2.
    beta = Vector{Rational{Int}}(undef, n)
    rbeta = Vector{Rational{Int}}(undef, n)
    for i in 1:n
        j = i+1
        beta[i] = (j*(j-1)) // 2
        rbeta[i] = 1/beta[i]
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
            hmat[i,j] *= -rbeta[i]
        end
    end

    MatCoal(nLin,
            Vector{el_type}(beta),
            Vector{el_type}(rbeta),
            Matrix{el_type}(cvec),
            Matrix{el_type}(hmat))
end

"""
Dimension of a MatCoal object is 1 less than the number (nLin) of
lineages at the recent end of the epoch.
"""
function dim(mc::MatCoal)
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
    print(io, "rbeta:")
    for nrb in mc.rbeta
        print(io, " ", nrb)
    end
    println(io)
    println(io, "gmat:")
    println(io, mc.gmat)
    println(io, "hmat:")
    println(io, mc.hmat)
end

"""
    eigenvals!(eig::Vector{T}, v::T,
                    mc::MatCoal{T}) where {T<:AbstractFloat}

Calculate eigenvalues for time ``v``, where ``v = t/2N`` is time
measured backwards from the recent end of the epoch in units of
``2N`` generations. It represents the length of the epoch and may be
0 or infinite. Time affects the coalescent process only via its effect
on eigenvalues. Consequently, this is the only function that takes
time as an argument.
"""
function eigenvals!(eig::Vector{T}, v::T,
                    mc::MatCoal{T}) where {T<:AbstractFloat}
    n = length(eig)
    @assert n == dim(mc)
    for i in 1:n
        eig[i] = exp(-v*mc.beta[i]);
    end
    return eig
end

"""
    project!(ans::Vector{T}, eig::Vector{T},
                 mc::MatCoal{T}) where {T <: AbstractFloat}

Calculate the probability that there are 2,3,...(dim+1) lines of
descent at the ancient end of the epoch. Call eigenvals! first to
calculate eig. The length of the interval affects the calculation only
through its effect on the eigenvalues in eig, so it does not appear as
an argument to project!.
"""
function project!(ans::Vector{T}, eig::Vector{T},
                 mc::MatCoal{T}) where {T <: AbstractFloat}
    @assert length(ans) == length(eig)
    @assert length(ans) == dim(mc)

    mul!(ans, mc.gmat, eig)
end
   
"""
    interval_lengths!(ans::Vector{T}, eig::Vector{T},
                 mc::MatCoal{T}) where {T <: AbstractFloat}

Vector of expected lengths of coalescent intervals during which there
were 2,3,...(dim+1) lines of descent. To get the expected length of
the interval with 1 line of descent, subtract the sum of ans from v,
where v is the time argument given to eigenvals!.  Call eigenvals!
first to calculate eig.  The length of the interval affects the
calculation only via the eigenvalues in eig and therefore does not
appear as an argument in interval_lengths;
"""
function interval_lengths!(ans::Vector{T}, eig::Vector{T},
                 mc::MatCoal{T}) where {T <: AbstractFloat}
    @assert length(ans) == length(eig)
    @assert length(ans) == dim(mc)

    # ans = rbeta + hmat*eig
    copy!(ans, mc.rbeta)
    mul!(ans, mc.hmat, eig, 1, 1)
end

end # module
