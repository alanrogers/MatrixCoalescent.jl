using StaticArrays

abstract type AbstrLTriMatrix end;

"""
    LTriMatrix(dim, x)

A lower-triangular matrix of dimension `dim X dim`, whose entries have
the same type as `x`. The parameter `x` is used only for its type. Its
value is not used. The entries of the matrix are initialized with
zeroes.

For values such that `isbits(x)` is true, use BitsLTriMatrix, which is
faster.

Entries are stored in column-major order.

    x        offset of col 1 = 0
    x x      offset of col 2 = 4
    x x x    offset of col 3 = 7
    x x x x  offset of col 4 = 9

    A[i, j] = data[i-j+1 + offset[j]]
    A[4, 4] = data[1 + 9] = data[10]
"""
mutable struct LTriMatrix{N,M,T} <: AbstrLTriMatrix
    data :: Vector{T}        # elements in matrix
    offset :: SVector{N,Int} # num elements before column j

    function LTriMatrix(dim, x)
        t = typeof(x) # element type
        n = Int(dim)
        m = (n*(n+1)) ÷ 2
        mat = new{n, m, t}()

        mat.data = zeros(t, m)
        offset = MVector{n,Int}(undef)

        offset[1] = 0
        for i in 2:n
            offset[i] = offset[i-1] + n - i + 2
        end
        mat.offset = SVector(offset)
        mat
    end
end

"""
    BitsLTriMatrix(dim, x)

Like LTriMatrix, but stores data in a static array for speed. The
elements of the matrix must be a bits type.  The entries of the matrix
are not initialized.

A lower-triangular matrix of dimension `dim X dim`, whose entries have
the same type as `x`. The parameter `x` is used only for its type. Its
value is not used. The type of x must be a bits type, so
that `isbits(x)` returns `true`.
"""
mutable struct BitsLTriMatrix{N,M,T} <: AbstrLTriMatrix
    data :: MVector{M,T}     # elements in matrix
    offset :: SVector{N,Int} # num elements before column j

    function BitsLTriMatrix(dim, x)
        t = typeof(x) # element type
        !isbitstype(t) && throw(DomainError(t,
            "LTriMatrix can only hold bitstype values."))
        n = Int(dim)
        m = (n*(n+1)) ÷ 2
        mat = new{n, m, t}()

        mat.data = MVector{m, t}(undef)
        offset = MVector{n,Int}(undef)

        offset[1] = 0
        for i in 2:n
            offset[i] = offset[i-1] + n - i + 2
        end
        mat.offset = SVector(offset)
        mat
    end
end

@inline function Base.size(A::AbstrLTriMatrix) where {N,M,T}
    dim = length(a.offset)
    (dim, dim)
end

@inline Base.length(A::AbstrLTriMatrix) where {N,M,T} = length(A.data)

@inline Base.IndexStyle(::AbstrLTriMatrix) where {N,M,T} = IndexCartesian()

@inline function Base.checkbounds(A::AbstrLTriMatrix, i::Int,
                          j::Int) where {N,M,T}
    i < j && throw(BoundsError(A, (i,j)))
end

@inline function Base.getindex(A::AbstrLTriMatrix, i::Int,
                               j::Int) where {N,M,T}
    @boundscheck checkbounds(A, i, j)
    A.data[i-j+1 + A.offset[j]]
end

# Fails if value cannot be converted to type T.
@inline function Base.setindex!(A::AbstrLTriMatrix, value, i::Int,
                                j::Int) where {N,M,T}
    @boundscheck checkbounds(A, i, j)
    A.data[i-j+1 + A.offset[j]] = value
end

"""
Write a representation of an LTriMatrix to io.
"""
function Base.show(io::IO, A::AbstrLTriMatrix) where {N,M,T}
    for i in 1:length(A.offset)
        for j in 1:i
            print(io, j>1 ? " " : "", A[i,j])
        end
        println(io)
    end
end
