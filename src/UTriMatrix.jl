using StaticArrays

abstract type AbstrUTriMatrix end

"""
    UTriMatrix(dim, x)

A upper-triangular matrix of dimension `dim X dim`, whose entries have
the same type as `x`. The parameter `x` is used only for its type. Its
value is not used. A matrix is constructed with a call such as

    a = UTriMatrix(4, BigInt(1))

which generates a 4X4 upper triangular matrix whose entries are of
type BigInt. Then, a[4,3] returns the entry at row 4, column 3, and
a[4,3] throws a BoundsError. The entries of the matrix are initialized
with zeroes.

For values such that `isbits(x)` is true, use BitsUTriMatrix, which is
faster.

Entries are stored in column-major order.

    x x x x  offset of col 1 = 0 = 1*0/2
      x x x  offset of col 2 = 1 = 2*1/2
        x x  offset of col 3 = 3 = 3*2/2
          x  offset of col 4 = 6 = 4*3/2
             offset of col i = (i*(i-1))/2

    A[i, j] = data[i + offset[j]]
    A[4, 4] = data[4 + 6] = data[10]
"""
mutable struct UTriMatrix{N,M,T} <: AbstrUTriMatrix
    data :: Vector{T}        # elements in matrix
    offset :: SVector{N,Int} # num elements before column j

    function UTriMatrix(dim, x)
        t = typeof(x) # element type
        n = Int(dim)
        m = (n*(n+1)) ÷ 2
        mat = new{n, m, t}()

        mat.data = zeros(t, m)
        offset = MVector{n,Int}(undef)

        for i in 1:n
            offset[i] = (i*(i-1)) ÷ 2
        end
        mat.offset = SVector(offset)
        mat
    end
end

"""
    BitsUTriMatrix(dim, x)

Like UTriMatrix, but stores data in a static array for speed. The
elements of the matrix must be a bits type.  The entries of the matrix
are not initialized.

An uppeer-triangular matrix of dimension `dim X dim`, whose entries have
the same type as `x`. The parameter `x` is used only for its type. Its
value is not used. The type of x must be a bits type, so
that `isbits(x)` returns `true`.

A matrix is constructed with a call such as

    a = UTriMatrix(4, 1.0)

which generates a 4X4 upper triangular matrix whose entries are of
type Float64. Then, a[3,4] returns the entry at row 3, column 4, and
a[4,3] throws a BoundsError. 
"""
mutable struct BitsUTriMatrix{N,M,T} <: AbstrUTriMatrix
    data :: MVector{M,T}     # elements in matrix
    offset :: SVector{N,Int} # num elements before column j

    function BitsUTriMatrix(dim, x)
        t = typeof(x) # element type
        !isbitstype(t) && throw(DomainError(t,
            "BitsUTriMatrix can only hold bitstype values."))
        n = Int(dim)
        m = (n*(n+1)) ÷ 2
        mat = new{n, m, t}()

        mat.data = MVector{m, t}(undef)
        offset = MVector{n,Int}(undef)

        for i in 1:n
            offset[i] = (i*(i-1)) ÷ 2
        end
        mat.offset = SVector(offset)
        mat
    end
end

@inline function Base.size(A::AbstrUTriMatrix) where {N,M,T}
    dim = length(A.offset)
    (dim, dim)
end

@inline Base.length(A::AbstrUTriMatrix) where {N,M,T} = length(A.data)

@inline Base.IndexStyle(::AbstrUTriMatrix) where {N,M,T} = IndexCartesian()

@inline function Base.checkbounds(A::AbstrUTriMatrix, i::Int,
                          j::Int) where {N,M,T}
    i > j && throw(BoundsError(A, (i,j)))
end

@inline function Base.getindex(A::AbstrUTriMatrix, i::Int,
                               j::Int) where {N,M,T}
    @boundscheck checkbounds(A, i, j)
    A.data[i + A.offset[j]]
end

# Fails if value cannot be converted to type T.
@inline function Base.setindex!(A::AbstrUTriMatrix, value, i::Int,
                                j::Int) where {N,M,T}
    @boundscheck checkbounds(A, i, j)
    A.data[i + A.offset[j]] = value
end

"""
Write a representation of an UTriMatrix to io.
"""
function Base.show(io::IO, A::AbstrUTriMatrix) where {N,M,T}
    n = length(A.offset)
    for i in 1:n
        for j in 1:(i-1)
            print(io, j>1 ? " " : "", zero(typeof(A[1,1])))
        end
        for j in i:n
            print(io, j>1 ? " " : "", A[i,j])
        end
        println(io)
    end
end
