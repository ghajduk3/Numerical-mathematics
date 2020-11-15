import Base:size,getindex, setindex!, firstindex, lastindex,*
using LinearAlgebra
"""
    RazprsenaMatrika(V,I,n)

Stores non-zero elements in the value matrix `V`, indexes of the elements from the original matrix in the index matrix `I`.

... 
# Arguments
- `V`: matrix that stores non-zero elements.
- `I`: matrix in which indexes of the non-zero elements are stored.
- `n`: size of the matrix.
...    
"""
struct RazprsenaMatrika{T}
    V::Array{Array{T,1},1}
    I::Array{Array{Int64,1},1}
    n::Int64
end

"""
    RazprsenaMatrika(RM)

Constructor which transforms and stores input `n x n` sparse-matrix into the structure of type `RazprsenaMatrika`.

... 
# Arguments
- `RM`: Input sparse matrix of size n x n.
... 

# Examples
a = [5,0,2,1,0]
    [0,5,3,2,0]
    [0,0,0,0,0]
    [1,0,6,0,0]
    [0,5,7,0,0]   
```jldoctest
julia> a =  RazprsenaMatrika([[5,0,2,1,0],[0,5,3,2,0],[0,0,0,0,0],[1,0,6,0,0],[0,5,7,0,0]])
RazprsenaMatrika{Array{Array{Int64,1},1}}(Array{Int64,1}[[5, 2, 1], [5, 3, 2], [], [1, 6], [5, 7]], Array{Int64,1}[[1, 3, 4], [2, 3, 4], [], [1, 3], [2, 3]], 5)
```
"""
function RazprsenaMatrika(matrika::Array{Array{T,1},1}) where {T}
    n = length(matrika)
    v = Array{Array{T,1},1}(undef,length(matrika))
    I = Array{Array{Int64,1},1}(undef,length(matrika))
    for i=1:n
        temp_v = []
        temp_i = []
        for j=1:n
            if matrika[i][j]!=0
                push!(temp_v,matrika[i][j])
                push!(temp_i,j)
            end
        end
        v[i] = temp_v
        I[i] = temp_i       
    end 
    return RazprsenaMatrika(v,I,n)
end

export RazprsenaMatrika,sor

"""
    getindex(RM,i,j)

Retrieve the value from the matrix stored at positions `[i,j]`. The syntax `RM[i,j]` is converted by the compiler to `getindex(RM,i,j)`. 

... 
# Arguments
- `RM`: Input matrix of type RazprsenaMatrika.
- `i`: Index of the row.
- `j`: Index of the column.
... 

# Examples
```jldoctest
julia> M[1,3]
2
```
"""
function getindex(RM::RazprsenaMatrika,i::Int64,j::Int64)
    val = findfirst(x->x==j,RM.I[i])
    if val==nothing
        return 0
    else 
        return RM.V[i][val]
    end
end

"""
    setindex!(RM,val,i,j)

Store value `val` within matrix `RM` at positions `[i,j]`. The syntax `RM[i,j]=val` is converted by the compiler to `setindex!(RM,val,i,j)`.  

... 
# Arguments
- `RM`: Input matrix of type RazprsenaMatrika.
- `val`: Value that is to be stored within M.
- `i`: Index of the row.
- `j`: Index of the column.
... 

# Examples
```jldoctest
julia> RM[1,3]=5
5
```
"""
function setindex!(RM::RazprsenaMatrika,val::Int,i::Int,j::Int)
    v = findfirst(x->x==j,RM.I[i])
    if v==nothing
        push!(RM.I[i],j)
        push!(RM.V[i],val)
    else 
        RM.V[i][v] = val
    end
end

"""
    firstindex(RM)

Return the first index of matrix RM of type `RazprsenaMatrika`. The syntax `RM[begin,begin]` is converted by the compiler to `firstindex(RM)`. 

# Examples
```jldoctest
julia> firstindex(RM)
1
```
"""
function firstindex(RM::RazprsenaMatrika)
    return 1
end

"""
    lastindex(RM)

Return the last index of matrix RM of type `RazprsenaMatrika`. The syntax `RM[end,end]` is converted by the compiler to `lastindex(RM)`. 

# Examples
```jldoctest
julia> lastindex(RM)
5
```
"""

function lastindex(RM::RazprsenaMatrika)
    return RM.n
end

"""
    *(RM::RazprsenaMatrika,vec::Vector)

Multiplicate matrix of type RazprsenaMatrika with a vector. Throw an error if matrix and vector are not of compatible sizes. 

... 
# Arguments
- `RM`: Input matrix of type RazprsenaMatrika.
- `vec`: Input vector.
... 

# Examples
```jldoctest
julia> RM * ones(5)
5-element Array{Int64,1}:
  8
 10
  0
  7
 12
```
"""

function *(RM::RazprsenaMatrika,vec::Array{T,1}) where {T}
    if RM.n != length(vec) 
        throw(ArgumentError("Matrix and vector are not of compatible sizes"))
    end
    result = Array{T,1}(undef,RM.n)
    for i=1:RM.n
        result[i] = sum(RM.V[i] .* [vec[x] for x in RM.I[i]])
    end
    return result
end

"""
    sor(RM::RazprsenaMatrika, b::Array{T,1}, x0::Array{T,1}, omega, tol=1e-10,maxit=500)

Solve the system of equations with the successive-over-relaxation algorithm until the treshold is reached.

... 
# Arguments
- `RM`: Input matrix of type RazprsenaMatrika.
- `b`: Solution of the equations.
- `x0`: Initial estimation.
- `omega` : Factor of relaxation.
- `tol` : Stoping criteria/treshold.
- `maxiter` : Maximum number of iterations.
... 

# Examples
```jldoctest
julia> sor(SM2,[-4.6,5.2,2.0,5.5],[1.0,1.0,1.0,1.0],1.5,1e-10)
4-element Array{Float64,1}:
 -1.2096969697029012 
  0.8248484848578334 
  1.0                
  0.33616161615179363
```
"""
function sor(RM::RazprsenaMatrika, b::Array{T,1}, x0::Array{T,1}, omega, tol=1e-10,maxiter=500) where {T}
    if RM.n != length(b) &&  RM.n != length(x0)
        throw(ArgumentError("Matrix and vector are not of compatible sizes"))
    end
    is_diagonally_dominant(RM)
    iterations = 0
    while norm(RM*x0-b,Inf)>tol
        for i=1:length(RM.V)
            privr = sum([RM.V[i][ind]*x0[j] for (ind,j) in enumerate(RM.I[i]) if i!=j])
            x0[i] = (1-omega)*x0[i] + (omega/RM[i,i])*(b[i]-privr)
        end
        iterations +=1 
        if iterations > maxiter
            break
        end
    end
    return x0,iterations
end

"""
    Internal hepler function which finds the optimal relaxation factor omega
"""
function optimal_omega(RM::RazprsenaMatrika, b::Array{T,1},tol=1e-10) where {T}
    omega = collect(1:0.01:2)
    iterations = [sor(RM,b,zeros(RM.n),om,tol)[2] for om in omega]
    minimum = argmin(iterations)
    return iterations[minimum],omega[minimum], iterations    
end

"""
Internal helper function. Developed for testing purposes.
Fill the zero matrix with values. Return the matrix and a converted matrix of type RazprsenaMatrika. 
"""
function fill_sample_matrix!(A, n)
    N = n^2
    for i=1:n
        A[i][i] = -2*4.0
        for j=1:floor(Int,n*0.3):n
            (i+j <= n) && (A[i][i+j] = 1*j)
            (i-j >= 1) && (A[i][i-j] = 1*(j+1))
        end
    end
    return hcat(A...)',RazprsenaMatrika(A)
end

"""
    Helper function that checks if the matrix is diagonally dominant
"""
function is_diagonally_dominant(RM::RazprsenaMatrika)
    for i=1:RM.n
        privr = 0
        for j=1:length(RM.V[i])
            privr+=abs(RM[i,j])
        end
        if privr - abs(RM[i,i]) > abs(RM[i,i])
            throw(ArgumentError("Matrix is not diagonally dominant"))
        end
    end  

end    