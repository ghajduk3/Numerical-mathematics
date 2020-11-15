import Base:size,getindex, setindex!, firstindex, lastindex,*, \ 


"""
    Concrete data type which inherits from AbstractArray and represents band matrix consists from:
    1) array representation of main diagonal
    2) array of arrays which represents upper bands
    3) array of arrays which represents lower bands
    
    PM = [3 2 5 1]
         [6 6 15 3]
         [0 4 13 1]
         [0 0 15 5]    

    Usage : 
                                        Diagonal    Upper_Band          Lower_Band
            julia> PM = PasovnaMatrika([3,6,13,5],[[2,15,1],[5,3],[1]],[[6,4,15]])
                PasovnaMatrika([3.0, 6.0, 13.0, 5.0], Array{Float64,N} where N[[2.0, 15.0, 1.0], [5.0, 3.0], [1.0]], Array{Float64,N} where N[[6.0, 4.0, 15.0]])

"""
struct PasovnaMatrika 
    diagonal::Array{Float64}
    upper_band::Array{Array{Float64}}
    lower_band::Array{Array{Float64}}
end

"""
    Concrete data type which inherits from AbstractArray and represents lower triangular band matrix consists from :
        1) array representation of main diagonal
        2) array of arrays which represents lower bands

    Usage:
        julia> low = SpodnjePasovnaMatrika([3,6,13,5],[[6,4,15]])
        SpodnjePasovnaMatrika([3.0, 6.0, 13.0, 5.0], Array{Float64,N} where N[[6.0, 4.0, 15.0]])
"""
struct SpodnjePasovnaMatrika 
    diagonal::Array{Float64}
    lower_band::Array{Array{Float64}}
end

"""
    Concrete data type which inherits from AbstractArray and represents upper triangular band matrix consists from :
        1) array representation of main diagonal
        2) array of arrays which represents upper bands

    Usage:
        julia> upper = ZgornjePasovnaMatrika([3,6,13,5],[[2,15,1],[5,3],[1]])
        ZgornjePasovnaMatrika([3.0, 6.0, 13.0, 5.0], Array{Float64,N} where N[[2.0, 15.0, 1.0], [5.0, 3.0], [1.0]])
"""
struct ZgornjePasovnaMatrika 
    diagonal::Array{Float64}
    upper_band::Array{Array{Float64}}
end

export PasovnaMatrika,SpodnjePasovnaMatrika,ZgornjePasovnaMatrika

"""
    Overrides the builtin getindex method and returns an element on the [i,j] position.
    Accepts input argument of data type PasovnaMatrika.
    PM = [3 2 5 1]
         [6 6 15 3]
         [0 4 13 1]
         [0 0 15 5] 

    Usage:
        julia> getindex(PM,3,2)
        4.0
    or
        julia> PM[3,2]
        4.0        
"""
function getindex(PM::PasovnaMatrika,i::Int64,j::Int64)
    n = length(PM.diagonal)
    if i > n || i < 1 || j > n || j < 1
        throw(BoundsError())
    elseif i<j 
        if (j-i)>length(PM.upper_band)
            return 0
        end
        return PM.upper_band[j-i][i]

    elseif i>j
        if (i-j)>length(PM.lower_band)
            return 0
        end
        return PM.lower_band[i-j][j]
    elseif i==j
        return PM.diagonal[i]
    end
end
"""
Overrides the builtin getindex method and returns an element on the [i,j] position.
Accepts input argument of data type ZgornjePasovnaMatrika.

PM = [3 2 5 1]
     [0 6 15 3]
     [0 0 13 1]
     [0 0 0 5] 

Usage:
    julia> getindex(PM,3,4)
    1.0
or
    julia> PM[3,4]
    1.0     
"""
function getindex(PM::ZgornjePasovnaMatrika,i::Int64,j::Int64)
    n = length(PM.diagonal)
    if i > n || i < 1 || j > n || j < 1 || j < i
        throw(BoundsError())
    elseif i<j 
        if (j-i)>length(PM.upper_band)
            return 0
        end
        return PM.upper_band[j-i][i]
    elseif i==j
        return PM.diagonal[i]
    end
end
"""
Overrides the builtin getindex method and returns an element on the [i,j] position.
Accepts input argument of data type ZgornjePasovnaMatrika.
PM = [3 0 0 0]
     [5 6 0 0]
     [6 2 13 0]
     [0 4 7 5] 

Usage:
    julia> getindex(PM,3,2)
    2.0
or
    julia> PM[3,2]
    2.0     
"""

function getindex(PM::SpodnjePasovnaMatrika,i::Int64,j::Int64)
    n = length(PM.diagonal)
    if i > n || i < 1 || j > n || j < 1 || j > i
        throw(BoundsError())
    elseif i>j
        if (i-j)>length(PM.lower_band)
            return 0
        end
        return PM.lower_band[i-j][j]
    elseif i==j
        return PM.diagonal[i]
    end
end

"""
    Overrides the builtin setindex! method and let us replace an element on the [i,j] position with a new element.
    Accepts input matrix of type PasovnaMatrika.
    PM = [3 2 5 1]
         [6 6 15 3]
         [0 4 13 1]
         [0 0 15 5] 

    Usage:
        julia> setindex!(PM,5,3,2)
        5.0
    or
        julia> PM[3,2]=5
        5.0        
"""

function setindex!(PM::PasovnaMatrika,elem,i::Int64,j::Int64)
    n = length(PM.diagonal)
    if i > n || i < 1 || j > n || j < 1
        throw(BoundsError())
    elseif i<j && (j - i) <= length(PM.upper_band)
            PM.upper_band[j-i][i] = elem
    elseif i>j && (i - j) <= length(PM.lower_band)
            PM.lower_band[(i-j)][j] = elem
    elseif i==j
        PM.diagonal[i] = elem 
    end
    return PM
end

"""
    Overrides the builtin setindex! method and let us replace an element on the [i,j] position with a new element.
    Accepts input matrix of type ZgornjePasovnaMatrika.
    PM = [3 2 5 1]
         [0 6 15 3]
         [0 0 13 1]
         [0 0 0 5] 

    Usage:
        julia> setindex!(PM,5,3,4)
        5.0
    or
        julia> PM[3,4]=5
        5.0        
"""

function setindex!(PM::ZgornjePasovnaMatrika,elem,i::Int64,j::Int64)
    n = length(PM.diagonal)
    if i > n || i < 1 || j > n || j < 1 || j < i
        throw(BoundsError())
    elseif i<j && (j - i) <= length(PM.upper_band)
            PM.upper_band[j-i][i] = elem
    elseif i==j
        PM.diagonal[i] = elem  
    end
    return PM
end

"""
    Overrides the builtin setindex! method and let us replace an element on the [i,j] position with a new element.
    Accepts input matrix of type SpodnjePasovnaMatrika.
    PM = [3 0 0 0]
         [4 6 0 0]
         [6 6 13 0]
         [2 2 7 5] 

    Usage:
        julia> setindex!(PM,5,2,1)
        5.0
    or
        julia> PM[2,1]=5
        5.0        
"""

function setindex!(PM::SpodnjePasovnaMatrika,elem,i::Int64,j::Int64)
    n = length(PM.diagonal)
    if i > n || i < 1 || j > n || j < 1 || j > i
        throw(BoundsError())
    elseif i>j && (i - j) <= length(PM.lower_band)
            PM.lower_band[(i-j)][j] = elem
    elseif i==j
        PM.diagonal[i] = elem         
    end
    return PM
end

"""
    Multiplication of matrix with a vector. It checks if the sizes of matrix and vector are compatible otherwise error is rosen.
    PM = [3 2 5 1]     *  vect = [2]
         [6 6 15 3]              [3]
         [0 4 13 1]              [4]   
         [0 0 15 5]              [5]
    
    Usage:
        julia> PM*[2,3,4,5]
        4-element Array{Float64,1}:
        37.0
        105.0
        69.0
        85.0

"""

function *(PM::PasovnaMatrika,v::Vector)
    n = length(PM.diagonal)
    if length(PM.diagonal) != length(v)
        throw(ArgumentError("Sizes of matrix : ",length(PM.diagonal)," and vecotr ",length(v),"are not compatible"))
    end
    privr = PM.diagonal .* v
    #add upper band
    for i=1:length(PM.upper_band)
        privr[1:end-i] += PM.upper_band[i] .* v[i+1:end]
    end
    #add lower band
    for i=1:length(PM.lower_band)
        privr[i+1:end] += PM.lower_band[i] .* v[1:end-i]
    end
    return privr
end

"""
    Solves the system Ax=b, it multiplicates A*b -- A is inversed.
    The input matrix is decomposed into L and U matrix, 
    then the direct substitution is applied onto L in order to compute y from  L*y=b.

    Reverse substitution is applied onto U to solve U*x=y, returning x which is result.

    PM = [3 -2 1]    b= [1]
         [1 -3 2]       [2]
         [-1 2 4]       [3]

    Usage:
        julia> PM \\ [1.0,2.0,3.0]
        3-element Array{Float64,1}:
        -0.027027027027026973
        -0.13513513513513498 
        0.810810810810811
"""
function \(PM::PasovnaMatrika, v::Vector)
    L,U = lu(PM)
    vec = deepcopy(v)
    n=length(vec)
    vec[1] = vec[1]/L[1,1]
    #Direct substitution
    for i=2:n   
        for j=max(1,i-length(PM.lower_band)):(i-1)
            vec[i] = vec[i]-L[i,j]*vec[j]
        end
        vec[i]=vec[i]/L[i,i]
    end
    # Reverse substitution
    vec[n] = vec[n] / U[n,n]
    for i=n-1:-1:1
        for j=i+1:min(n,length(PM.upper_band)+i)
            vec[i]=vec[i] - U[i,j] * vec[j]
        end
        vec[i] = vec[i] / U[i,i]
    end
    return vec
end

"""
    Performs LU decomposition of a given matrix of type :PasovnaMatrika.
    Decomposition is conducted only if the matrix is diagonally dominant.
    Decomposes input matrix PM into U - upper triangular matrix of type ZgornjePasovnaMatrika
    and L - lower triangular matrix of type SpodnjePasovnaMatrika.

    Algorithm 2.4,from the book
    
    PM = [3 -2 1]
         [1 -3 2]
         [-1 2 4]

    Usage:
    julia> NumMat.lu(PM)
    (SpodnjePasovnaMatrika([1.0, 1.0, 1.0], Array{Float64,N} where N[[0.3333333333333333, -0.5714285714285715], [-0.3333333333333333]]),
    ZgornjePasovnaMatrika([3.0, -2.3333333333333335, 5.285714285714286], Array{Float64,N} where N[[-2.0, 1.6666666666666667], [1.0]]))

"""

function lu(P::PasovnaMatrika)
    is_diagonally_dominant(P)
    PM=deepcopy(P)
    n = length(PM.diagonal)
    for k=1:n-1
        for i=k+1:min(n,k+length(PM.lower_band))
            PM[i,k] /=  PM[k,k]
            for j= k + 1:min(n,k+length(PM.upper_band))
                PM[i,j] -= PM[i,k] * PM[k,j]
            end
        end
    end
    return SpodnjePasovnaMatrika(ones(n),PM.lower_band),ZgornjePasovnaMatrika(PM.diagonal,PM.upper_band)
end

"""
    Helper function that checks if the matrix is diagonally dominant
"""
function is_diagonally_dominant(PM::PasovnaMatrika)
    n = length(PM.diagonal)
    for i=1:n
        privr = 0
        for j=max(1,i-length(PM.lower_band)):min(n,length(PM.upper_band)+i)
            privr+=abs(PM[i,j])
        end
        if privr - abs(PM[i,i]) > abs(PM[i,i])
            throw(ArgumentError("Matrix is not diagonally dominant"))
        end
    end  
end

"""
    Helper function which returns zero's matrix of type PasovnaMatrika with dimensions 'n x n'
    with s lower diagonals and z upper diagonals.
    Function is developed for testing purposes.
"""
function pasovna_nicle(n, s, z)
    diag = zeros(n)
    zgornje = []
    spodnje = []
    for i=1:z
        push!(zgornje,zeros(n-i))
    end
    for i=1:s
        push!(spodnje,zeros(n-i))
    end
    return PasovnaMatrika(diag,zgornje,spodnje)
end

"""
    Helper function developed for testing purposes.
    Fills up input matrix with values.
    Function is borrowed from assistant's benchmarking example.
"""
function napolni_vzorcno_matriko!(A, n)
    N = n^2
    for i=1:n*n
        A[i,i] = -2*n
        for j=1:n
            (i+j <= n) && (A[i, i+j] = 1)
            (i-j >= 1) && (A[i, i-j] = 1)
        end
    end
    return A
end
