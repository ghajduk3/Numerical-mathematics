# Band-restricted diagonally dominant matrices
Band matrix is a sparse matrix whose non-zero entries are confined to a diagonal band, comprising the main diagonal and zero or more diagonals on either side.
Such a representation allows us to perform more time and memory efficient operations on matrices.

Three datatypes are implemented in order to represent band-restricted matrices:
 - PasovnaMatrika consists of:
     - Array representing main diagonal
     - Array of arrays representing upper bands 
     - Array of arrays representing lower bands
 - ZgornjePasovnaMatrika consists of:
     - Array representing main diagonal
     - Array of arrays representing upper bands 
 - SpodnjePasovnaMatrika consists of:
     - Array representing main diagonal
     - Array of arrays representing lower bands 
#### Examples and usage:
    PM = [3 2 5 1]
         [6 6 15 3]
         [0 4 13 1]
         [0 0 15 5]   
```sh
julia> PM = PasovnaMatrika([3,6,13,5],[[2,15,1],[5,3],[1]],[[6,4,15]])
PasovnaMatrika([3.0, 6.0, 13.0, 5.0], Array{Float64,N} where N[[2.0, 15.0, 1.0], [5.0, 3.0], [1.0]], Array{Float64,N} where N[[6.0, 4.0, 15.0]])
```

___
### Functions implemented
For the above metioned datatypes following functions are implemented :
#### Base.getindex
Overrides the builtin `getindex` method and returns an element on the [i,j] position.
```julia
julia> getindex(PM,3,2)
4.0
julia> PM[3,2]
4.0 
```
#### Base.setindex!
Overrides the builtin `setindex!` method and let us replace an element on the [i,j] position with a new element.
```julia
julia> setindex!(PM,5,3,2)
5.0
julia> PM[3,2]=5
5.0 
```
#### Base.*
Overides the builtin multiplication method and represents multiplication of banded matrix with a vector. It checks if the sizes of matrix and vector are compatible otherwise error is rosen.
```julia
julia> PM*[2,3,4,5]
4-element Array{Float64,1}:
37.0
105.0
69.0
85.0
```
#### Base.\
Overrides the builting method. Method multiplicates matrix with a vector from the "left". It solves the system `A*x=b`.
```julia
julia> PM \\ [1.0,2.0,3.0]
3-element Array{Float64,1}:
-0.027027027027026973
-0.13513513513513498 
0.810810810810811
```

#### LU
Performs `LU` decomposition of a given matrix of type `:PasovnaMatrika`. Decomposition is conducted only if the matrix is diagonally dominant. Decomposes input matrix PM into `U - upper triangular matrix` of type ZgornjePasovnaMatrika
and `L - lower triangular matrix` of type SpodnjePasovnaMatrika. 

```julia
julia> NumMat.lu(PM)
(SpodnjePasovnaMatrika([1.0, 1.0, 1.0], Array{Float64,N} where N[[0.3333333333333333, -0.5714285714285715], [-0.3333333333333333]]),
ZgornjePasovnaMatrika([3.0, -2.3333333333333335, 5.285714285714286], Array{Float64,N} where N[[-2.0, 1.6666666666666667], [1.0]]))
```
--- 
### Function optimization and benchmark results
In theory, optimized operations on band-restricted matrices which hold non-zero elements on the diagonals should have better time and memory performace compared to operations on the regular sparse matrices. Benchmark results on my implementation of `lu` decomposition compared to the `lu` decomposition from Linear algebra package on sparse matrices  show worse results event though function is optimized for band-restricted matrices in a sense that it operates only on the elements stored into the diagonals. The cause of such result might be type of the structure used to store bands, which as a consequence has lower performance of functions `getindex` and `setindex!`. 


```sh
julia> include("test/benchmarks/domaca_1.jl")
BenchmarkTools.TrialRatio: 
  time:             16.51513708174147
  gctime:           41.28103238572986
  memory:           18.479881311515413
  allocs:           517647.55384615384
```


