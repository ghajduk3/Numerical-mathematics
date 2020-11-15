using BenchmarkTools
using NumMat
using SparseArrays
import LinearAlgebra

"""
    napolni_vzorcno_matriko(A, n)

Napolni vzorčno matriko z vrednostmi.
"""
function napolni_vzorcno_matriko!(A, n)
    N = n^2
    for i=1:n*n
        A[i,i] = -2*n
        for j=1:n
            (i+j <= N) && (A[i, i+j] = 1)
            (i-j >= 1) && (A[i, i-j] = 1)
        end
    end
    return A
end
"""
    sparse_vzorcna_matrika(n)

Vrne vzorčno matriko dimenzije n x n tipa SparseMatrixCSC.
"""
function sparse_vzorcna_matrika(n)
    N = n^2
    P = spzeros(n^2, n^2)
    print(P)
    napolni_vzorcno_matriko!(P, n)
    
    return P
end

"""
    pasovna_nicle(n, s, z)

Vrne ničelno matriko tipa `PasovnaMatrika` dimenzije `n x n`,
ki ima `s` spodnjih obdiagonal in `z` zgornjih obdiagonal.
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
    pasovna_vzorcna_matrika(n)

Vrne vzorčno matriko dimenzije n x n tipa PasovnaMatrika.
"""
function pasovna_vzorcna_matrika(n)
    N = n^2
    P = pasovna_nicle(N, n, n)
    napolni_vzorcno_matriko!(P, n)
    return P
end

a = sparse_vzorcna_matrika(50)

b = pasovna_vzorcna_matrika(50)

benchmark_sparse = @benchmarkable LinearAlgebra.lu(P) setup=(P=a)
benchmark_pasovna = @benchmarkable NumMat.lu(P) setup=(P=b)
ratio(median(run(benchmark_pasovna)), median(run(benchmark_sparse)))