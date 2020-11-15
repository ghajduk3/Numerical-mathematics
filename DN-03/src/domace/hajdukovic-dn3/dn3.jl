using LinearAlgebra
"""
    InterpolationPolynom(X,Y,Cx,Cy,W)

Stores the data needed for calculating interpolation with Barycentric fomrula.

... 
# Arguments
- `X`: array that stores function's x coordinates.
- `Y`: array that stores function values.
- `Cx`: array that stores Chebyshev x coordinates.
- `Cy`: array that stores Chebyshev y coordinates.
- `W`: array that stores weights.
... 

"""

struct InterpolationPolynom{T}
    X::Array{T,1}
    Y::Array{T,1}
    Cx::Array{T,1}
    Cy::Array{T,1}
    W::Array{T,1}
end    

export InterpolationPolynom,inter_polynom

"""
    cheb_nodes(a,b,cheb)

Transforms Chebyshev point to any interval [a,b].

... 
# Arguments
- `a`: Lower bound of interval.
- `b`: Upper bound of interval.
- `cheb`: Chebyshev point.
... 
"""

function cheb_nodes(a,b,cheb)
    return 0.5*(a+b) + 0.5*(b-a) * cheb

end

"""
    weight(n)

Returns an array of size n of the Barycentric weights.

... 
# Arguments
- `n`: Number of Chebyshev points - size of the weights array.
... 
"""

function weight(n)
    weights = ones(n+1)
    for i=1:2:n
        weights[i+1] = -1
    end
    weights[1] = 0.5
    weights[end] = ((-1)^(n))*0.5
    return weights
end

"""
    inter_polynom(f,a,b,n)

Calculates the interpolation polynom for the continous function on the interval [a,b] for the n Chebyshev points. 
Interpolation polynom is calculated for the fixed number of base points from interval [a,b].

... 
# Arguments
- `f`: Function that is to be interpolated.
- `a`: Lower bound of interval.
- `b`: Upper bound of interval.
- `n`: Number of Chebyshev points.
...
# Examples
```jldoctest
julia> inter_polynom(x -> sin(x)/x, 0.01, 10, 5)
InterpolationPolynom{Float64}([0.01, 0.5357894736842105, 1.061578947368421, 1.5873684210526315, 2.113157894736842, 2.6389473684210523, 3.1647368421052633, 3.6905263157894734, 4.216315789473684, 4.742105263157894, 5.267894736842105, 5.793684210526316, 6.319473684210527, 6.8452631578947365, 7.371052631578947, 7.8968421052631586, 8.422631578947367, 8.948421052631579, 9.474210526315789, 10.0], [0.9999833334166665, 0.9528370096707686, 0.8224789122394635, 0.629886970616157, 0.4053138318859372, 0.18255207375475702, -0.007312495008156571, -0.1413830462553275, -0.2085845670602122, -0.21078370223560844, -0.16128527031196976, -0.08115485996779152, 0.005741049808744027, 0.07785610389421581, 0.12015089457117747, 0.12651660211568277, 0.10004340481455126, 0.05124308653078304, -0.0052154673795106994, -0.05440211108893698], [10.0, 9.046039886902863, 6.548539886902862, 3.4614601130971376, 0.963960113097138, 0.009999999999999787], [-0.05440211108893698, 0.040874047963127, 0.040047317178941305, -0.09084049373667034, 0.852168395492975, 0.9999833334166665], [0.5, -1.0, 1.0, -1.0, 1.0, -0.5])
```
"""

function inter_polynom(f,a,b,n)
    x_values = collect(LinRange(a,b,50))
    y_values = f.(x_values)
    weights = weight(n)
    cheb_x = [cheb_nodes(a,b,cos(((i-1)*pi) / n)) for i=1:n+1]
    cheb_y = f.(cheb_x)
    return InterpolationPolynom(x_values,y_values,cheb_x,cheb_y,weights)
end

"""
    interpolation_value(inter_pol, x)

Calculates interpolation value of a point x based on the Lagrange Barycentric formula.

... 
# Arguments
- `inter_pol`: Interpolation polynom of type InterpolationPolynom.
- `x`: X coordinate of the point to be interpolated.
...
# Examples
```jldoctest
julia> inter_value(InterpolationPolynom{Float64}([0.01, 0.5357894736842105, 1.061578947368421, 1.5873684210526315, 2.113157894736842, 2.6389473684210523, 3.1647368421052633, 3.6905263157894734, 4.216315789473684, 4.742105263157894, 5.267894736842105, 5.793684210526316, 6.319473684210527, 6.8452631578947365, 7.371052631578947, 7.8968421052631586, 8.422631578947367, 8.948421052631579, 9.474210526315789, 10.0], [0.9999833334166665, 0.9528370096707686, 0.8224789122394635, 0.629886970616157, 0.4053138318859372, 0.18255207375475702, -0.007312495008156571, -0.1413830462553275, -0.2085845670602122, -0.21078370223560844, -0.16128527031196976, -0.08115485996779152, 0.005741049808744027, 0.07785610389421581, 0.12015089457117747, 0.12651660211568277, 0.10004340481455126, 0.05124308653078304, -0.0052154673795106994, -0.05440211108893698], [10.0, 9.046039886902863, 6.548539886902862, 3.4614601130971376, 0.963960113097138, 0.009999999999999787], [-0.05440211108893698, 0.040874047963127, 0.040047317178941305, -0.09084049373667034, 0.852168395492975, 0.9999833334166665], [0.5, -1.0, 1.0, -1.0, 1.0, -0.5])
       ,0.5357894736842105)
0.9753528229939753
```
"""

function inter_value(inter_pol::InterpolationPolynom, x)
    indeks = findfirst(val->val==x, inter_pol.Cx)
    if indeks != nothing
        return inter_pol.Cy[indeks]
    else
        n = length(inter_pol.Cx)
        nom = Array{Float64}(undef,n)
        denom = Array{Float64}(undef,n)
        for i=1:n
            nom[i] = (inter_pol.Cy[i]*inter_pol.W[i]) / (x-inter_pol.Cx[i])
            denom[i] = inter_pol.W[i] / (x-inter_pol.Cx[i])
        end    
        return sum(nom)/sum(denom)
    end        
end


"""
    evaluate_functions(funs,intervals,iter, threshold)

Evaluates lagrange interpolation for a desired functions on a given interval. Interpolation stopping criteria is when error reaches below the threshold.

... 
# Arguments
- `funs`: Array of functions to be interpolated.
- `intervals`: Array of intervals for desired functions.
- `iter`: Maximum polynomial degree.
- `threshold`: Stopping criteria, threshold for the interpolation error.
...
# Examples
Evaluate functions: x->exp((-x)^2) on interval [-1,1], x->sin(x)/x on interval [0,10], x->abs(x^2 - 2*x) on interval [1,3]
```jldoctest
julia>evaluate_functions([x->exp((-x)^2),x->sin(x)/x,x->abs(x^2 - 2*x)],[(-1,1),(0.001,10),(1,3)],1000,10e-6)
Polynomial degree : 9
Polynomial degree : 10
Polynomial degree : 224
```
"""
# # - `funs`: Array which consists of functions to be interpolated.
# funs = [x->exp((-x)^2),x->sin(x)/x,x->abs(x^2 - 2*x)]
# # - `intervals`: Intervals on which the interpolation performed.
# intervals = [(-1,1),(0.001,10),(1,3)]

function evaluate_functions(funs,intervals,iter=1000,threshold=10e-6)
    # funs = [x->exp((-x)^2),x->sin(x)/x,x->abs(x^2 - 2*x)]
    # intervals = [(-1,1),(0.001,10),(1,3)]
    for n=1:length(funs)
        for i=1:iter
            i_pol = inter_polynom(funs[n],intervals[n][1],intervals[n][2],i)
            values = Array{Float64}(undef,50)
            for j=1:50
                values[j] = inter_value(i_pol,i_pol.X[j])
            end
            if  sum(abs.(i_pol.Y.-values))/length(values) < threshold
                println("Polynomial degree : ",i-1)
                break
            end    
        end  
    end    
end




