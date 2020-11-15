const g = 9.80665
"""
    runge_kutta4(f,h,t,y)

Estimates the solution of Ordinary Differential Equations at time step `t+1`. Returns the time step `t+1` along with new estimates of `y` at time step `t+1`.   

... 
# Arguments
- `f`: Differential of the model function to be estimated.    
- `h`: step size.
- `t`: current time step.
- `y`: scalar or vector with values of the system of ODE at time step `t`
... 

# Examples
```jldoctest
julia> runge_kutta4((t,(theta,dtheta))->[dtheta, -9.80665*sin(theta)],0.04,2,[pi/6,0])
(2.04, [0.5196805600301858, -0.19568843149086004])
```
"""

function runge_kutta4(f,h,t,y)

    k1 = h*f(t, y)
	k2 = h*f(t+h/2, y + k1/2)
	k3 = h*f(t+h/2, y + k2/2)
    k4 = h*f(t+h, y + k3)

	return t+h,y + (k1 + 2*k2 + 2*k3 + k4)/6    
end  


"""
    pendulum(l,t,theta0,dtheta0,n)

Models the oscilations of mathematical pendulum and returns pendulum's angle offset in radians at time step `t+1`.    
... 
# Arguments
- `l`: length of a string or wire.  
- `t`: time of pendulum's oscilations.
- `theta0`: angle offset of pendulum.
- `dtheta0`: angular velocity.
- `n`: number of time steps to divide interval [0,t] on. 
... 

# Examples
```jldoctest
julia> pendulum(1,10,pi/6,0,100)
0.4212055089473042
```
"""
function pendulum(l,t,theta0,dtheta0,n)
    t1 = 0
    f(t,(theta,dtheta)) = [dtheta, -(g/l)*sin(theta)]
    h  = t/n

    theta = Array{Float64}(undef,n+1)
    dtheta = Array{Float64}(undef,n+1)
    theta[1] = theta0
    dtheta[1] = dtheta0
    for i=1:n
        t1,(theta[i+1],dtheta[i+1]) = runge_kutta4(f,h,t1,[theta[i],dtheta[i]])
    end
    return theta[end],dtheta[end],theta
end

"""
    harmonic_pendulum(l,t,theta0,n)

Models the oscilations of mathematical pendulum at small angles. It simulates the motion as simple harmonic oscilator..    
... 
# Arguments
- `l`: length of a string or wire.  
- `t`: time of pendulum's oscilations.
- `theta0`: angle offset of pendulum.
- `n`: number of time steps to divide interval [0,t] on. 
... 

# Examples
```jldoctest
julia> harmonic_pendulum(1,10,pi/6,100)
0.5209643557637554
```
"""
function harmonic_pendulum(l,t,theta0,n)
    h = t/n
    theta = Array{Float64}(undef,n+1)
    time = [0:h:t;]
    theta = theta0 * cos.(time * sqrt(g/l))
    return theta[end]

end

"""
    period(theta0,dtheta0,t,n)

Calculates the period of harmonic pendulum. First step of calculation is finding period's approximate. 
At second step Newton method is used in order to find better aproximate of the period.      
... 
# Arguments  
- `theta0`: angle offset of pendulum.
- `dtheta0`: angular velocity.
- `t`: time of pendulum's oscilations.
- `n`: number of time steps to divide interval [0,t] on. 
... 

# Examples
```jldoctest
julia> period(pi/3,0,10,100)
2.153241027194634
```
"""
function period(theta0,dtheta0,t,n)
    l=1
    h = t/n
    t0 = 2
    times = Array{Float64}(undef,n)
    difs = Array{Float64}(undef,n)
    (theta,dtheta,allthetas) = pendulum(1,t,theta0,dtheta0,n)
    time = [0:h:t;]
    maks = false
    local_maxims = []
    if allthetas[1]>allthetas[2]
        maks = true
    end   
    for i=2:length(allthetas)-1
        if allthetas[i-1]<allthetas[i]>allthetas[i+1]
            push!(local_maxims,i)
        end
    end
    if maks
        t0 = time[local_maxims[1]]
    else
        t0 = time[local_maxims[2]] - time[local_maxims[1]]
    end
    (theta,dtheta,allthetas) = pendulum(1,t0,theta0,dtheta0,n)
    for j=1:n 
        tn = t0 - (theta-theta0)/dtheta
        t1,(theta,dtheta) = runge_kutta4((t,(x,y))->[y,-9.80665*sin(x)],tn-t0,tn,[theta,dtheta])
        t0 = tn
        difs[j] = abs(theta-theta0)
        times[j] = tn
        if (abs(theta-theta0) < 1e-10)
            return tn
        end
    end
    return times[argmin(difs)]

end

"""
    period_start_angle(theta0,theta1)

Helper function which evaluates periods of the pendulum based on different starting angular displacements. It continously calculate the period on the interval [a,b] with step 0.1.    
"""
function period_start_angle(theta0,theta1)
    
    thetas = [theta0:0.1:theta1;]
    periods = zeros(length(thetas))
    for i=1:length(thetas)-2
        periods[i] = period(thetas[i],0.1,10,100) 
    end
    return periods
end    

"""
    period_start_velocity(dtheta0,dtheta1)

Helper function which evaluates periods of the pendulum based on different starting angular velocities. It continously calculate the period on the interval [a,b] with step 0.1.    
"""    
function period_start_velocity(dtheta0,dtheta1)
    
    dthetas = [dtheta0+0.1:0.1:dtheta1;]
    periods = zeros(length(dthetas))
    
    for i=1:length(dthetas)
        periods[i] = period(0,dthetas[i],10,100) 
    end
    return periods
end    




