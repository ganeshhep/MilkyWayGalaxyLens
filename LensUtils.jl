module LensUtils

include("PlummerLens.jl")
include("ExponentialDiskLens.jl")
include("NFWLens.jl")

using Optim
using .PlummerLens
using .ExponentialDiskLens
using .NFWLens
using LinearAlgebra
using LensFactory
using Contour: contour, lines, coordinates

export Reducedχ2Fit, VcTotal, Images

function Reducedχ2Fit(VcO, R, error, q0, p0, VcTotal)
    # function to fit the observed rotation curve data VcO with the theoretical rotation curve VcTotal
    # VcO is the observed circular velocity data in units of km/s
    # R is the radial distance from the center of the disc in kpc
    # error is the measurement error in the observed circular velocity data in units of km/s
    # q0 is a vector containing the fixed parameters for the bulge, disk and halo
    # p0 is a vector containing the initial guess for the free parameters of the disk component

    # number of data points and number of free parameters
    N = length(VcO)
    k = length(p0) 
    
    # degrees of freedom
    ν = N - k
    
    # objective function for the reduced chi-square
    function Reducedχ2(p)
        if any(p .<= 0.0)
            return Inf
        end
        
        χ2T = 0.0
        for i in 1:N
            VcC = VcTotal(R[i], q0, p)
            
            residual = VcO[i] - VcC
            χ2T += (residual/error[i])^2
        end
        
        # reduced chi-square
        return χ2T/ν
    end
    
    # optimization using Nelder-Mead
    result = optimize(Reducedχ2, p0, NelderMead())
    
    # the best fit parameters and the minimum reduced chi-square value
    pFit = Optim.minimizer(result)
    χ2νMin = Optim.minimum(result)
    
    return pFit, χ2νMin
end

function VcTotal(R, q0, p)
    # function to compute the total circular velocity Vc for model 3 of the Milky Way galaxy
    # q0 is a vector containing the fixed parameters for the bulge, disk and halo components
    # p is a vector containing the free parameters for the disk components
    # R is the radial distance from the center of the disk in kpc

    r = R
    G, Mp, a, ξdt, ξdT, Mh, rs = q0
    Σ0td, Σ0Td = p
    
    VcBulge = PlummerLens.Vc(G, Mp, a, r)
    Vctd = ExponentialDiskLens.Vc(G, Σ0td, ξdt, R)
    VcTd = ExponentialDiskLens.Vc(G, Σ0Td, ξdT, R)
    VcHalo = NFWLens.Vc(G, Mh, rs, r)
    VcTotal = sqrt(VcBulge^2 + Vctd^2 + VcTd^2 + VcHalo^2)
    
    return VcTotal
end

function Intersections(x1, y1, x2, y2)
    # x1, y1 and x2, y2 are the vectors containing the x and y coordinates of the first curve and the second curve, respectively.
    # The function returns a vector of tuples containing the coordinates of the intersection points between the two curves. 

    """
    # Based on : https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections
    # Python implementation : https://github.com/sukhbinder/intersection/tree/master/intersect
    """
    # Check if curve1 is consistent
    if !(length(x1) > 1) || !(length(y1) > 1) || length(x1) == length(y1) || 
        throw(ArgumentError("Incompatible input axes for input vectors."))
    end
    
    # Check if curve2 is consistent
    if !(length(x2) > 1) || !(length(y2) > 1) || length(x2) == length(y2) || 
        throw(ArgumentError("Incompatible input axes for input vectors."))
    end    

    # Number of line segments in both curves
    n1 = length(x1) - 1
    n2 = length(x2) - 1

    xy1 = hcat(x1, y1)
    xy2 = hcat(x2, y2)

    dxy1 = diff(xy1, dims = 1)
    dxy2 = diff(xy2, dims = 1)

    ijc = Array{Any}(undef, n2)

    minx1 = mvmin(x1)
	maxx1 = mvmax(x1)
	miny1 = mvmin(y1)
	maxy1 = mvmax(y1)

    for k in 1:n2
        k1 = k + 1
        ijc[k] = findall((minx1 .≤ max(x2[k], x2[k1])) .& (maxx1 .≥ min(x2[k], x2[k1])) .&
                         (miny1 .≤ max(y2[k], y2[k1])) .& (maxy1 .≥ min(y2[k], y2[k1])))

        # Second column
        ijc[k] = [[ij[1], k] for ij in ijc[k]]
    end

    ij = vcat(ijc...)
    i = [vec[1] for vec in ij]
    j = [vec[2] for vec in ij]
    
    n = length(i)
    T = zeros(4, n)

    AA = zeros(4, 4, n)
    AA[[1, 2], 3, :] .= - 1
    AA[[3, 4], 4, :] .= - 1 
    AA[[1, 3], 1, :] .= dxy1[i, :]'
    AA[[2, 4], 2, :] .= dxy2[j, :]'
    B = - hcat(x1[i], x2[j], y1[i], y2[j])'

    for k in 1:n        
        try
            T[:, k] .= AA[:, :, k] \ B[:, k]
        catch err        
            T[:, k] .= Inf        
        end
    end

    inRange = ((T[1, :] .≥ 0.0) .& (T[2, :] .≥ 0.0) .& (T[1, :] .≤ 1.0) .& (T[2, :] .≤ 1.0))'
    
    xy0 = []

    for ll in axes(inRange, 2)
        if inRange[ll]
            push!(xy0,(T[3, ll], T[4, ll]))
        end
    end
    
   # Find unique solutions
   xy0Unique::Vector{NTuple{2, Real}} = []
   l::Int = 0

   for i in axes(xy0, 1)
      l = 0
      for j in i + 1:length(xy0)
         norm = sqrt((xy0[i][1] - xy0[j][1])^2 + (xy0[i][2] - xy0[j][2])^2)
         if norm < 1.0e-12
               l += 1
         end
      end
      if l == 0
         push!(xy0Unique, xy0[i])
      end
   end
   return xy0Unique
end

@inline function mvmin(x)
    return minimum.(eachrow(hcat(x[1:end - 1, :], x[2:end, :])))
end

@inline function mvmax(x)
    return maximum.(eachrow(hcat(x[1:end - 1, :], x[2:end, :])))
end

function Images(θMin, θMax, N, α, βSrc)
    # The function returns the image positions for a given source position βSrc in a N-sized grid with limits θMin and θMax, and with deflection angles α1 and α2. 
    # The function uses contour plotting to find the intersection points of the curves defined by the lens equation.

    h = (θMax - θMin)/(N - 1)
    θ1, θ2 = Lenses.get_meshgrid(θMax, θMax, h)

    β1Src = βSrc[1]
    β2Src = βSrc[2]

    α1 = α[1]
    α2 = α[2]

    cont1 = @. θ1 - α1 - β1Src
    cont2 = @. θ2 - α2 - β2Src

    θ1vals = θ2vals = range(θMin, θMax, N)

    c1 = contour(θ1vals, θ2vals, cont1, 0.0)
    c2 = contour(θ1vals, θ2vals, cont2, 0.0)

    θ1Images = Float64[]
    θ2Images = Float64[]

    for line1 in lines(c1)
        for line2 in lines(c2)
            θ1c1, θ2c1 = coordinates(line1)
            θ1c2, θ2c2 = coordinates(line2)

            θs = Intersections(θ1c1, θ2c1, θ1c2, θ2c2)

            l = length(θs)

            for i in 1:l
                θ = θs[i]

                push!(θ1Images, θ[1])
                push!(θ2Images, θ[2])
            end
        end
    end

    return θ1Images, θ2Images
end
    
end