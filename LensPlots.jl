module LensPlots

using Contour: contour, lines, coordinates
using GLMakie
using Interpolations: Gridded, Linear, interpolate
using LaTeXStrings: @L_str

export CriticalCaustics

function CriticalCaustics(α, J, θMin, θMax, N, d)
    # function to compute the critical curves and the caustics for the Milky Way galaxy lens
    ψ1 = d^2 * α[1]
    ψ2 = d^2 * α[2]

    ψ11 = d^2 * J[1]
    ψ12 = d^2 * J[2]
    ψ22 = d^2 * J[3]

    # compute the convergence κ, the shear components γ₁ and γ₂, the total shear γ, and the determinant of the Jacobian matrix det = (1 - κ)^2 - γ^2
    κ = @. 0.5 * (ψ11 + ψ22)
    γ₁ = @. 0.5 * (ψ11 - ψ22)
    γ₂ = @. ψ12
    γ = @. sqrt(γ₁^2 + γ₂^2)
    det = @. (1 - κ)^2 - γ^2

    figCrit = GLMakie.Figure(size = (800, 600))
    axCrit = GLMakie.Axis(figCrit[1, 1], xlabel = L"θ_{1} \ [arcs]", ylabel = L" θ_{2} \ [arcs]", xgridvisible = false, ygridvisible = false,
                          limits = (θMin, θMax, θMin, θMax), aspect = DataAspect()) 

    figCaus = GLMakie.Figure(size = (800, 600))
    axCaus = GLMakie.Axis(figCaus[1, 1], xlabel = L"β_{1} \ [arcs]", ylabel = L"β_{2} \ [arcs]", xgridvisible = false, ygridvisible = false,
                          limits = (θMin, θMax, θMin, θMax), aspect = DataAspect())

    θ1Vals = θ2Vals = range(θMin, θMax, N)

    DetZ = contour(θ1Vals, θ2Vals, det, 0.0)
    
    # bilinear interpolation
    itp1 = interpolate((θ1Vals, θ2Vals), ψ1, Gridded(Linear()))
    itp2 = interpolate((θ1Vals, θ2Vals), ψ2, Gridded(Linear()))

    for line in lines(DetZ)
            θ1, θ2 = coordinates(line) 

            α1 = [itp1(t1, t2) for (t1, t2) in zip(θ1, θ2)]
            α2 = [itp2(t1, t2) for (t1, t2) in zip(θ1, θ2)]
            
            # caustics from the lens equation
            β1 = @. θ1 - α1
            β2 = @. θ2 - α2 

            GLMakie.lines!(axCrit, θ1, θ2)
            GLMakie.lines!(axCaus, β1, β2)
    end

    return figCrit, axCrit, figCaus, axCaus
end

end