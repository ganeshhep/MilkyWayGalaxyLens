module CoredHaloLens

using StaticArrays

export Deflections, Jacobians, Vc

# cored isothermal halo potential ψ(θ1, θ2) for the halo component of the Milky Way galaxy
# the parameter (dimensionless) ρh is, ρh = ρc * (πrc/Σcr)
# ρc is the core density of the halo in units of MSun/kpc^3, rc is the core radius of the halo in kpc
# ρc = Mc/(4π * rc^3), where Mc is the mass of the halo in units of MSun
# the parameter λ (dimensionless) is, λ = ξ0/rc


function Deflections(ρh, λ, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    α = 2.0 * ρh * θ/(1.0 + sqrt(1.0 + λ^2 * θ^2))

    α1 = α * θ1/θ
    α2 = α * θ2/θ

    α1 /= d
    α2 /= d

    return [α1, α2]
end

function Jacobians(ρh, λ, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    κ = ρh/sqrt(1.0 + λ^2 * θ^2)
    α = 2.0 * ρh * θ/(1.0 + sqrt(1.0 + λ^2 * θ^2))

    ψ11 = 2.0 * κ * θ1^2/θ^2 - α * (θ1^2 - θ2^2)/θ^3
    ψ12 = 2.0 * (κ - α/θ) * θ1 * θ2/θ^2
    ψ22 = 2.0 * κ * θ2^2/θ^2 + α * (θ1^2 - θ2^2)/θ^3

    ψ11 /= d^2
    ψ12 /= d^2
    ψ22 /= d^2

    return [ψ11, ψ12, ψ22]    
end

function Vc(G, M, rc, r)
    # function to compute the circular velocity Vc for the cored isothermal halo model
    # M is the mass of the halo in units of MSun, rc is the core radius of the halo in kpc and r is the radial distance from the center of the halo in kpc

    if r == 0.0
        Vc2 = 0.0

    else
        Vc2 = G * M/rc * (1.0 - atan(r/rc)/(r/rc))
    end

    Vc = sqrt(Vc2)

    return Vc
end

end