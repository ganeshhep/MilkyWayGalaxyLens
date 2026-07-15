module NFWLens

using StaticArrays

export Deflections, Jacobians, Vc

# NFW potential ψ(θ1, θ2) for the halo of the Milky Way galaxy
# the parameter (dimensionless) ks is, ks = ρs * rs/Σcr
# ρs is the characteristic density of the halo in units of MSun/kpc^3, rs is the scale radius of the halo in kpc
# ρs = Mc/(4π * rs^3 * (log(1.0 + λ) - λ/(1.0 + λ))), where Mc is the mass of the halo in units of MSun, λ = ξ0/rs

function χ(θ)
    if θ > 1.0
        arg1 = sqrt(θ^2 - 1.0)
        return atan(arg1)/arg1

    elseif θ < 1.0
        arg2 = sqrt(1.0 - θ^2)
        return atanh(arg2)/arg2
    end
end

function Deflections(κs, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    α = 4.0 * κs * (log(0.5 * θ) + χ(θ))/θ

    α1 = α * θ1/θ
    α2 = α * θ2/θ

    α1 /= d
    α2 /= d

    return [α1, α2]
end

function Jacobians(κs, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    κ = 2.0 * κs * (1.0 - χ(θ))/(θ^2 - 1.0)
    α = 4.0 * κs * (log(0.5 * θ) + χ(θ))/θ

    ψ11 = 2.0 * κ * θ1^2/θ^2 - α * (θ1^2 - θ2^2)/θ^3
    ψ12 = 2.0 * (κ - α/θ) * θ1 * θ2/θ^2
    ψ22 = 2.0 * κ * θ2^2/θ^2 + α * (θ1^2 - θ2^2)/θ^3

    ψ11 /= d^2
    ψ12 /= d^2
    ψ22 /= d^2

    return [ψ11, ψ12, ψ22]    
end

function Vc(G, M, rs, r)
    # function to compute the circular velocity Vc for the NFW halo model
    # M is the mass of the halo in units of MSun, rs is the scale radius of the halo in kpc and r is the radial distance from the center of the halo in kpc

    Vc2 = G * M/r * (log(1.0 + r/rs) - (r/rs)/(1.0 + r/rs))
    Vc = sqrt(Vc2)

    return Vc
end

end