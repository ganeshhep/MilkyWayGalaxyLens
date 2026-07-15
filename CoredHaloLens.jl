module CoredHaloLens

export Deflections, Jacobians, Vc

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
    if r == 0.0
        Vc2 = 0.0

    else
        Vc2 = G * M/rc * (1.0 - atan(r/rc)/(r/rc))
    end

    Vc = sqrt(Vc2)

    return Vc
end

end
