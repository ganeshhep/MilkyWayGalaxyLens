module HernquistLens

export Deflections, Jacobians, Vc

function χ(θ)
    if θ > 1.0
        arg1 = sqrt(θ^2 - 1.0)
        return atan(arg1)/arg1

    elseif θ < 1.0
        arg2 = sqrt(1.0 - θ^2)
        return atanh(arg2)/arg2
    end
end

function Deflections(m, ξ0, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    α = m/ξ0 * θ * (1.0 - χ(θ))/(θ^2 - 1.0)

    α1 = α * θ1/θ
    α2 = α * θ2/θ

    α1 /= d
    α2 /= d

    return [α1, α2]
end

function Jacobians(m, ξ0, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    κs = 0.5 * m/ξ0^2
    κ = κs/(θ^2 - 1.0)^2 * (- 3.0 + (2.0 + θ^2) * χ(θ))
    α = m/ξ0 * θ * (1.0 - χ(θ))/(θ^2 - 1.0)

    ψ11 = 2.0 * κ * θ1^2/θ^2 - α * (θ1^2 - θ2^2)/θ^3
    ψ12 = 2.0 * (κ - α/θ) * θ1 * θ2/θ^2
    ψ22 = 2.0 * κ * θ2^2/θ^2 + α * (θ1^2 - θ2^2)/θ^3

    ψ11 /= d^2
    ψ12 /= d^2
    ψ22 /= d^2

    return [ψ11, ψ12, ψ22]    
end

function Vc(G, M, r0, r)
    Vc2 = G * M * r/(r + r0)^2
    Vc = sqrt(Vc2)

    return Vc
end

end
