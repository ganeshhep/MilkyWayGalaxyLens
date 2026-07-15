module PlummerLens

export Deflections, Jacobians, Vc

function Deflections(m, a0, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    α = m/a0^2 * θ/(1.0 + θ^2/a0^2)

    α1 = α * θ1/θ
    α2 = α * θ2/θ

    α1 /= d
    α2 /= d

    return [α1, α2]
end

function Jacobians(m, a0, θ1, θ2, d)
    θ = sqrt(θ1^2 + θ2^2)
    θ /= d

    κ = m/a0^2 * 1.0/(1.0 + θ^2/a0^2)^2
    α = m/a0^2 * θ/(1.0 + θ^2/a0^2)

    ψ11 = 2.0 * κ * θ1^2/θ^2 - α * (θ1^2 - θ2^2)/θ^3
    ψ12 = 2.0 * (κ - α/θ) * θ1 * θ2/θ^2
    ψ22 = 2.0 * κ * θ2^2/θ^2 + α * (θ1^2 - θ2^2)/θ^3

    ψ11 /= d^2
    ψ12 /= d^2
    ψ22 /= d^2

    return [ψ11, ψ12, ψ22]    
end

function Vc(G, M, a, r)
    Vc2 = G * M * r^2/(a^2 + r^2)^1.5
    Vc = sqrt(Vc2)

    return Vc
end

end
