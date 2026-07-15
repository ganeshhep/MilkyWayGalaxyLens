module NFWLens

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
    Vc2 = G * M/r * (log(1.0 + r/rs) - (r/rs)/(1.0 + r/rs))
    Vc = sqrt(Vc2)

    return Vc
end

end
