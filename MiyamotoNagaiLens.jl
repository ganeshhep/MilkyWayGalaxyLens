module MiyamotoNagaiLens

using StaticArrays

# Miyamoto-Nagai disk potential ψ(θ1, θ2) for the disk component of the Milky Way galaxy
# the parameter md (dimensionless) is, md = Md/(π * Σcr * ξ0^2)
# Md is the mass of the disk in units of MSun
# the parameter a (dimensionless) controls the θ1-θ2 asymmetry of the surface mass density, a = 0 is the circularly symmetric surface mass density (disk component)
# the parameter b (dimensionless) controls the sharpness of the ridge along the θ1 axis of the surface mass density (disk component)

export Deflections, Jacobians, Vc

function Deflections(m, a, b, θ1, θ2, d)
    θ1 /= d
    θ2 /= d

    term1 = θ1
    term2 = θ1^2 + (a + sqrt(θ2^2 + b^2))^2
    α1 = m * term1/term2
    α1 /= d
     
    term3 = θ2 * (a + sqrt(θ2^2 + b^2))
    term4 = sqrt(θ2^2 + b^2) * (θ1^2 + (a + sqrt(θ2^2 + b^2))^2)
    α2 = m * term3/term4
    α2 /= d

    return [α1, α2]
end

function Jacobians(m, a, b, θ1, θ2, d)
    θ1 /= d
    θ2 /= d

    term1 = (a + sqrt(θ2^2 + b^2))^2 - θ1^2
    term2 = ((a + sqrt(θ2^2 + b^2))^2 + θ1^2)^2
    ψ11 = m * term1/term2
    ψ11 /= d^2

    term3 = - 2.0 * θ1 * θ2 * (a + sqrt(θ2^2 + b^2))
    term4 = sqrt(θ2^2 + b^2) * (θ1^2 + (a + sqrt(θ2^2 + b^2))^2)^2
    ψ12 = m * term3/term4
    ψ12 /= d^2

    term5 = a * sqrt(θ2^2 + b^2) + 2.0 * θ2^2 + b^2
    term6 = (θ1^2 + (a + sqrt(θ2^2 + b^2))^2) * (θ2^2 + b^2)
    term7 = θ2^2 * (a + sqrt(θ2^2 + b^2)) * (4.0 * a * sqrt(θ2^2 + b^2) + 3.0 * (θ2^2 + b^2) + (θ1^2 + a^2))
    term8 = (θ1^2 + (a + sqrt(θ2^2 + b^2))^2)^2 * (θ2^2 + b^2)^1.5
    ψ22 = m * (term5/term6 - term7/term8)
    ψ22 /= d^2

    return [ψ11, ψ12, ψ22]
end

function Vc(G, M, A, B, r)
    # function to compute the circular velocity Vc for the Miyamoto-Nagai disk model
    # M is the mass of the disk in units of MSun, A is the scale radius of the disk in kpc, B is the scale height of the disk in kpc and r is the radial distance from the center of the disk in kpc

    Vc2 = G * M * r^2/(r^2 + (A + B)^2)^1.5
    Vc = sqrt(Vc2)

    return Vc
end

end