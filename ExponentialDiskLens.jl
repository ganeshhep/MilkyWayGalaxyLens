module ExponentialDiskLens

using QuadGK
using StaticArrays
using SpecialFunctions: besseli, besselk

export DeflecsJacobs, Vc

function DeflecsJacobs(θ1, θ2, q, xd, d, κ0; rtol = 1e-5)
    # function to compute the deflection angles and the Jacobians for the exponential disk
    # q is the axis ratio of the disk
    # xd is the scale length of the disk in units of ξ0, xd = ξd/ξ0
    # κ0 is dimensionless central surface mass density, κ0 = Σ0/Σcr

    θ1 /= d
    θ2 /= d

    # the integrand function
    function integrand(u)
        dr = 1.0 - (1.0 - q^2) * u
        ξq2 = u * (θ1^2 + θ2^2/dr)
        ξ = sqrt(ξq2)
        κ = κ0/q * exp(-ξ/xd)
        κ′ = - κ/(2.0 * xd * ξ)
         
        # integrands
        integ1 = κ/dr^0.5
        integ2 = κ/dr^1.5
        integ3 = u * κ′/dr^0.5
        integ4 = u * κ′/dr^1.5
        integ5 = u * κ′/dr^2.5
        
        return SVector(integ1, integ2, integ3, integ4, integ5)
    end

    # adaptive Gauss-Kronrod quadrature
    result, err = quadgk(integrand, 0.0, 1.0, rtol = rtol)
    J0, J1, K0, K1, K2 = result

    # deflection angles
    α1 = q * θ1 * J0
    α2 = q * θ2 * J1

    # jacobians
    ψ11 = 2.0 * q * θ1^2 * K0 + q * J0
    ψ12 = 2.0 * q * θ1 * θ2 * K1
    ψ22 = 2.0 * q * θ2^2 * K2 + q * J1
    
    # the deflection angles in angular coordinates
    α1 /= d
    α2 /= d

    ψ11 /= d^2
    ψ12 /= d^2
    ψ22 /= d^2

    return [α1, α2, ψ11, ψ12, ψ22]
end

function Vc(G, Σ0, ξd, R)
    # function to compute the circular velocity Vc for the exponential disk model
    # Σ0 is the central surface mass density of the disk in units of MSun/kpc^2 and ξd is the scale radius of the disk in units of kpc
    # R is the radial distance from the center of the disk in units of kpc

    y = 0.5 * R/ξd
    Vc2 = 4.0 * π * G * Σ0 * ξd * y^2 * (besseli(0, y) * besselk(0, y) - besseli(1, y) * besselk(1, y))
    Vc = sqrt(Vc2)  

    return Vc    
end

end