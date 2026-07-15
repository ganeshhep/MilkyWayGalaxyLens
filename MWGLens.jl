include("PlummerLens.jl")
include("HernquistLens.jl")
include("MiyamotoNagaiLens.jl")
include("ExponentialDiskLens.jl")
include("CoredHaloLens.jl")
include("NFWLens.jl")
include("LensUtils.jl")
include("LensPlots.jl")

using .PlummerLens
using .HernquistLens
using .MiyamotoNagaiLens
using .ExponentialDiskLens
using .CoredHaloLens
using .NFWLens
using .LensUtils
using .LensPlots
using LensFactory
using GLMakie
using DelimitedFiles
using YAML
using Latexify
using Base.Threads: @threads

parms = YAML.load_file("MWGParameters.yaml")

data = readdlm("RotationCurveGaiaDR3.txt", skipstart = 1)

# DR3+ data, Labini et al. 2023
R = data[:,1]
Vc = data[:,2]
σVc = data[:,3]

# the lens and source redshifts
zl = parms["Redshift"]["lens"]
zs = parms["Redshift"]["source"]

# initialize cosmology
H0 = parms["Cosmology"]["H0"]
Ωm0 = parms["Cosmology"]["Omegam0"]
Ωr0 = parms["Cosmology"]["Omegar0"]
Ωw0 = parms["Cosmology"]["Omegaw0"]

cosmo = Cosmology.init_cosmology(Omega_m0 = Ωm0, Omega_r0 = Ωr0, Omega_w0 = Ωw0)

# the distances are in the units of meters
Dd = Cosmology.angular_diameter_distance(cosmo, 0., zl)
Ds = Cosmology.angular_diameter_distance(cosmo, 0., zs)
Dds = Cosmology.angular_diameter_distance(cosmo, zl, zs)

# critical surface mass density
Σcr = Lenses.get_critical_density(Dd, Dds, Ds; unit = :msun_pc2) # in units of MSun/pc^2
Σcr *= 1.0e6 # in units of MSun/kpc^2

# conversion factors
arc2rad = Constants.ANGLE_ARCSEC # arcseconds to radians
kpc2m = Constants.DIST_KPC # kpc to meters
MSun2kg = Constants.MASS_SUN # MSun to kg

# constants
G = Constants.CONST_G * MSun2kg/kpc2m * 1e-6  # gravitation constant in units of kpc/MSun * (km/s)^2
H = H0/1000.0 # Hubble constant in units of km/s/kpc

# distances in units of kpc
Dd /= kpc2m
Ds /= kpc2m
Dds /= kpc2m

# source position
β = parms["Source"]["coordinates"]
β1 = β[1]
β2 = β[2]

# grid parameters
N = parms["Grid"]["N"] 
θMin = parms["Grid"]["thetaMin"]
θMax = parms["Grid"]["thetaMax"]
h = (θMax - θMin)/(N - 1)

θ1, θ2 = Lenses.get_meshgrid(θMax, θMax, h)

# grid parameters for the inclined disc model
Nd = parms["Grid"]["Nd"]
θdMin = parms["Grid"]["thetadMin"]
θdMax = parms["Grid"]["thetadMax"]
hd = (θdMax - θdMin)/(Nd - 1)

θd1, θd2 = Lenses.get_meshgrid(θdMax, θdMax, hd)

# grid parameters for rotation curve
Nr = parms["Grid"]["Nr"]
rMin = parms["Grid"]["rMin"]
rMax = parms["Grid"]["rMax"] 
r = range(rMin, rMax, Nr)

# ...............................................................................
# model 1
# parameters based on Milky Way galaxy model from Shin and Evans 2006 (https://arxiv.org/abs/astro-ph/0611134).
Mb1 = parms["Lenses"]["Model 1"]["Hernquist Bulge"]["M"]
r0 = parms["Lenses"]["Model 1"]["Hernquist Bulge"]["r0"]
Md = parms["Lenses"]["Model 1"]["Miyamoto Nagai Disk"]["M"]
Ad = parms["Lenses"]["Model 1"]["Miyamoto Nagai Disk"]["A"]
Bd = parms["Lenses"]["Model 1"]["Miyamoto Nagai Disk"]["B"] 
Mc = parms["Lenses"]["Model 1"]["Cored Isothermal Halo"]["M"]
rc = parms["Lenses"]["Model 1"]["Cored Isothermal Halo"]["rc"]
ξ01 = r0 
ρc = Mc/(4π * rc^3)
Mh1 = 4.0 * π * ρc * rc^3

# dimensionless parameters
ad = Ad/ξ01
bd = Bd/ξ01
ρh = ρc * π * rc/Σcr
λ = r0/rc
md = Md/(π * Σcr * ξ01^2)
mb1 = Mb1/(π * Σcr * ξ01^2)

# scale parameter for lensing potential in arcseconds
d1 = (ξ01/Dd)/arc2rad

# deflection angles and jacobians
α1s  = zeros(Float64, N, N)
α2s  = zeros(Float64, N, N)
ψ11s = zeros(Float64, N, N)
ψ12s = zeros(Float64, N, N)
ψ22s = zeros(Float64, N, N)

@threads for idx in CartesianIndices(α1s)
    res1 = HernquistLens.Deflections(mb1, ξ01, θ1[idx], θ2[idx], d1)
    res2 = MiyamotoNagaiLens.Deflections(md, ad, bd, θ1[idx], θ2[idx], d1)
    res3 = CoredHaloLens.Deflections(ρh, λ, θ1[idx], θ2[idx], d1)

    res4 = HernquistLens.Jacobians(mb1, ξ01, θ1[idx], θ2[idx], d1)
    res5 = MiyamotoNagaiLens.Jacobians(md, ad, bd, θ1[idx], θ2[idx], d1)
    res6 = CoredHaloLens.Jacobians(ρh, λ, θ1[idx], θ2[idx], d1)
    
    α1s[idx], α2s[idx] = res1 + res2 + res3
    ψ11s[idx], ψ12s[idx], ψ22s[idx] = res4 + res5 + res6 
end

# deflection angles and Jacobian matrices
α1 = [α1s, α2s]
J1 = [ψ11s, ψ12s, ψ22s]

# critical curves and caustics
figCrit1, axCrit1, figCaus1, axCaus1 = LensPlots.CriticalCaustics(α1, J1, θMin, θMax, N, d1)
axCrit1.title = "Milky Way Galaxy Lens (Model 1) : Critical Curves"
axCaus1.title = "Milky Way Galaxy Lens (Model 1) : Caustic Curves"

# image positions
θ1Img1, θ2Img1 = LensUtils.Images(θMin, θMax, N, α1, β, d1)

# circualr velocities
Vcb1 = @. HernquistLens.Vc(G, Mb1, r0, r)
Vcd1 = @. MiyamotoNagaiLens.Vc(G, Md, Ad, Bd, r)
Vch1 = @. CoredHaloLens.Vc(G, Mh1, rc, r)
Vc1 = @. sqrt(Vcb1^2 + Vcd1^2 + Vch1^2)
# ................................................................................

# ................................................................................
# model 2
# parameters based on the Milky Way galaxy model (first model) from Labini et al. 2023 (https://arxiv.org/abs/2302.01379).
Mb2 = parms["Lenses"]["Model 2"]["Plummer Bulge"]["M"]
a = parms["Lenses"]["Model 2"]["Plummer Bulge"]["a"]
Mtd = parms["Lenses"]["Model 2"]["Miyamoto Nagai Thin Disk"]["M"]
MTd = parms["Lenses"]["Model 2"]["Miyamoto Nagai Thick Disk"]["M"]
Atd = parms["Lenses"]["Model 2"]["Miyamoto Nagai Thin Disk"]["A"]
ATd = parms["Lenses"]["Model 2"]["Miyamoto Nagai Thick Disk"]["A"]
Btd = parms["Lenses"]["Model 2"]["Miyamoto Nagai Thin Disk"]["B"]
BTd = parms["Lenses"]["Model 2"]["Miyamoto Nagai Thick Disk"]["B"]
r200 = parms["Lenses"]["Model 2"]["NFW Halo"]["r200"]
c = parms["Lenses"]["Model 2"]["NFW Halo"]["c"]
ξ02 = a
ρc = 3.0 * H^2/(8.0 * π * G)
ρ0 = 200.0/3.0 * ρc * c^3/(log(1.0 + c) - c/(1.0 + c))

# dimensionless parameters
a0 = a/ξ02
atd = Atd/ξ02
aTd = ATd/ξ02
btd = Btd/ξ02
bTd = BTd/ξ02
mtd = Mtd/(π * Σcr * ξ02^2)
mTd = MTd/(π * Σcr * ξ02^2)
mb2 = Mb2/(π * Σcr * ξ02^2)
rs = r200/c
κs = ρ0 * rs/Σcr

Mh2 = 4.0 * π * ρ0 * rs^3 

# scale parameter for lensing potential in arcseconds
d2 = (ξ02/Dd)/arc2rad

# deflection angles and jacobians
α1l  = zeros(Float64, N, N)
α2l  = zeros(Float64, N, N)
ψ11l = zeros(Float64, N, N)
ψ12l = zeros(Float64, N, N)
ψ22l = zeros(Float64, N, N)

@threads for idx in CartesianIndices(α1l)
    res1 = PlummerLens.Deflections(mb2, a0, θ1[idx], θ2[idx], d2)
    res2 = MiyamotoNagaiLens.Deflections(mtd, atd, btd, θ1[idx], θ2[idx], d2)
    res3 = MiyamotoNagaiLens.Deflections(mTd, aTd, bTd, θ1[idx], θ2[idx], d2)
    res4 = NFWLens.Deflections(κs, θ1[idx], θ2[idx], d2)

    res5 = PlummerLens.Jacobians(mb2, a0, θ1[idx], θ2[idx], d2)
    res6 = MiyamotoNagaiLens.Jacobians(mtd, atd, btd, θ1[idx], θ2[idx], d2)
    res7 = MiyamotoNagaiLens.Jacobians(mTd, aTd, bTd, θ1[idx], θ2[idx], d2)
    res8 = NFWLens.Jacobians(κs, θ1[idx], θ2[idx], d2)
    
    α1l[idx], α2l[idx] = res1 + res2 + res3 + res4
    ψ11l[idx], ψ12l[idx], ψ22l[idx] = res5 + res6 + res7 + res8 
end

# deflection angles and Jacobian matrices
α2 = [α1l, α2l]
J2 = [ψ11l, ψ12l, ψ22l]

# critical curves and caustics
figCrit2, axCrit2, figCaus2, axCaus2 = LensPlots.CriticalCaustics(α2, J2, θMin, θMax, N, d2)
axCrit2.title = "Milky Way Galaxy Lens (Model 2) : Critical Curves"
axCaus2.title = "Milky Way Galaxy Lens (Model 2) : Caustic Curves"

# image positions
θ1Img2, θ2Img2 = LensUtils.Images(θMin, θMax, N, α2, β, d2)

# circualr velocities
Vcb2 = @. PlummerLens.Vc(G, Mb2, a, r)
Vctd2 = @. MiyamotoNagaiLens.Vc(G, Mtd, Atd, Btd, r)
VcTd2 = @. MiyamotoNagaiLens.Vc(G, MTd, ATd, BTd, r)
Vcd2 = @. sqrt(Vctd2 + VcTd2)
Vch2 = @. NFWLens.Vc(G, Mh2, rs, r)
Vc2 = @. sqrt(Vcb2^2 + Vcd2^2 + Vch2^2)
# ...............................................................................

# .................................................................................
# model 3
# the inclined Milky Way Galaxy model
# the fixed parametrs, free parameters, best fit parameters, the initial guesses and the minimum reduced chi-square value
ξdt = parms["Lenses"]["Model 3"]["Exponential Thin Disk"]["Rd"]
ξdT = parms["Lenses"]["Model 3"]["Exponential Thick Disk"]["Rd"]
zdt = parms["Lenses"]["Model 3"]["Exponential Thin Disk"]["zd"]
zdT = parms["Lenses"]["Model 3"]["Exponential Thick Disk"]["zd"]
q = 0.5 # the axis ratio
qMin = min(zdt/ξdt, zdT/ξdT) 
qtrunc = trunc(q, digits = 2)
qMintrunc = trunc(qMin, digits = 2)
Σ0tdInitial = 15.0
Σ0TdInitial = 20.0
q0 = [G, Mb2, a, ξdt, ξdT, Mh2, rs]
p0 = [Σ0tdInitial, Σ0TdInitial]

pFit, χ2ν = LensUtils.Reducedχ2Fit(Vc, R, σVc, q0, p0, VcTotal)

Σ0td, Σ0Td = pFit

ξ03 = ξ02
xdt = ξdt/ξ03
xdT = ξdT/ξ03
κ0td = Σ0td/Σcr
κ0Td = Σ0Td/Σcr
d3 = (ξ03/Dd)/arc2rad

# deflection angles and jacobians
α1d  = zeros(Float64, Nd, Nd)
α2d  = zeros(Float64, Nd, Nd)
ψ11d = zeros(Float64, Nd, Nd)
ψ12d = zeros(Float64, Nd, Nd)
ψ22d = zeros(Float64, Nd, Nd)

@threads for idx in CartesianIndices(α1d)
    res1 = DeflecsJacobs(θd1[idx], θd2[idx], qMin, xdt, d3, κ0td; rtol = 1e-5)
    res2 = DeflecsJacobs(θd1[idx], θd2[idx], qMin, xdT, d3, κ0Td; rtol = 1e-5)

    α1d[idx], α2d[idx], ψ11d[idx], ψ12d[idx], ψ22d[idx] = res1 + res2

    res3 = PlummerLens.Deflections(mb2, a0, θd1[idx], θd2[idx], d3)
    res4 = NFWLens.Deflections(κs, θd1[idx], θd2[idx], d3)
    α1r, α2r = res3 + res4
    α1d[idx] += α1r
    α2d[idx] += α2r

    res5 = PlummerLens.Jacobians(mb2, a0, θd1[idx], θd2[idx], d3)
    res6 = NFWLens.Jacobians(κs, θd1[idx], θd2[idx], d3)
    ψ11r, ψ12r, ψ22r = res5 + res6
    ψ11d[idx] += ψ11r
    ψ12d[idx] += ψ12r
    ψ22d[idx] += ψ22r
end

α3 = [α1d, α2d]
J3 = [ψ11d, ψ12d, ψ22d]

# critical curves and caustics
figCrit3, axCrit3, figCaus3, axCaus3 = LensPlots.CriticalCaustics(α3, J3, θdMin, θdMax, Nd, d3)
axCrit3.title = "Milky Way Galaxy Lens (Model 3, axis ratio = $qMintrunc) : Critical Curves"
axCaus3.title = "Milky Way Galaxy Lens (Model 3, axis ratio = $qMintrunc) : Caustic Curves"

# image positions
θ1Img3, θ2Img3 = LensUtils.Images(θdMin, θdMax, Nd, α3, β, d3)

# circular velocities
Vctd3 = @. ExponentialDiskLens.Vc(G, Σ0td, ξdt, r)
VcTd3 = @. ExponentialDiskLens.Vc(G, Σ0Td, ξdT, r)
Vcd3 = @. sqrt(Vctd3 + VcTd3)
Vc3 = @. sqrt(Vcb2^2 + Vcd3^2 + Vch2^2)
# .............................................................................

# source and image positions on the source plane and the lens plane, respectively for model 1 and model 2
GLMakie.scatter!(axCrit1, θ1Img1, θ2Img1, color = :red, marker = :star5, markersize = 10)
GLMakie.scatter!(axCaus1, β1, β2, color = :blue, marker = :star5, markersize = 10)

GLMakie.scatter!(axCrit2, θ1Img2, θ2Img2, color = :red, marker = :star5, markersize = 10)
GLMakie.scatter!(axCaus2, β1, β2, color = :blue, marker = :star5, markersize = 10)

GLMakie.scatter!(axCrit3, θ1Img3, θ2Img3, color = :red, marker = :star5, markersize = 10)
GLMakie.scatter!(axCaus3, β1, β2, color = :blue, marker = :star5, markersize = 10)

# rotation curve figures for model 1 and model 2
figVcb = GLMakie.Figure(size = (800, 600))
axVcb = GLMakie.Axis(figVcb[1, 1], xlabel = L"r \ [kpc]", ylabel = L"V_{c} \ [km/s]", title = "Milky Way Galaxy Lens : Bulge Rotation Curve", xgridvisible = false, ygridvisible = false)

figVcd = GLMakie.Figure(size = (800, 600))
axVcd = GLMakie.Axis(figVcd[1, 1], xlabel = L"r \ [kpc]", ylabel = L"V_{c} \ [km/s]", title = "Milky Way Galaxy Lens : Disk Rotation Curve", xgridvisible = false, ygridvisible = false)

figVch = GLMakie.Figure(size = (800, 600))
axVch = GLMakie.Axis(figVch[1, 1], xlabel = L"r \ [kpc]", ylabel = L"V_{c} \ [km/s]", title = "Milky Way Galaxy Lens : Dark Halo Rotation Curve", xgridvisible = false, ygridvisible = false)

figVc = GLMakie.Figure(size = (800, 600))
axVc = GLMakie.Axis(figVc[1, 1], xlabel = L"r \ [kpc]", ylabel = L"V_{c} \ [km/s]", title = "Milky Way Galaxy Lens : Total Rotation Curve", xgridvisible = false, ygridvisible = false)

GLMakie.lines!(axVcb, r, Vcb1, label = "Model 1")
GLMakie.lines!(axVcd, r, Vcd1, label = "Model 1")
GLMakie.lines!(axVch, r, Vch1, label = "Model 1")
GLMakie.lines!(axVc, r, Vc1, label = "Model 1")

GLMakie.lines!(axVcb, r, Vcb2, label = "Model 2")
GLMakie.lines!(axVcd, r, Vcd2, label = "Model 2")
GLMakie.lines!(axVch, r, Vch2, label = "Model 2")
GLMakie.lines!(axVc, r, Vc2, label = "Model 2")

GLMakie.lines!(axVcd, r, Vcd3, label = "Model 3")
GLMakie.lines!(axVc, r, Vc3, label = "Model 3")

# DR3+ rotation curve and error bars
GLMakie.scatter!(axVc, R, Vc, color = :black, marker = :circle, markersize = 5, label = "DR3+ Data")
errorbars!(axVc, R, Vc, σVc, color = :red, label = "Error Bar")

# legends for figures
axislegend(axVcb, position = :rb)  
axislegend(axVcd, position = :rb)
axislegend(axVch, position = :rb)
axislegend(axVc, position = :rb)

# display figures
display(GLMakie.Screen(), figCrit1)
display(GLMakie.Screen(), figCaus1)
display(GLMakie.Screen(), figCrit2)
display(GLMakie.Screen(), figCaus2)
display(GLMakie.Screen(), figCrit3)
display(GLMakie.Screen(), figCaus3)
display(GLMakie.Screen(), figVcb)
display(GLMakie.Screen(), figVcd)
display(GLMakie.Screen(), figVch)
display(GLMakie.Screen(), figVc)

readline()
