module Projection

export CP1, Spinor, SU2, spinor, project, I

struct CP1
    z::ComplexF64
end

struct CP1Atlas
    ζ::ComplexF64
    ξ::ComplexF64
end

# -----------------------
# Spinor = point on S³ ⊂ C²
# -----------------------
struct Spinor
    α::ComplexF64
    β::ComplexF64
end

function Spinor(θ::Float64, φ::Float64)
    α = cos(θ/2)
    β = exp(im * φ) * sin(θ/2)
    Spinor(ComplexF64(α), ComplexF64(β))
end

# projection S³ → CP¹
function project(s::Spinor)
    CP1(s.β / s.α)
end

struct SU2
    U::Matrix{ComplexF64}
end

function SU2(α::ComplexF64, β::ComplexF64)
    n = sqrt(abs2(α) + abs2(β))
    SU2([])
end

# identity element
I = SU2(1.0 + 0im, 0.0 + 0im)

function matrix(u::SU2)
    α, β = u.α, u.β
    return [α -conj(β);β conj(α)]
end

#---------------
# Riemann Sphere
#---------------
struct S2
    x::Float
    y::Float
    z::Float
end

function S2(θ::Float,φ::Float)
    θ,φ=mod(θ,2*pi),mod(φ,2*pi)
    x=cos(φ)*cos(θ)
    y=cos(φ)*sin(θ)
    z=sin(φ)
    return S2(x,y,z)
end

function StereoProj(v::S2)
    ζ = (v.x+v.y*im) / (1-v.z)
    ξ = (v.x-v.y*im) / (1+v.z)
    return CP1Atlas(ζ,ξ)
end

function InvStereoProj(v::ComplexF64)
    m = abs2(v)
    x = (2*v.real)/(1+m)
    y = (2*v.imag)/(1+m)
    z = (1-m)/(1+m)
    return S2(x,y,z)
end

function InvStereoProj(v::CP1Atlas)
    return InvStereoProj(v.ζ), ΙnvStereoProj(v.ξ)
end

function MobiusRotCoef(z::ComplexF64)
    norm = sqrt(1 + real(z)^2 + imag(z)^2)
    a = Complex(1/norm,0)
    b = z/norm
    c = conj(z)/norm
    d = Complex(1/norm,0)
    return (a=a,b=b,c=c,d=d)
end

function MobiusTrans(z::Complex,coef::(Complex,Complex,Complex,Complex))
    a,b,c,d = coef
    return z_trans = ((a*z) + b) / ((c*z) + d)
end

end
