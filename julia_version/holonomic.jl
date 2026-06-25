module Projection

export CP1, Spinor, SU2, spinor, project, I

struct S2
    x::Float64
    y::Float64
    z::Float64
end

function theta(v::S2)
    return acos(clamp(v.z, -1.0, 1.0))
end

function phi(v::S2)
    return atan(v.y, v.x)
end

function S2_from_angles(θ::Float,φ::Float)
    x=cos(φ)*cos(θ)
    y=cos(φ)*sin(θ)
    z=sin(φ)
    n = sqrt(x^2 + y^2 + z^2)
    return S2(x/n,y/n,z/n)
end


struct CP1
    z::ComplexF64
end

struct CP1Atlas
    ζ::ComplexF64
    ξ::ComplexF64
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
    n = sqrt(x^2 + y^2 + z^2)
    return S2(x/n,y/n,z/n)
end

# -----------------------
# Spinor = point on S³ ⊂ C²
# -----------------------
struct UnitSpinor
    α::ComplexF64
    β::ComplexF64
    function UnitSpinor(v::S2)
        θ = theta(v)
        φ = phi(v)
        α = cos(θ/2)
        β = exp(im * φ) * sin(θ/2)
        new(ComplexF64(α), ComplexF64(β))
    end
end

function project(s::UnitSpinor)
    CP1(s.β / s.α)
end

struct SU2
    U::Matrix{ComplexF64}
    function SU2(s::UnitSpinor)
        new([s.α s.β;-conj(s.β) conj(s.α)])
    end
end

function MobiusRotCoef(v::S2)
    U = SU2(UnitSpinor(v)).U
    a = U[1,1]
    b = U[1,2]
    c = U[2,1]
    d = U[2,2]
    return (a=a,b=b,c=c,d=d)
end

function MobiusTrans(z::Complex,coef::(Complex,Complex,Complex,Complex))
    a,b,c,d = coef
    return z_trans = ((a*z) + b) / ((c*z) + d)
end

end

