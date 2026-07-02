module SU2
using LinearAlgebra
using .Spin
export SU2,SU2_from_UnitSpinor,σy,σz,iσx,iσy,iσz

struct SU2
    U::Matrix{ComplexF64}
end

I   = SU2(ComplexF64[1 0; 0 1])
σx  = SU2(ComplexF64[0 1; 1im 0])
σy  = SU2(ComplexF64[0 -1im; 1im 0])
σz  = SU2(ComplexF64[1 0; 0 -1])
iσx = im * σx.U
iσy = im * σy.U
iσz = im * σz.U

function SU2_from_UnitSpinor(s::UnitSpinor)
    return SU2([s.α s.β;-conj(s.β) conj(s.α)])
end

function SU2_from_axis_angle(axis::S2, angle::Float64)
  x,y,z = axis.x,axis.y,axis.z
  dotsum = x*σx.U + y*σy.U + z*σz.U
  U = cos(angle/2)*I.U - im*sin(angle/2)*dotsum
  return SU2(U)
end

function MobiusRotTrans(z::ComplexF64,coef::SU2)
    a=coef.U[1,1]
    b=coef.U[1,2]
    c=coef.U[2,1]
    d=coef.U[2,2]
    return ((a*z) + b) / ((c*z) + d)
end

#$\textfrak{su}(2) is the liealgebra constructed from smooth adjacent elements on the liegroup SU2$
#the torsion constructor assumes a and b to be adjacent elements from a cumulative rolling connection created by a closed #curve on S2
function torsion_su2(a::SU2,b::SU2)
    Δ = adjoint(a.U) * b.U
    return log(Δ)
end

struct UnitQuaternion
    w::Float
    i::Float
    j::Float
    k::Float

    function UnitQuaternion(dir::S², ang::Float)
        s = sin(ang / 2)
        w = cos(ang / 2)
        c = toVec(dir)

        i = c.x * s
        j = c.y * s
        k = c.z * s

        n = sqrt(w*w + i*i + j*j + k*k)

        new(w/n, i/n, j/n, k/n)
    end
end



end
