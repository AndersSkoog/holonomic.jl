module SU2
using StaticArrays
using .C2:C2
using .S2:S2
using .Spin:UnitSpinor
export SU2, su2, I, σx,σy,σz,iσx,iσy,iσz,SU2_from_UnitSpinor,SU2_from_axis_angle,MobiusRotTrans,torsion_su2

const SU2 = SMatrix{2,2,ComplexF64,4}
const su2 = SMatrix{2,2,ComplexF64,4}

Base.:*(U::SU2, p::C2) = begin
    v = U * [p.z1; p.z2]
    C2(v[1], v[2])
end

I   = SU2(ComplexF64[1 0; 0 1])
σx  = SU2(ComplexF64[0 1; 1im 0])
σy  = SU2(ComplexF64[0 -1im; 1im 0])
σz  = SU2(ComplexF64[1 0; 0 -1])
iσx = su2(im * σx)
iσy = su2(im * σy)
iσz = su2(im * σz)

function SU2_from_UnitSpinor(s::UnitSpinor)
    return SU2([s.α s.β;-conj(s.β) conj(s.α)])
end

function SU2_from_axis_angle(axis::S2, angle::Float64)
  x,y,z = axis.x,axis.y,axis.z
  dotsum = σx*x + σy*x + σz*x
  U = cos(angle/2)*I.U - im*sin(angle/2)*dotsum
  return SU2(U)
end

function MobiusRotTrans(z::ComplexF64,coef::SU2)
    a=coef[1,1]
    b=coef[1,2]
    c=coef[2,1]
    d=coef[2,2]
    return ((a*z) + b) / ((c*z) + d)
end

#$\textfrak{su}(2) is the liealgebra constructed from smooth adjacent elements on the liegroup SU2$
#the torsion constructor assumes a and b to be adjacent elements from a cumulative rolling connection created by a closed #curve on S2
function torsion_su2(a::SU2,b::SU2)
    Δ = adjoint(a.U) * b.U
    return su2(log(Δ))
end

end
