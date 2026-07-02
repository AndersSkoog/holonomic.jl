module S3
using LinearAlgebra
using .CP1
using .S2
using .SU2
using .Holomorph

struct C2
  z1::ComplexF64
  z2::ComplexF64
end

struct S3
  fiber1::Vector{C2}
  fiber2::Vector{C2}
end

basefiber_xy = [C2(cis(a),0+0im) for a in range(0,2*pi,360)]

function hopf_link(h::Holomorph2, ind::Int64)
  S = InvStereoProj(h.refpoint)
  R = SO3_from_axis_angle(S,h.torangle)
  a = InvStereoProj(h.result[ind])
  b = R ⊗ a
  U1 = SU2_from_axis_angle(a,pi)
  U2 = SU2_from_axis_angle(b,pi)
  fib1 = basefiber_xy ⊗ U1
  fib2 = basefiber_xy ⊗ U2
  return S3(fib1,fib2)
end



#ambigious but saved for later reference
function hopf_link(pts::Array{ComplexF64},frames::Array{SU2},ind::Int64)
  F1 = ind > 1 ? frames[ind-1] : init_frame
  F2 = frames[ind]
  K = torsion_su2(F1,F2)
  P=pts[ind] #selected Complex number returned by Mobius Transform of Roll contact points
  S=inv_stereo_proj(P) #P converted to S2
  U1 = SU2(UnitSpinor(S)) # Frame Associated to S2
  U2 = U1 * exp(K)
  fib1 = basefiber_xy * U1
  fib2 = basefiber_xy * U2
  return S3(fib1,fib2)
end

end
