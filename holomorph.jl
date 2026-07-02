module Holomorph
using .SU2
using .S2
using .Connection:RollConnection
export Holomorph2

invariant_scale = tan(acos(4/5)/2)

struct Holomorph2
  index::Int64
  referece_point::ComplexF64
  torsion_angle::Float64
  frame::SU2
  result::Array{ComplexF64}
  function Holomorph2(c::RollConnection,ind::Int64)
    frame=c.framesSU2[ind]
    ref=c.contacts[ind]
    result = map(p-> MobiusRotTrans(ref,frame),c.contacts)
    new(ind,ref,c.torangles[ind],frame,result)
  end
end





end
