module Holomorph
using .Structs:RollConnection,Holomorph2
export MobiusRotCoef,Holomorph2
invariant_scale = tan(acos(4/5)/2)

function MobiusRotCoef
  m = abs2(z)
  a = 1 / sqrt(1 + m)
  b = z / sqrt(1 + m)
  c = -conj(z) / sqrt(1 + m)
  d = a
  return (a=a,b=b,c=c,d=d)
end

function Holomorph2(c::RollConnection,ind::Int64)
  frame=c.framesSU2[ind]
  ref=c.contacts[ind]
  result = map(p-> MobiusRotTrans(ref,frame),c.contacts)
  return Holomorph2(ind,ref,c.torangles[ind],frame,result)
end
end
