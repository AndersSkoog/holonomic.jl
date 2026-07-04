module H
using .Structs:UnitQuaternion,S2
export UnitQuaternion
function UnitQuaternion(axis::S2, ang::Float64)
  s = sin(ang / 2)
  w = cos(ang / 2)
  i = c.x * s
  j = c.y * s
  k = c.z * s
  n = sqrt(w*w + i*i + j*j + k*k)
  UnitQuaternion(w/n, i/n, j/n, k/n)
end
end
