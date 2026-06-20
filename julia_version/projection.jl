module Projection
using .Primitives: S², C²
export stereographic, inv_stereographic 

    function stereographic(c::S²)
        c = toVec(c)
        x,y,z = c[0],c[1],c[2]
        ζ = (x + y*im) / (1 - z)
        ξ = (x - y*im) / (1 + z)
        return C²(ζ,ξ)
    end

    function inv_stereographic(v::ComplexF64)
      m = abs(v)^2
      x = (2*v.real)/(1+m)
      y = (2*v.imag)/(1+m)
      z = (1-m)/(1+m)
      return [x,y,z]
    end
end






struct StereographicProjection 
    zeta::ComplexF64
    xi::ComplexF64

    function (dir::S²)
        c = toVec(dir)
        zeta = (c.x + c.y*im) / (1 - c.z)
        xi   = (c.x - c.y*im) / (1 + c.z)
        new(zeta, xi)
    end
end




end