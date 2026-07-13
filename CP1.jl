module CP1
module CP1
using .Structs:S2 CP1Atlas
export StereoProj, InvStereoProj

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

end
