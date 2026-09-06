
const default_angles = range(0,2pi,length=360)
const CR = sin(pi/5)

function Plane(polar::Float64,azi::Float64)
  u = @SVector [-sin(polar),cos(polar),0]
  v = @SVector [-(sin(azi) * cos(polar)),-(sin(azi) * sin(polar)),cos(azi))
  return u, v
end

function S1(polar::Float64)
  return @SVector [cos(polar),sin(polar)]
end

function S1(r::Float64,ang::Float64)
    x=r*cos(ang)
    y=r*sin(ang)
    return @SVector [x,y]
end

function S2(azi::Float64,polar::Float64)
    pr=sin(polar)
    x=pr*cos(azi)
    y=pr*sin(azi)
    z=cos(polar)
    return @SVector [x,y,z]
end

function S2(r::Float64,azi::Float64,polar::Float64)
  pr=sin(polar)
  x=r*pr*cos(azi)
  y=r*pr*sin(azi)
  z=r*cos(polar)
  return @Svector [x,y,z]
end

function S3(azi::Float64,polar::Float64,orb::Float64)
  z1=Complex(cos(orb),sin(orb))*cos(azi/2)
  z2=Complex(cos(orb-polar),sin(orb-polar))*sin(azi/2)
  return z1,z2
end

function S3_R3(z1::ComplexF64,z2::ComplexF64)
  d=1-z2.im
  x,y,z=(2*z1.re)/d,(2*z1.im)/d,(2*z2.re)/d
  return @SVector [x,y,z]
end

function S3_R3(azi::Float65,polar::Float64,orb::Float64)
  z1,z2 = S3(azi,polar,orb)
  return S3_R_3(z1,z2)
end

function StereoProj(x::Float64,y::Float64,z::Float64)
  ζ = (x + y * im) / (1 - z)
  ξ = (x - y * im) / (1 + z)
  return ζ,ξ
end

StereoProj(v::S2) = StereoProj(v[1],v[2],v[3])

function InvStereoProj(v::ComplexF64) :: SVector{3,Float64}
    m = abs2(v)
    x = (2*v.re)/(1+m)
    y = (2*v.im)/(1+m)
    z = (1-m)/(1+m)
    return @SVector [x,y,z]
end

function HopfFibre(azi::Float64,polar::Float64,res::Int=360)
  ts=res==360 ? default_angles : range(0,2π,length=res)
  z1=[Complex(cos(t),sin(t))*cos(azi/2) for t in ts]
  z2=[Complex(cos((t)-polar),sin((t)-polar))*sin(azi/2) for t in ts]
  return z1,z2
end

function HopfFibre(v::S2,res::Int=360)
  x,y,z=v[1],v[2],v[3]
  azi=atan(y,x)
  polar=acos(z)
  return HopfFibre(azi,polar,res)
end


function SO3(x::Float64,y::Float64,z::Float64,ang::Float64)
 c,s=cos(ang),sin(ang)
 d=1-c
 return @SMatrix [(c+x^2*d) (x*y*d)-(z*s) (x*z*d)+(y*s);(y*x*d)+(z*s) c+(y^2*d) (y*z*d)-(x*s)]
end

function SO3(axis::S2,ang::Float64)
  return SO3(axis[1],axis[2],axis[3],ang)
end

function SO3(azi::Float64,polar::Float64,ang::Float64)
  axis=S2(azi,polar)
  return SO3(axis[1],axis[2],axis[3],ang)
end

function MöbiusRotCoef(v::ComplexF64)
  m =abs2(v)
  dn=sqrt(1+m)
  a,b,c=1/dn,v/dn,-conj(v)/dn
  d=a
  return a,b,c,d
end

function ConeCircle(azi::Float64,polar::Float64,res::Int = 360)
  ts=res==360 ? default_angles : range(0,2π,length=res)
  u,v = Plane3(azi,polar)
  c=0.8*S2(azi,polar)
  return [c+CR*cos(t)*u+CR*sin(t)*v for t in ts]
end

function TorsionAngle(T1::SVector{3,Float64},T2::SVector{3,Float64},T3::SVector{3,Float64},T4::SVector{3,Float64})
  d1 = T2-T1
  d2 = T3-T2
  d3 = T4-T3
  c1 = cross(d1, d2)
  c2 = cross(d2, d3)
  # Degenerate cases
  if norm(c1) < 1e-12 || norm(c2) < 1e-12
        return 0.0
  end
  n1 = normalize(c1)
  n2 = normalize(c2)
  axis = normalize(d2)
  x = dot(n1, n2)
  y = dot(axis, cross(n1, n2))
  return atan(y, x)
end

function TorsionAngle(T::Vector{R3},n::Int)::Float64
  return TorsionAngle(T[n],T[n+1],T[n+2],T[n+3])
end






