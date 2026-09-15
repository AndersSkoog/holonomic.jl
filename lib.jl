
const ISO3::TSO3 = @SMatrix[1 0 0;0 1 0;0 0 1]
const ISE3::TSE3 = @SMatrix[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]
const Angles :: Vector{TAng} = collect(range(-pi,pi,length=360))
const Zeros :: Vector{TNum} = zeros(TNum,360)
const CR = sin(pi/5)
const RA::TAng = pi/2

Angle(v::TNum) = clamp(v,-pi,pi) :: TAng
Angle(x::TNum,y::TNum) = angle(Complex(x,y)) :: TAng
function Plane3(azi::TAng,polar::TAng) :: TPlane3
  u = TVec3([-sin(azi),cos(azi),0.0])
  v = TVec3([-(sin(polar) * cos(azi)),-(sin(polar) * sin(azi)),cos(polar)])
  return u, v
end

S1(a::TAng) = TS1(cos(a),sin(a))
S1(r::TNum,a::TAng) = TS1(r*cos(a),r*sin(a))

function S2(azi::TAng,polar::TAng) :: TVec3
    pr=sin(polar)
    x=pr*cos(azi)
    y=pr*sin(azi)
    z=cos(polar)
    return TVec3([x,y,z])
end

function S2(r::Float64,azi::TAng,polar::TAng) :: TVec3
  pr=sin(polar)
  x=r*pr*cos(azi)
  y=r*pr*sin(azi)
  z=r*cos(polar)
  return TVec3([x,y,z])
end

function S3(azi::TAng,polar::TAng,orb::TAng) :: TCplx2
  z1=Complex(cos(orb),sin(orb))*cos(azi/2)
  z2=Complex(cos(orb-polar),sin(orb-polar))*sin(azi/2)
  return TCplx2((z1,z2))
end

function S3_R3(p::TCplx2) :: TVec3
  z1,z2 = p[1],p[2]
  d=1-z2.im
  x,y,z=(2*z1.re)/d,(2*z1.im)/d,(2*z2.re)/d
  return TVec3([x,y,z])
end

function S3_R3(azi::TAng,polar::TAng,orb::TAng) :: TVec3
  z1,z2 = S3(azi,polar,orb)
  return S3_R3(TCplx2((z1,z2)))
end

function StereoProj(x::TNum,y::TNum,z::TNum) :: Tupple{ComplexF64,ComplexF64}
  ζ = (x + y * im) / (1 - z)
  ξ = (x - y * im) / (1 + z)
  return ζ,ξ
end

StereoProj(v::TVec3) = StereoProj(v[1],v[2],v[3])

function InvStereoProj(v::TVec2) :: TVec3
  m = abs2(v)W
  x = (2*v[1])/(1+m)
  y = (2*v[2])/(1+m)
  z = (1-m)/(1+m)
  return TVec3([x,y,z])
end

function InvStereoProj(v::TCplx) :: TVec3
    m = abs2(v)
    x = (2*v.re)/(1+m)
    y = (2*v.im)/(1+m)
    z = (1-m)/(1+m)
    return TVec3([x,y,z])
end


function SO3(axis::TVec3,ang::TAng) :: TSO3
  x,y,z=axis[1],axis[2],axis[3]
  c,s=cos(ang),sin(ang)
  d=1-c
  return @SMatrix [
    (c+(x^2)*d)   (((x*y)*d)-(z*s))   ((x*z*d) + (y*s));
    ((y*x*d)+(z*s))  (c + ((y^2)*d))  ((y*z*d) - (x*s));
    ((z*x*d) - (y*s))  ((z*y*d) + (x*s))   (c + (z^2)*d)
  ]
end

function SO3(azi::TAng,polar::TAng,ang::TAng) :: TSO3
  axis=S2(azi,polar)
  return SO3(axis,ang)
end

function MöbiusRotCoef(v::TCplx) :: TMöbiusCoef
  m =abs2(v)
  dn=sqrt(1+m)
  a,b,c=1/dn,v/dn,-conj(v)/dn
  d=a
  return (a,b,c,d)
end

function HolomorphicTransform(dev::TArrCplx,n::Int) :: TArrCplx
  a,b,c,d = MöbiusRotCoef(dev[n])
  return [((a*z)+b)/((c*z)+d) for z in dev]
end

function ConeCircle(azi::TAng,polar::TAng,angles::Vector{TAng}=Angles) :: TArrVec3
  u,v = Plane3(azi,polar)
  c=0.8*S2(azi,polar)
  return [c+CR*cos(t)*u+CR*sin(t)*v for t in angles]
end

function ConeCircle(p::TVec3,angles::TArrAng=Angles) :: TArrVec3
 azi=Angle(p[1],p[2])
 pol=acos(p[3])
 u,v = Plane3(azi,pol)
 c = 0.8*p
 return [c+CR*cos(t)*u+CR*sin(t)*v for t in angles]
end

function TorsionAngle(T1::TVec3,T2::TVec3,T3::TVec3,T4::TVec3) :: TAng
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
  return Angle(x,y)
end

function TorsionAngle(T::TArrVec3,n::Int)::TAng
  return TorsionAngle(T[n],T[n+1],T[n+2],T[n+3])
end


function HopfFibre(azi::TAng,polar::TAng,res::Int=360) :: THopfFibre
  ts=res==360 ? Angles : range(-pi,pi,length=res)
  return THopfFibre([TC2(az,polar,t) for t in ts])
end

HopfFibre(v::TVec3,res::Int=360) = HopfFibre(Angle(v[1],v[2]),acos(v[3]),res) :: THopfFibre
HopfLink(a::TVec3,b::TVec3,res::Int=360) = THopfLink(HopfFibre(a,res),HopfFibre(b,res)) :: THopfLink
HopfLink(a::TVec3,tor::TAng,res::Int=360) = THopfLink(HopfFibre(a,res),HopfFibre(SO3(a,tor)*a,res)) :: THopfLink


function SE3(m::TSO3,p::TVec3) :: TSE3
  return TSE3([m[1,1] m[1,2] m[1,3] p[1]; m[2,1] m[2,2] m[2,3] p[2]; m[3,1] m[3,2] m[3,3] p[3]; 0 0 0 1])
end

CInv(v::ComplexF64) = abs(v) <= 1.0 ? v : 1 / conj(v)

#function fourier_series(a::TArrNum,b::TArrNum) :: TArrAng
#  li = lastindex(a)
#  tau = 2pi
#  if li == 0
#    throw(ErrorException("coef lists cannot be empty"))
#  end
#  if li != lastindex(b)
#    throw(ErrorException("coef lists must be of equal size"))
#  end
  #res = lastindex(angles)
#  out = zeros(TAng,360)
#  for h in 1:li
#    out .+= (a[h] .* cos.(Angles .* h) .+ (b[h] .* sin.(Angles .* h)))
#  end
#  return out
#end

function fourier_series(c1::Vector{Float64},c2::Vector{Float64}) :: Vector{Float64}
  li = lastindex(c1)
  #tau = 2pi
  if li == 0
    throw(ErrorException("coef lists cannot be empty"))
  end
  if li != lastindex(c2)
    throw(ErrorException("coef lists must be of equal size"))
  end
  #res = lastindex(angles)
  out = zeros(TAng,360)
  for h in 1:li
    out .+= (c1[h] .* cos.(Angles .* h) .+ (c2[h] .* sin.(Angles .* h)))
  end
  return out
end


function adjpairs(coll)
  li = lastindex(coll)
  l1 = coll[1:li]
  l2 = append!(coll[2:li],[coll[1]])
  return collect(zip(l1,l2))
end


















