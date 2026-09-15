const FILE = "connections.jls"

connections = isfile(FILE) ? deserialize(FILE) : Dict{String,TSConnection}()

function SPt(azi,pol)
  return @SVector [sin(pol)*cos(azi),sin(pol)*sin(azi),cos(pol)]
end


function ΔG(adjcontacts::Tuple{TVec3,TVec3}) :: Tuple{TSO3,TVec3}
  c1,c2 = adjcontacts
  h=normalize(cross(c1,c2))
  psi=acos(dot(c1,c2))
  #print(psi)
  ΔO=SO3(h,psi)
  n=ΔO[3,:]
  #print(n)
  #print(normalize(n))
  #print(psi/normalize(n))
  ΔP=(psi/norm(n))*n
  #print(ΔP)
  return (ΔO,ΔP)
end

function scurve(a1::TArrNum,a2::TArrNum,b1::TArrNum,b2::TArrNum)
  harms = lastindex(a1)
  if harms == 0
    throw(ErrorException("coef list cannot be empty"))
  end
  if !allequal(lastindex,[a1,a2,b2,b2])
    throw(ErrorException("coef lists cannot be empty"))
  end
  azi = fourier_series(a1,a2)
  pol = fourier_series(b1,b2)
  pts = SPt.(azi,pol)
  return TSCurve((a1,a2,b1,b2),azi,pol,pts)
end



function sconnection(sc::TSCurve) :: TSConnection
  li = lastindex(sc.pts)
  initO = ISO3
  initP = @SVector [0.0,0.0,0.0]
  form::Vector{Tuple{TSO3,TVec3}} = ΔG.(adjpairs(sc.pts))
  orients::TArrSO3 = accumulate(*,getindex.(form,1),init=initO)
  pos::TArrVec3 = accumulate(+,getindex.(form,2),init=initP)
  dev::TArrCplx = Complex.(getindex.(pos,1),getindex.(pos,2))
  return TSConnection(sc,form,pos,orients,dev)
end

sconnection(coef::TSCurveCoef) = sconnection(scurve(coef))

function get_sconnection(key::String)
  if haskey(connections,key)
    return connections[key]
  else
    throw(ErrorException("no matching key"))
  end
end



function get_sconnection(key::String, a1::TArrNum,a2::TArrNum,b1::TArrNum,b2::TArrNum)
  if haskey(connections,key)
    return connections[key]
  end
  sc = scurve(a1,a2,b1,b2)
  conn = sconnection(sc)
  #@show Base.summarysize(conn)
  connections[key] = conn
  serialize(FILE,connections)
  return conn
end



















