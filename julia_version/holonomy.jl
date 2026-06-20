module Holonomy
using LinearAlgebra
using .Primitives:S²,toVec,vec3,vec2,toArray
using .Orientation:SO3
export ArcS², ArcRot, SphereCurveRollData

    function toS²(v::Array{Float})
        x,y,z = v[0],v[1],v[2]
        r = sqrt(x^2+y^2+z^2)
        θ = atan2(y,x)
        φ = acos(z/r)
        return S²(θ,φ)
    end

    struct ArcS²
        from::S²
        to::S²
    end

    function ArcRot(v::ArcS²)
        v1,v1 = toVec(v.from)),toVec(v.to))
        ang = acos(clamp(dot(v1, v2), -1, 1))
        ax  = normalize(cross(v1, v2))
        dir = toS²(ax)
        return SO3(dir,ang)
    end

    function RollTrans(spherecurve::Array{S²}, i::Int, contact::Array{Float,2}, R::SO3)
        r=1.0
        arc = ArcS²(curve[i],curve[i+1])
        ang = acos(clamp(dot(toArray(arc.from), toArray(arc.to)), -1, 1))
        R_inc = rotation(arc)
        R_new = R * R_inc
        n = [R_new[1,3], R_new[2,3], R_new[3,3]]
        d = [0.0,0.0,-1.0]
        move_dir = cross3(n, d)
        mdn = norm3(move_dir)

        if mdn < 1e-12
            disp = [0.0, 0.0, 0.0]
        else
            move_dir = normalize(move_dir)
            disp = [move_dir[0] * r * ang,
                    move_dir[1] * r * ang,
                    move_dir[2] * r * ang]
        end

        contact_new = contact .+ disp

        return contact_new, R_new
    end

    function inital_R(spherecurve::Array{S²}):
        arc = ArcS²(S²(0,pi/2),curve[0])
        return rotation(arc)
    end

    function SphereCurveRollData(spherecurve::Array{S²})
        prev_contact = vec3(0.0, 0.0, 0.0)
        prev_R = initial_R(curve)
        curveC = vec3[] #planar contact points of rolling sphere
        curveO = SO3[]  #moving frame
        curveT = vec3[] #torsion curve  

        for i in 1:length(curve)-2
            contact, R = rolltranslation(spherecurve,i,prev_contact,prev_R)
            tor_z = R.mtx * [0.0, 0.0, 1.0]
            tor_pt = vec3(contact.x,contact.y,tor_z[3])
            push!(curveC, contact)
            push!(curveO, R)
            push!(curveT, tor_pt)
            prev_contact = contact
            prev_R = R
        end

        return curveC, curveO, curveT
    end

    function HolonomicView3(sphere_curve::Array{S²},seli::Int)
      curveC, curveO, curveT = SphereCurveRollData(sphere_curve)
      sel=Complex(curveC[seli].x,curveC[seli].y)
      a,b,c,d = MobiusRotCoef(sel)
      

      
    
    end











end







