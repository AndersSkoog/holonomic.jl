module SO3
using StaticArrays
using .S2:S2
export SO3, SO3_from_axis_angle

SO3 = SMatrix{3,3,Float64,9}


function SO3_from_axis_angle(axis::S2,ang::Float64)
    x,y,z = axis.x,axis.y,axis.z
    c = cos(ang)
    s = sin(ang)
    C = 1 - c
    R = [
        c + x*x*C   x*y*C - z*s   x*z*C + y*s;
        y*x*C + z*s  c + y*y*C     y*z*C - x*s;
        z*x*C - y*s  z*y*C + x*s   c + z*z*C
        ]
    return SO3(R)
end
