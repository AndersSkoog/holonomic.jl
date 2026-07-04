module SO3
using .Structs:SO3,S2
export SO3, SO3_from_axis_angle

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
end
