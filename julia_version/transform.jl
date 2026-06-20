module Transform
using .Primitives:S², toVec

    function MobiusRotCoef
        (z::Complex)
        norm = sqrt(1 + real(z)^2 + imag(z)^2)
        a = Complex(1/norm,0)
        b = z/norm
        c = conj(z)/norm
        d = Complex(1/norm,0)
        return (a=a,b=b,c=c,d=d)
    end
    
end
