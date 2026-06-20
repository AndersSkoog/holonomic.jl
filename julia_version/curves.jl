module Curves
using Random

    function random_closed_sphere_curve(n::Int=360, k::Int=5)
        t = range(θ, 2π, length=n+1)[1:end-1]
        θ = zeros(Float64, n)
        φ = zeros(Float64, n)
        for i in 1:k
        
            Aθ, Bθ = randn(), randn()
            Aφ, Bφ = randn(), randn()
        
            phaseθ = 2π * rand()
            phaseφ = 2π * rand()
        
            θ .+= Aθ .* cos.(i .* t .+ phaseθ) .+
                Bθ .* sin.(i .* t .+ phaseθ)
        
            φ .+= Aφ .* cos.(i .* t .+ phaseφ) .+
                Bφ .* sin.(i .* t .+ phaseφ)
        end
    
        [S²(θ[j],φ[j]) for j in eachindex(θ)]
    end





end