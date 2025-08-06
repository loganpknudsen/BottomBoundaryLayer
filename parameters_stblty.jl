using Parameters

f = 1e-4
θ = 0.01

@with_kw struct DenseParams

    S01gammau = (S = 0.1
        N² = (0.1)^2*f^2/theta^2,
        ϕ = pi/2
        δ = 0.5,
        T = 1,
        Nz = 256,
        Lz = 200,
        V∞ = 0.01,
        γ = (1+0.1)^(-1),
    )

end