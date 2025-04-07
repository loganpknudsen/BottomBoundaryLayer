using Parameters

@with_kw struct DenseParams
    δ = 0.5
    T = 30

    # Ht: shelf region height, Lt: slopelength, Lb: length of sfc cooling flux
    f1e4N21e5threta029gammau = (f=1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.29,
        const S∞ = (N²*tand(θ)^2)/(f^2),
        γ = (cosd(θ)*(1+S∞*(1-δ)))^(-1)
    )
    f1e4N21e5threta029gammal = (f=1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.29,
        const S∞ = (N²*tand(θ)^2)/(f^2),
        γ = (cosd(θ)*(1+S∞*(1-δ)))^(-1)
    )