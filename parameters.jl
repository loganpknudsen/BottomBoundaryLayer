using Parameters

@with_kw struct DenseParams
    δ = 0.5
    T = 30

    # Ht: shelf region height, Lt: slopelength, Lb: length of sfc cooling flux
    f1e4N21e5theta029gammau = (f = 1e-4,
        N²= 1e-5,               # maximum wave velocity
        θ = 0.29,
        δ = δ,
        T = T,
        S∞ = (1e-5*tand(0.29)^2)/((1e-4)^2),
        γ = (cosd(0.29)*(1+(1e-5*tand(0.29)^2)/((1e-4)^2)*(1-δ)))^(-1),
    )

    f1e4N21e5theta029gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.29,
        δ = δ,
        T = T,
        S∞ = (1e-5*tand(0.29)^2)/((1e-4)^2),
        γ = (3-(1e-5*tand(0.29)^2)/((1e-4)^2))*(cosd(0.29)*(3*(1+(1e-5*tand(0.29)^2)/((1e-4)^2))-2*(1e-5*tand(0.29)^2)/((1e-4)^2)))^(-1),
    )
end