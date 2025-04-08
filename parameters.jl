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
        V∞ = 0.05,
        γ = (cosd(0.29)*(1+(1e-5*tand(0.29)^2)/((1e-4)^2)*(1-δ)))^(-1),
    )

    f1e4N21e5theta029gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.29,
        δ = δ,
        T = T,
        V∞ = 0.05,
        S∞ = (1e-5*tand(0.29)^2)/((1e-4)^2),
        γ = (3-(1e-5*tand(0.29)^2)/((1e-4)^2))*(cosd(0.29)*(3*(1+(1e-5*tand(0.29)^2)/((1e-4)^2))-2*(1e-5*tand(0.29)^2)/((1e-4)^2)))^(-1),
    )
    f1e4N21e6theta1gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1,
        δ = δ,
        T = T,
        V∞ = 0.2,
        S∞ = (1e-5*tand(1)^2)/((1e-4)^2),
        γ = (cosd(1)*(1+(1e-5*tand(1)^2)/((1e-4)^2)*(1-δ)))^(-1),
    )
    f1e4N21e6theta6gammau = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 6,
        δ = δ,
        T = T,
        V∞ = 0.05,
        S∞ = (1e-6*tand(6)^2)/((1e-4)^2),
        γ = (cosd(6)*(1+(1e-6*tand(6)^2)/((1e-4)^2)*(1-δ)))^(-1),
    )
    f1e4N21e7theta05gammau = (f = 1e-4,
        N² = 1e-7,               # maximum wave velocity
        θ = 0.5,
        δ = δ,
        T = T,
        V∞ = 0.001,
        S∞ = (1e-7*tand(0.5)^2)/((1e-4)^2),
        γ = (cosd(0.5)*(1+(1e-7*tand(0.5)^2)/((1e-4)^2)*(1-δ)))^(-1),
    )
    f1e4N21e5theta0009gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.009,
        δ = δ,
        T = T,
        V∞ = 0.002,
        S∞ = (1e-5*tand(0.009)^2)/((1e-4)^2),
        γ = (cosd(0.009)*(1+(1e-5*tand(0.009)^2)/((1e-4)^2)*(1-δ)))^(-1),
    )
end