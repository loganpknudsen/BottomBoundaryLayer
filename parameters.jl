using Parameters

@with_kw struct DenseParams

    # Ht: shelf region height, Lt: slopelength, Lb: length of sfc cooling flux
    f1e4theta029N21e5delta05Vinf005gammau = (f = 1e-4,
        N²= 1e-5,               # maximum wave velocity
        θ = 0.29,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        S∞ = (1e-5*tand(0.29)^2)/((1e-4)^2),
        V∞ = 0.05,
        γ = (cosd(0.29)*(1+(1e-5*tand(0.29)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )

    f1e4theta029N21e5delta05Vinf005gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.29,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-5*tand(0.29)^2)/((1e-4)^2),
        γ = (3-(1e-5*tand(0.29)^2)/((1e-4)^2))*(cosd(0.29)*(3*(1+(1e-5*tand(0.29)^2)/((1e-4)^2))-0.5*4*(1e-5*tand(0.29)^2)/((1e-4)^2)))^(-1),
    )
    f1e4theta1N21e5delta05Vinf02gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(1)^2)/((1e-4)^2),
        γ = (cosd(1)*(1+(1e-5*tand(1)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta6N21e6delta05Vinf005gammau = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 6,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(6)^2)/((1e-4)^2),
        γ = (cosd(6)*(1+(1e-6*tand(6)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta05N21e7delta05Vinf0001gammau = (f = 1e-4,
        N² = 1e-7,               # maximum wave velocity
        θ = 0.5,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.001,
        S∞ = (1e-7*tand(0.5)^2)/((1e-4)^2),
        γ = (cosd(0.5)*(1+(1e-7*tand(0.5)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta0009N21e5delta05Vinf0002gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.009,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.002,
        S∞ = (1e-5*tand(0.009)^2)/((1e-4)^2),
        γ = (cosd(0.009)*(1+(1e-5*tand(0.009)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f175e7theta05N21e5delta05Vinf02gammau = (f = 175e-7,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.5,
        δ = 0.5,
        T = 30,
        Nz = 384,
        Lz = 300,
        V∞ = 0.2,
        S∞ = (1e-5*tand(0.5)^2)/((175e-7)^2),
        γ = (cosd(0.5)*(1+1e-5*(tand(0.5)^2)/((175e-7)^2)*(1-0.5)))^(-1),
    )
    f175e7theta05N21e5delta05Vinf02gammam = (f = 175e-7,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.5,
        δ = 0.5,
        T = 30,
        Nz = 384,
        Lz = 300,
        V∞ = 0.2,
        S∞ = (1e-5*tand(0.5)^2)/((175e-7)^2),
        γ = ((cosd(0.5)*(1+1e-5*(tand(0.5)^2)/((175e-7)^2)*(1-0.5)))^(-1)+(3-1e-5*tand(0.5)^2*(175e-7)^(-2))*(cosd(0.5)*(3*(1+1e-5*tand(0.5)^2*(175e-7)^(-2))-4*0.5*1e-5*tand(0.5)^2*(175e-7)^(-2)))^(-1))/2,
    )
    f175e7theta05N21e5delta05Vinf02gammal = (f = 175e-7,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.5,
        δ = 0.5,
        T = 30,
        Nz = 384,
        Lz = 300,
        V∞ = 0.2,
        S∞ = (1e-5*tand(0.5)^2)/((175e-7)^2),
        γ = (3-1e-5*tand(0.5)^2*(175e-7)^(-2))*(cosd(0.5)*(3*(1+1e-5*tand(0.5)^2*(175e-7)^(-2))-4*0.5*1e-5*tand(0.5)^2*(175e-7)^(-2)))^(-1),
    )
    f1e4theta05N21e5delta05Vinf02gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2)^2)/((1e-4)^2),
        γ = (cosd(2)*(1+(1e-5*tand(2)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta05N21e5delta05Vinf02gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2)^2)/((1e-4)^2),
        γ = ((cosd(2)*(1+0.5*1e-5*(tand(2))^2*10^(8)))^(-1)+(3-1e-5*tand(2)^2*10^(8))*(cosd(2)*(3*(1+1e-5*tand(2)^2*10^(-8))-4*1e-5*0.5*tand(2)^2*10^(8)))^(-1))/2,
    )
    f1e4theta05N21e5delta05Vinf02gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2)^2)/((1e-4)^2),
        γ = (3-1e-5*tand(2)^2*10^(8))*(cosd(2)*(3*(1+1e-5*tand(2)^2*10^(8))-4*0.5*tand(2)^2*10^(8)*1e-5))^(-1),
    )
    f1e4theta05N21e5delta025Vinf02gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2,
        δ = 0.25,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2)^2)/((1e-4)^2),
        γ = (cosd(2)*(1+(1e-5*tand(2)^2)/((1e-4)^2)*(1-0.25)))^(-1),
    )
    f1e4theta05N21e5delta075Vinf02gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2,
        δ = 0.75,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2)^2)/((1e-4)^2),
        γ = (cosd(2)*(1+(1e-5*tand(2)^2)/((1e-4)^2)*(1-0.75)))^(-1),
    )
    f1e4theta029N21e5delta05Vinf005gammam = (f = 1e-4,
        N²= 1e-5,               # maximum wave velocity
        θ = 0.29,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        S∞ = (1e-5*tand(0.29)^2)/((1e-4)^2),
        V∞ = 0.05,
        γ = ((cosd(0.29)*(1+0.5*(1e-5*tand(0.29))^2*10^(8)))**(-1)+(3-1e-5*tand(0.29)^2*10^(8))*(cosd(0.29)*(3*(1+1e-5*tand(0.29)^2*10^(8))-4*0.5*1e-5*tand(0.29)^2*10^(8)))^(-1))/2,
    )
    f1e5theta029N21e6delta05Vinf005gammau = (f = 1e-5,
        N² = 1e-6,               # maximum wave velocity
        θ = 0.29,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(0.29)^2)/((1e-5)^2),
        γ = (cosd(0.29)*(1+(1e-6*tand(0.29)^2)/((1e-5)^2)*(1-0.5)))^(-1),
    )
    f1e5theta029N21e6delta05Vinf005gammam = (f = 1e-5,
        N² = 1e-6,               # maximum wave velocity
        θ = 0.29,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(0.29)^2)/((1e-5)^2),
        γ = ((cosd(0.29)*(1+0.5*(1e-6*tand(0.29))^2*10^(10)))^(-1)+(3-1e-6*tand(0.29)^2*10^(10))*(cosd(0.29)*(3*(1+1e-6*tand(0.29)^2*10^(10))-4*0.5*1e-6*tand(0.29)^2*(10)^(10)))^(-1))/2
    )
    f1e5theta029N21e6delta05Vinf005gammal = (f = 1e-5,
        N² = 1e-6,               # maximum wave velocity
        θ = 0.29,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(0.29)^2)/((1e-5)^2),
        γ = (3-1e-6*tand(0.29)^2*10^(10))*(cosd(0.29)*(3*(1+1e-6*tand(0.29)^2*10^(10))-4*1e-6*0.5*tand(0.29)^2*10^(10)))^(-1)
    )
    f1e5theta029N21e6delta025Vinf005gammau = (f = 1e-5,
        N² = 1e-6,               # maximum wave velocity
        θ = 0.29,
        δ = 0.25,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(0.29)^2)/((1e-5)^2),
        γ = (cosd(0.29)*(1+(1e-6*tand(0.29)^2)/((1e-5)^2)*(1-0.25)))^(-1),
    )
    f1e5theta029N21e6delta075Vinf005gammau = (f = 1e-5,
        N² = 1e-6,               # maximum wave velocity
        θ = 0.29,
        δ = 0.75,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(0.29)^2)/((1e-5)^2),
        γ = (cosd(0.29)*(1+(1e-6*tand(0.29)^2)/((1e-5)^2)*(1-0.75)))^(-1),
    )

end