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
    f1e4theta10N21e5delta05Vinf02gammau = (f = 1e-4,
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
    f1e4theta60N21e6delta05Vinf005gammau = (f = 1e-4,
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
        T = 90,
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
        T = 60,
        Nz = 256,
        Lz = 200,
        V∞ = 0.002,
        S∞ = (1e-5*tand(0.009)^2)/((1e-4)^2),
        γ = (cosd(0.009)*(1+(1e-5*tand(0.009)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta250N21e5delta05Vinf02gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.5,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2.5)^2)/((1e-4)^2),
        γ = (cosd(2.5)*(1+1e-5*(tand(2.5)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta250N21e5delta05Vinf02gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.5,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2.5)^2)/((1e-4)^2),
        γ = ((cosd(2.5)*(1+1e-5*(tand(2.5)^2)/((1e-4)^2)*(1-0.5)))^(-1)+(3-1e-5*tand(2.5)^2*(1e-4)^(-2))*(cosd(2.5)*(3*(1+1e-5*tand(2.5)^2*(1e-4)^(-2))-4*0.5*1e-5*tand(2.5)^2*(1e-4)^(-2)))^(-1))/2,
    )
    f1e4theta250N21e5delta05Vinf02gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.5,
        δ = 0.5,
        T = 60,
        Nz = 384,
        Lz = 300,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2.5)^2)/((1e-4)^2),
        γ = (3-1e-5*tand(2.5)^2*(1e-4)^(-2))*(cosd(2.5)*(3*(1+1e-5*tand(2.5)^2*(1e-4)^(-2))-4*0.5*1e-5*tand(2.5)^2*(1e-4)^(-2)))^(-1),
    )
    f1e4theta20N21e5delta05Vinf02gammau = (f = 1e-4,
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
    f1e4theta20N21e5delta05Vinf02gammam = (f = 1e-4,
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
    f1e4theta20N21e5delta05Vinf02gammal = (f = 1e-4,
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
    f1e4theta20N21e5delta025Vinf02gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2,
        δ = 0.25,
        T = 60,
        Nz = 256,
        Lz = 200,
        V∞ = 0.2,
        S∞ = (1e-5*tand(2)^2)/((1e-4)^2),
        γ = (cosd(2)*(1+(1e-5*tand(2)^2)/((1e-4)^2)*(1-0.25)))^(-1),
    )
    f1e4theta20N21e5delta075Vinf02gammau = (f = 1e-4,
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
        γ = ((cosd(0.29)*(1+0.5*1e-5*(tand(0.29))^2*10^(8)))^(-1)+(3-1e-5*tand(0.29)^2*10^(8))*(cosd(0.29)*(3*(1+1e-5*tand(0.29)^2*10^(8))-4*0.5*1e-5*tand(0.29)^2*10^(8)))^(-1))/2,
    )
    f1e4theta30N21e6delta05Vinf005gammau = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 3,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(3)^2)/((1e-4)^2),
        γ = (cosd(3)*(1+(1e-6*tand(3)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta30N21e6delta05Vinf005gammam = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 3,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(3)^2)/((1e-4)^2),
        γ = ((cosd(3)*(1+0.5*1e-6*(tand(3))^2*10^(8)))^(-1)+(3-1e-6*tand(3)^2*10^(8))*(cosd(3)*(3*(1+1e-6*tand(3)^2*10^(8))-4*0.5*1e-6*tand(3)^2*(10)^(8)))^(-1))/2
    )
    f1e4theta30N21e6delta05Vinf005gammal = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 3,
        δ = 0.5,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(3)^2)/((1e-4)^2),
        γ = (3-1e-6*tand(3)^2*10^(8))*(cosd(3)*(3*(1+1e-6*tand(3)^2*10^(8))-4*1e-6*0.5*tand(3)^2*10^(8)))^(-1)
    )
    f1e4theta30N21e6delta025Vinf005gammau = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 3,
        δ = 0.25,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(3)^2)/((1e-4)^2),
        γ = (cosd(3)*(1+(1e-6*tand(3)^2)/((1e-4)^2)*(1-0.25)))^(-1),
    )
    f1e4theta30N21e6delta075Vinf005gammau = (f = 1e-4,
        N² = 1e-6,               # maximum wave velocity
        θ = 3,
        δ = 0.75,
        T = 30,
        Nz = 256,
        Lz = 200,
        V∞ = 0.05,
        S∞ = (1e-6*tand(3)^2)/((1e-4)^2),
        γ = (cosd(3)*(1+(1e-6*tand(3)^2)/((1e-4)^2)*(1-0.75)))^(-1),
    )
    f1e4theta01812N21e5delta05Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.1812,
        δ = 0.5,
        T = 30,
        Nz = 512,
        Lz = 400,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.1812)^2)/((1e-4)^2),
        γ = (cosd(0.1812)*(1+(1e-5*tand(0.1812)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta01812N21e5delta05Vinf01gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.1812,
        δ = 0.5,
        T = 30,
        Nz = 512,
        Lz = 400,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.1812)^2)/((1e-4)^2),
        γ = ((cosd(0.1812)*(1+0.5*1e-5*tand(0.1812)^2*10^(8)))^(-1)+(3-1e-5*tand(0.1812)^2*10^(8))*(cosd(0.1812)*(3*(1+1e-5*tand(0.1812)^2*10^(8))-4*0.5*1e-5*tand(0.1812)^2*(10)^(8)))^(-1))/2,
    ) 
    f1e4theta01812N21e5delta05Vinf01gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.1812,
        δ = 0.5,
        T = 30,
        Nz = 512,
        Lz = 400,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.1812)^2)/((1e-4)^2),
        γ = (3-1e-5*tand(0.1812)^2*10^(8))*(cosd(0.1812)*(3*(1+1e-5*tand(0.1812)^2*10^(8))-4*0.5*1e-5*tand(0.1812)^2*(10)^(8)))^(-1),
    )
    f1e4theta09059N21e5delta05Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.9059,
        δ = 0.5,
        T = 30,
        Nz = 192,
        Lz = 150,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.9059)^2)/((1e-4)^2),
        γ = (cosd(0.9059)*(1+(1e-5*tand(0.9059)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta09059N21e5delta05Vinf01gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.9059,
        δ = 0.5,
        T = 30,
        Nz = 192,
        Lz = 150,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.9059)^2)/((1e-4)^2),
        γ = ((cosd(0.9059)*(1+0.5*1e-5*tand(0.9059)^2*10^(8)))^(-1)+(3-1e-5*tand(0.9059)^2*10^(8))*(cosd(0.9059)*(3*(1+1e-5*tand(0.9059)^2*10^(8))-4*0.5*1e-5*tand(0.9059)^2*(10)^(8)))^(-1))/2,
    ) 
    f1e4theta09059N21e5delta05Vinf01gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.9059,
        δ = 0.5,
        T = 30,
        Nz = 192,
        Lz = 150,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.9059)^2)/((1e-4)^2),
        γ = (3-1e-5*tand(0.9059)^2*10^(8))*(cosd(0.9059)*(3*(1+1e-5*tand(0.9059)^2*10^(8))-4*0.5*1e-5*tand(0.9059)^2*(10)^(8)))^(-1),
    )
    f1e4theta09059N21e5delta025Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.9059,
        δ = 0.25,
        T = 30,
        Nz = 192,
        Lz = 150,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.9059)^2)/((1e-4)^2),
        γ = (cosd(0.9059)*(1+(1e-5*tand(0.9059)^2)/((1e-4)^2)*(1-0.25)))^(-1),
    )
    f1e4theta09059N21e5delta075Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 0.9059,
        δ = 0.75,
        T = 30,
        Nz = 192,
        Lz = 150,
        V∞ = 0.1,
        S∞ = (1e-5*tand(0.9059)^2)/((1e-4)^2),
        γ = (cosd(0.9059)*(1+(1e-5*tand(0.9059)^2)/((1e-4)^2)*(1-0.75)))^(-1),
    )
    f1e4theta181130N21e5delta05Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1.8113,
        δ = 0.5,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(1.8113)^2)/((1e-4)^2),
        γ = (cosd(1.8113)*(1+(1e-5*tand(1.8113)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta181130N21e5delta05Vinf01gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1.8113,
        δ = 0.5,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(1.8113)^2)/((1e-4)^2),
        γ = ((cosd(1.8113)*(1+0.5*1e-5*tand(1.8113)^2*10^(8)))^(-1)+(3-1e-5*tand(1.8113)^2*10^(8))*(cosd(1.8113)*(3*(1+1e-5*tand(1.8113)^2*10^(8))-4*0.5*1e-5*tand(1.8113)^2*(10)^(8)))^(-1))/2,
    ) 
    f1e4theta181130N21e5delta05Vinf01gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1.8113,
        δ = 0.5,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(1.8113)^2)/((1e-4)^2),
        γ = (3-1e-5*tand(1.8113)^2*10^(8))*(cosd(1.8113)*(3*(1+1e-5*tand(1.8113)^2*10^(8))-4*0.5*1e-5*tand(1.8113)^2*(10)^(8)))^(-1),
    )
    f1e4theta181130N21e5delta025Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1.8113,
        δ = 0.25,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(1.8113)^2)/((1e-4)^2),
        γ = (cosd(1.8113)*(1+(1e-5*tand(1.8113)^2)/((1e-4)^2)*(1-0.25)))^(-1),
    )
    f1e4theta181130N21e5delta075Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 1.8113,
        δ = 0.75,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(1.8113)^2)/((1e-4)^2),
        γ = (cosd(1.8113)*(1+(1e-5*tand(1.8113)^2)/((1e-4)^2)*(1-0.75)))^(-1),
    )
    f1e4theta27160N21e5delta05Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.716,
        δ = 0.5,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(2.716)^2)/((1e-4)^2),
        γ = (cosd(2.716)*(1+(1e-5*tand(2.716)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta27160N21e5delta05Vinf01gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.716,
        δ = 0.5,
        T = 60,
        Nz = 192,
        Lz = 150,
        V∞ = 0.1,
        S∞ = (1e-5*tand(2.716)^2)/((1e-4)^2),
        γ = ((cosd(2.716)*(1+0.5*1e-5*tand(2.716)^2*10^(8)))^(-1)+(3-1e-5*tand(2.716)^2*10^(8))*(cosd(2.716)*(3*(1+1e-5*tand(2.716)^2*10^(8))-4*0.5*1e-5*tand(2.716)^2*(10)^(8)))^(-1))/2,
    ) 
    f1e4theta27160N21e5delta05Vinf01gammal = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.716,
        δ = 0.5,
        T = 120,
        Nz = 320,
        Lz = 250,
        V∞ = 0.1,
        S∞ = (1e-5*tand(2.716)^2)/((1e-4)^2),
        γ = (3-1e-5*tand(2.716)^2*10^(8))*(cosd(2.716)*(3*(1+1e-5*tand(2.716)^2*10^(8))-4*0.5*1e-5*tand(2.716)^2*(10)^(8)))^(-1),
    )
    f1e4theta27160N21e5delta025Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.716,
        δ = 0.25,
        T = 120,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(2.716)^2)/((1e-4)^2),
        γ = (cosd(2.716)*(1+(1e-5*tand(2.716)^2)/((1e-4)^2)*(1-0.25)))^(-1),
    )
    f1e4theta27160N21e5delta075Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 2.716,
        δ = 0.75,
        T = 60,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(2.716)^2)/((1e-4)^2),
        γ = (cosd(2.716)*(1+(1e-5*tand(2.716)^2)/((1e-4)^2)*(1-0.75)))^(-1),
    )
    f1e4theta36190N21e5delta05Vinf01gammau = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 3.619,
        δ = 0.5,
        T = 120,
        Nz = 128,
        Lz = 100,
        V∞ = 0.1,
        S∞ = (1e-5*tand(3.619)^2)/((1e-4)^2),
        γ = (cosd(3.619)*(1+(1e-5*tand(3.619)^2)/((1e-4)^2)*(1-0.5)))^(-1),
    )
    f1e4theta36190N21e5delta05Vinf01gammam = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 3.619,
        δ = 0.5,
        T = 180,
        Nz = 192,
        Lz = 250,
        V∞ = 0.1,
        S∞ = (1e-5*tand(3.619)^2)/((1e-4)^2),
        γ = ((cosd(3.619)*(1+0.5*1e-5*tand(3.619)^2*10^(8)))^(-1)+(3-1e-5*tand(3.619)^2*10^(8))*(cosd(3.619)*(3*(1+1e-5*tand(3.619)^2*10^(8))-4*0.5*1e-5*tand(3.619)^2*(10)^(8)))^(-1))/2,
    ) 
    f1e4theta36190N21e5delta05Vinf01gamma005 = (f = 1e-4,
        N² = 1e-5,               # maximum wave velocity
        θ = 3.619,
        δ = 0.5,
        T = 240,
        Nz = 512,
        Lz = 400,
        V∞ = 0.1,
        S∞ = (1e-5*tand(3.619)^2)/((1e-4)^2),
        γ = 0.05 
    )

end