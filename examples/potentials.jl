using StaticArrays


function Harmonic(X)
    return 0.5 * X[1]^2
end

function DoubleWell(X)
    return (X[1]^2 -1)^2
end

function Muller(x)

    aa = @SVector [-1, -1, -6.5, 0.7];
    bb = @SVector [0., 0., 11., 0.6];
    cc = @SVector [-10., -10., -6.5, 0.7];
    AA = @SVector [-200., -100., -170., 15.];
    XX = @SVector [1., 0., -0.5, -1.];
    YY = @SVector [0., 0.5, 1.5, 1.];

    return ( AA[1]*exp(aa[1]*(x[1]-XX[1])^2+bb[1]*(x[1]-XX[1])*(x[2]-YY[1])+cc[1]*(x[2]-YY[1])^2)
                 +AA[2]*exp(aa[2]*(x[1]-XX[2])^2+bb[2]*(x[1]-XX[2])*(x[2]-YY[2])+cc[2]*(x[2]-YY[2])^2)
                 +AA[3]*exp(aa[3]*(x[1]-XX[3])^2+bb[3]*(x[1]-XX[3])*(x[2]-YY[3])+cc[3]*(x[2]-YY[3])^2)
                 +AA[4]*exp(aa[4]*(x[1]-XX[4])^2+bb[4]*(x[1]-XX[4])*(x[2]-YY[4])+cc[4]*(x[2]-YY[4])^2));
end
