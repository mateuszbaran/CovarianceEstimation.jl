# see NOTE in reference_ledoitwolf.jl
# this comes from appendix D in http://www.econ.uzh.ch/static/wp/econwp264.pdf
# adapted to Julia and checked for match with Octave.
# NOTE: it was modified to make the case p < n with centering work:
# - lambda[max(1, p-n+2):p] instead of lambda[max(1, p-n+1):p]
# - vcat(dtilde0 * ones(p-n+1,1), dtilde1); instead of vcat(dtilde0 * ones(p-n+1,1), dtilde1);

using LinearAlgebra
using Statistics

function matlab_ledoitwolf_analytical_shrinkage(X)
    n, p = size(X);
    @assert n ≥ 12 "there must be more than 12 samples"
    sample = cov(SimpleCovariance(), X); # note centering is applied here
    lambda, u = eigen(sample)
    perm = sortperm(lambda)
    lambda = lambda[perm]

    u = u[:, perm]
    η = ifelse(p<n, n, n-1)
    lambda = lambda[max(1, p-η+1):p]
    L = repeat(lambda, outer=(1, min(p, η)))

    h = η^(-1/3)
    H = h*L';
    x = (L-L')./H;

    ftilde = (3/4/sqrt(5))*mean(max.(1 .- x.^2 ./5,0)./H, dims=2);
    Hftemp = (-3/10/pi)*x + (3/4/sqrt(5)/pi) * (1 .- x.^2 ./5) .* log.(abs.((sqrt(5) .- x) ./ (sqrt(5) .+ x)));
    Hftemp[abs.(x).==sqrt(5)] .= (-3/10/pi)*x[abs.(x).==sqrt(5)];
    Hftilde = mean(Hftemp./H, dims=2);

    if p <= η
        dtilde = lambda ./ ((pi*(p/n)*lambda.*ftilde).^2 .+ (1 .- (p/n) .- pi*(p/n)*lambda.*Hftilde).^2);
    else
        Hftilde0 = (1/pi)*(3/10/h^2 + 3/4/sqrt(5)/h*(1-1/5/h^2) * log((1+sqrt(5)*h)/(1-sqrt(5)*h))) * mean(1 ./ lambda);
        dtilde0 = 1/(pi*(p-η)/η*Hftilde0);
        dtilde1 = lambda ./ (pi^2*lambda.^2 .* (ftilde.^2 .+ Hftilde.^2));
        dtilde = vcat(dtilde0 * ones(p-η,1), dtilde1);
    end

    u * (dtilde[:] .* u')
end
