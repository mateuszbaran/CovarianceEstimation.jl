using Distributions
using PDMats
using PosDefManifold

Random.seed!(1234)

@testset "M-estimators" begin
    # create data drawn randomly from a multivariate t-student
    # distribution with `df` degrees of freedom
    # and check how far the estimated shape is different from `trueC`
    n, t, df = 30, 512, 3.0
    trueC = randP(n)
    trueC = trueC / tr(trueC)
    tdist = MvTDist(3.0, zeros(n), PDMat(Matrix(trueC)))
    X = rand(tdist, t)

    # run 100 simulations
    # and check the Fisher distance between true and estimated shape
    t = 1_000
    ntrials = 100
    dscm = Vector{Float64}(undef, ntrials)
    dtme = similar(dscm)
    dnrtme = similar(dscm)

    for i = 1:ntrials
        # println("trial ", i, " of ", ntrials)
        tdist = MvTDist(df, zeros(n), PDMat(Matrix(trueC)))
        X = rand(tdist, t)
        #println(cov(td))
        scm = cov(X')
        scm = scm / tr(scm)
        M = tme(X)
        nrM = nrtme(X)
        dscm[i] = distance(Fisher, Hermitian(trueC), Hermitian(scm))
        dtme[i] = distance(Fisher, Hermitian(trueC), Hermitian(M))
        dnrtme[i] = distance(Fisher, Hermitian(trueC), Hermitian(nrM))
    end
    # results in dB
    println(10 * log10(mean(dscm)))
    println(10 * log10(mean(dtme)))
    println(10 * log10(mean(dnrtme)))
end
