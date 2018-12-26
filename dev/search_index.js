var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#CovarianceEstimation.jl-1",
    "page": "Home",
    "title": "CovarianceEstimation.jl",
    "category": "section",
    "text": "Lightweight robust covariance estimation in Julia.A package for estimating covariance matrices."
},

{
    "location": "#Package-Features-1",
    "page": "Home",
    "title": "Package Features",
    "category": "section",
    "text": "Standard corrected and uncorrected covariance estimators,\nLinear and Nonlinear shrinkage estimators\nFocus on speed and lightweight dependencies"
},

{
    "location": "#Manual-outline-1",
    "page": "Home",
    "title": "Manual outline",
    "category": "section",
    "text": "Pages = [\"man/methods.md\", \"man/lshrink.md\", \"man/nlshrink.md\"]"
},

{
    "location": "#Library-Outline-1",
    "page": "Home",
    "title": "Library Outline",
    "category": "section",
    "text": "Pages = [\"lib/public.md\", \"lib/internals.md\"]"
},

{
    "location": "#main-index-1",
    "page": "Home",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"lib/public.md\"]"
},

{
    "location": "man/methods/#",
    "page": "Methods",
    "title": "Methods",
    "category": "page",
    "text": ""
},

{
    "location": "man/methods/#Methods-1",
    "page": "Methods",
    "title": "Methods",
    "category": "section",
    "text": ""
},

{
    "location": "man/methods/#Notations-1",
    "page": "Methods",
    "title": "Notations",
    "category": "section",
    "text": "In these docs, X denotes an ntimes p data matrix describing n observations with p features (variables) and possibly p  n. For most estimators, the matrix is assumed to have entries in mathbb R.note: Note\nWhile in the docs we assume that the rows of X correspond to observations and the columns to features, in the code you can specify this using the dims keyword with dims=1 being the default whereas dims=2 will take rows as features and columns as observations.We will write X_c the centered matrix i.e. where each column sums up to zero. We will also write X_s the standardised matrix i.e. where each column not only sums up to zero but is scaled to have sample variance one. Finally, we will write S the standard sample covariance estimator (see below) of size ptimes p and D the diagonal matrix matching the diagonal of S.Note that the sample variance and covariance can be corrected or uncorrected (and this can be specified in the code). In order to avoid having to specify this everywhere in the document, it is useful to introduce a last symbol: kappa which is either set to n (uncorrected case) or (n-1) (corrected case).With these notations we can write:begineqnarray\n    S = kappa^-1X_c^T X_c labelsimple-covariance\n    X_s = X_c D^-1\nendeqnarray"
},

{
    "location": "man/methods/#Simple-estimator-1",
    "page": "Methods",
    "title": "Simple estimator",
    "category": "section",
    "text": "The standard covariance estimator is easily obtained via \\eqref{simple-covariance}. It can be specified with the constructor Simple which can take a named argument corrected (either false (default) or true).using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 5, 7\nX = randn(n, p)\n# corrected covariance\nS = cov(X, Simple(corrected=true))\n# we can also manually compute it and compare\nXc = (X .- sum(X, dims=1)/n) # centering\nκ = n-1\nS ≈ (Xc\'*Xc)/κ"
},

{
    "location": "man/methods/#Linear-shrinkage-estimators-1",
    "page": "Methods",
    "title": "Linear shrinkage estimators",
    "category": "section",
    "text": "Linear shrinkage estimators correspond to covariance estimators of the formC = (1-lambda)S + lambda Fwhere F is a target matrix of appropriate dimensions and lambdain01 is a shrinkage intensity. There are several standard targets that can be used, a simple example being the identity matrix.The shrinkage intensity lambda can be specified manually or computed automatically. Depending on the target, different approaches are implemented to compute a good intensity such as, for example, the Ledoit-Wolfe optimal intensity (which is the default intensity if you don\'t specify it).You can read more on the targets that can be used and the corresponding automatic intensities here.Here is an example using the identity matrix as a target and automatic shrinkage intensity (Ledoit-Wolfe):using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 2, 3\nX = randn(n, p)\ntarget = DiagonalUnitVariance()\nshrinkage = :lw # Ledoit-Wolfe optimal shrinkage\nmethod = LinearShrinkageEstimator(target, shrinkage)\ncov(X, method)You can also specify the intensity manually:using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1) # hide\nn, p = 2, 3 # hide\nX = randn(n, p) # hide\ntarget = DiagonalUnitVariance() # hide\nshrinkage = 0.8\nmethod2 = LinearShrinkageEstimator(target, shrinkage)\ncov(X, method2)"
},

{
    "location": "man/methods/#Nonlinear-shrinkage-estimators-1",
    "page": "Methods",
    "title": "Nonlinear shrinkage estimators",
    "category": "section",
    "text": "More on linear shrinkage estimators..."
},

{
    "location": "man/lshrink/#",
    "page": "Linear shrinkage estimators",
    "title": "Linear shrinkage estimators",
    "category": "page",
    "text": ""
},

{
    "location": "man/lshrink/#lshrink-1",
    "page": "Linear shrinkage estimators",
    "title": "Linear shrinkage estimators",
    "category": "section",
    "text": ""
},

{
    "location": "man/lshrink/#blabla-1",
    "page": "Linear shrinkage estimators",
    "title": "blabla",
    "category": "section",
    "text": ""
},

{
    "location": "man/nlshrink/#",
    "page": "Nonlinear shrinkage estimators",
    "title": "Nonlinear shrinkage estimators",
    "category": "page",
    "text": ""
},

{
    "location": "man/nlshrink/#nlshrink-1",
    "page": "Nonlinear shrinkage estimators",
    "title": "Nonlinear shrinkage estimators",
    "category": "section",
    "text": ""
},

{
    "location": "lib/public/#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public/#Public-Documentation-1",
    "page": "Public",
    "title": "Public Documentation",
    "category": "section",
    "text": "Documentation for CovarianceEstimation.jl\'s public interface.See Internal Documentation for internal package docs."
},

{
    "location": "lib/public/#Contents-1",
    "page": "Public",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public/#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public/#Statistics.cov",
    "page": "Public",
    "title": "Statistics.cov",
    "category": "function",
    "text": "cov(x::AbstractVector, sc::Simple)\n\nCompute the sample variance of the vector x. The sum is scaled with n where n = length(x) if sc.corrected is false and with n-1 otherwise.\n\n\n\n\n\ncov(x::AbstractVector, y, sc::Simple)\n\nCompute the covariance of the vectors x and y using formula frac1nsum_i=1^n (x_i-bar x) (y_i-bar y)^* where * denotes complex conjugate. If sc.corrected is true then the fraction frac1n is replaced with frac1n-1.\n\n\n\n\n\ncov(X::AbstractMatrix, sc::Simple; dims::Int=1)\n\nCompute the covariance matrix associated with X along dimension dims. The sum is scaled with n where  n = length(x) if sc.corrected is false and with n-1 otherwise.\n\n\n\n\n\ncov(X, sc::Simple{<:AbstractWeights}; dims=1)\n\nCompute the weighted covariance matrix.\n\n\n\n\n\ncov(X, lse::LinearShrinkageEstimator; dims=1)\n\nLinear shrinkage covariance estimator for matrix X along dimension dims. Computed using the method described by lse.\n\n\n\n\n\ncov(X, ans::AnalyticalNonlinearShrinkage; dims=1)\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.CovarianceEstimator",
    "page": "Public",
    "title": "CovarianceEstimation.CovarianceEstimator",
    "category": "type",
    "text": "CovarianceEstimator\n\nBasic type for all covariance estimators.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.Simple",
    "page": "Public",
    "title": "CovarianceEstimation.Simple",
    "category": "type",
    "text": "Simple(corrected::Bool, weights::Union{AbstractWeights, Nothing} = nothing)\n\nSimple covariance estimator (scaled by n-1 if corrected and n otherwise where n is the number of samples). Optionally supports weights (see http://juliastats.github.io/StatsBase.jl/stable/cov.html).\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.LinearShrinkageEstimator",
    "page": "Public",
    "title": "CovarianceEstimation.LinearShrinkageEstimator",
    "category": "type",
    "text": "LinearShrinkageEstimator(target, shrinkage; corrected = false)\n\nLinear shrinkage estimator described by equation (1 - lambda) S + lambda F where S is standard covariance matrix, F is shrinkage target described by argument target and lambda is a shrinkage parameter, either given explicitly in shrinkage or automatically determined according to one of the supported methods.\n\nThe corrected estimator is used if corrected is true.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.DiagonalUnitVariance",
    "page": "Public",
    "title": "CovarianceEstimation.DiagonalUnitVariance",
    "category": "type",
    "text": "DiagonalUnitVariance\n\nTarget for linear shrinkage: unit matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.DiagonalCommonVariance",
    "page": "Public",
    "title": "CovarianceEstimation.DiagonalCommonVariance",
    "category": "type",
    "text": "DiagonalCommonVariance\n\nTarget for linear shrinkage: unit matrix multiplied by average variance of variables.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.DiagonalUnequalVariance",
    "page": "Public",
    "title": "CovarianceEstimation.DiagonalUnequalVariance",
    "category": "type",
    "text": "DiagonalUnequalVariance\n\nTarget for linear shrinkage: diagonal of covariance matrix.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.CommonCovariance",
    "page": "Public",
    "title": "CovarianceEstimation.CommonCovariance",
    "category": "type",
    "text": "CommonCovariance\n\nTarget for linear shrinkage: see target_C.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.PerfectPositiveCorrelation",
    "page": "Public",
    "title": "CovarianceEstimation.PerfectPositiveCorrelation",
    "category": "type",
    "text": "PerfectPositiveCorrelation\n\nTarget for linear shrinkage: see target_E.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.ConstantCorrelation",
    "page": "Public",
    "title": "CovarianceEstimation.ConstantCorrelation",
    "category": "type",
    "text": "ConstantCorrelation\n\nTarget for linear shrinkage: see target_F.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#CovarianceEstimation.AnalyticalNonlinearShrinkage",
    "page": "Public",
    "title": "CovarianceEstimation.AnalyticalNonlinearShrinkage",
    "category": "type",
    "text": "AnalyticalNonlinearShrinkage\n\nAnalytical nonlinear shrinkage estimator. See docs for analytical_nonlinear_shrinkage for details.\n\n\n\n\n\n"
},

{
    "location": "lib/public/#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": "cov\r\nCovarianceEstimator\r\nSimple\r\nLinearShrinkageEstimator\r\nDiagonalUnitVariance\r\nDiagonalCommonVariance\r\nDiagonalUnequalVariance\r\nCommonCovariance\r\nPerfectPositiveCorrelation\r\nConstantCorrelation\r\nAnalyticalNonlinearShrinkage"
},

{
    "location": "lib/internals/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals/#Internal-Documentation-1",
    "page": "Internals",
    "title": "Internal Documentation",
    "category": "section",
    "text": "Documentation for CovarianceEstimation.jl\'s internal methods.See Public Documentation for public package docs."
},

{
    "location": "lib/internals/#Contents-1",
    "page": "Internals",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"internals.md\"]"
},

{
    "location": "lib/internals/#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"internals.md\"]"
},

{
    "location": "lib/internals/#CovarianceEstimation.rescale",
    "page": "Internals",
    "title": "CovarianceEstimation.rescale",
    "category": "function",
    "text": "rescale(M, d)\n\nInternal function to scale the rows and the columns of a square matrix M according to the elements in d. This is useful when dealing with the standardised data matrix which can be written Xs=Xc*D where Xc is the centered data matrix and D=Diagonal(d) so that its simple covariance is D*S*D where S is the simple covariance of Xc. Such D*M*D terms appear often in the computations of optimal shrinkage λ.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.uccov",
    "page": "Internals",
    "title": "CovarianceEstimation.uccov",
    "category": "function",
    "text": "uccov(X)\n\nInternal function to compute X*X\'/n where n is the number of rows of X. This corresponds to the uncorrected covariance of X if X is centered. This operation appears often in the computations of optimal shrinkage λ.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.sumij",
    "page": "Internals",
    "title": "CovarianceEstimation.sumij",
    "category": "function",
    "text": "sumij(S)\n\nInternal function to compute the sum of elements of a square matrix S. A keyword with_diag can be passed to indicate whether to include or not the diagonal of S in the sum. Both cases happen often in the computations of optimal shrinkage λ.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.square",
    "page": "Internals",
    "title": "CovarianceEstimation.square",
    "category": "function",
    "text": "square(x)\n\nInternal function to compute the square of a real number x. Defined here so it can be used in the other internal function sumij2.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.sumij2",
    "page": "Internals",
    "title": "CovarianceEstimation.sumij2",
    "category": "function",
    "text": "sumij2(S)\n\nInternal function identical to sumij except that it passes the function square to the sum so that it is the sum of the elements of S squared which is computed. This is significantly more efficient than using sumij(S.^2) for large matrices as it allocates very little.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.sum_fij",
    "page": "Internals",
    "title": "CovarianceEstimation.sum_fij",
    "category": "function",
    "text": "sum_fij(Xc, S, n, κ)\n\nInternal function corresponding to _ijf_ij that appears in http://strimmerlab.org/publications/journals/shrinkcov2005.pdf p.11.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.linear_shrinkage",
    "page": "Internals",
    "title": "CovarianceEstimation.linear_shrinkage",
    "category": "function",
    "text": "linear_shrinkage(target, Xc, S, λ, n, p, corrected)\n\nPerforms linear shrinkage with target of type target for data matrix Xc of size n by p with covariance matrix S and shrinkage parameter λ. Calculates corrected covariance if corrected is true.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.analytical_nonlinear_shrinkage",
    "page": "Internals",
    "title": "CovarianceEstimation.analytical_nonlinear_shrinkage",
    "category": "function",
    "text": "analytical_nonlinear_shrinkage(X)\n\nBased on Matlab code in Olivier Ledoit and Michael Wolf. Analytical Nonlinear Shrinkage of Large-Dimensional Covariance Matrices. (Nov 2018) http://www.econ.uzh.ch/static/wp/econwp264.pdf\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": "CovarianceEstimation.rescale\r\nCovarianceEstimation.uccov\r\nCovarianceEstimation.sumij\r\nCovarianceEstimation.square\r\nCovarianceEstimation.sumij2\r\nCovarianceEstimation.sum_fij\r\nCovarianceEstimation.linear_shrinkage\r\nCovarianceEstimation.analytical_nonlinear_shrinkage"
},

]}
