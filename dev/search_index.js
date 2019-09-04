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
    "text": "Lightweight robust covariance estimation in Julia.A package for robustly estimating covariance matrices of real-valued data."
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
    "text": "Pages = [\"man/methods.md\", \"man/lshrink.md\", \"man/nlshrink.md\", \"man/msecomp.md\"]"
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
    "text": "In these docs, X denotes an ntimes p data matrix describing n observations with p features (variables) and possibly p  n. For most estimators, the matrix is assumed to have entries in mathbb R.note: Note\nWhile in the docs we assume that the rows of X correspond to observations and the columns to features, in the code you can specify this using the dims keyword with dims=1 being the default whereas dims=2 will take rows as features and columns as observations.We will write X_c the centered matrix i.e. where each column sums up to zero. We will also write X_s the standardised matrix i.e. where each column not only sums up to zero but is scaled to have sample variance one. Finally, we will write S the standard sample covariance estimator (see below) of size ptimes p and D the diagonal matrix matching the diagonal of S.Note that the sample variance and covariance can be corrected or uncorrected (and this can be specified in the code). In order to avoid having to specify this everywhere in the document, it is useful to introduce a last symbol: kappa which is either set to n (uncorrected case) or (n-1) (corrected case).With these notations we can write:begineqnarray\n    S = kappa^-1X_c^T X_c labelsimple-covariance\n    D_ij = S_ijmathbf 1_i=j\n    X_s = X_c D^-12\nendeqnarray"
},

{
    "location": "man/methods/#Simple-estimator-1",
    "page": "Methods",
    "title": "Simple estimator",
    "category": "section",
    "text": "The standard covariance estimator is easily obtained via \\eqref{simple-covariance}. It can be specified with the constructor SimpleCovariance which can take a named argument corrected (either false (default) or true).using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 5, 7\nX = randn(n, p)\n# corrected covariance\nS = cov(SimpleCovariance(corrected=true), X)\n# we can also manually compute it and compare\nXc = (X .- sum(X, dims=1)/n) # centering\nκ = n-1 # correction factor\nS ≈ (Xc\'*Xc)/κ"
},

{
    "location": "man/methods/#Linear-shrinkage-estimators-1",
    "page": "Methods",
    "title": "Linear shrinkage estimators",
    "category": "section",
    "text": "Linear shrinkage estimators correspond to covariance estimators of the formhatSigma = (1-lambda)S + lambda Fwhere F is a target matrix of appropriate dimensions, lambdain01 is a shrinkage intensity and S is the sample covariance estimator. There are several standard targets that can be used, a simple example being the identity matrix.The shrinkage intensity lambda can be specified manually or computed automatically. Depending on the target, different approaches are implemented to compute a good intensity such as, for example, the Ledoit-Wolf optimal intensity (which is the default intensity if you don\'t specify it).You can read more on the targets that can be used and the corresponding automatic intensities here.Here is an example using the identity matrix as a target and automatic shrinkage intensity (Ledoit-Wolfe):using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 2, 3\nX = randn(n, p)\ntarget = DiagonalUnitVariance()\nshrinkage = :lw # Ledoit-Wolf optimal shrinkage\nmethod = LinearShrinkage(target, shrinkage)\ncov(method, X)You can also specify the intensity manually:using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1) # hide\nn, p = 2, 3 # hide\nX = randn(n, p) # hide\ntarget = DiagonalUnitVariance() # hide\nshrinkage = 0.8\nmethod2 = LinearShrinkage(target, shrinkage)\ncov(method2, X)Read more on linear shrinkage estimators..."
},

{
    "location": "man/methods/#Nonlinear-shrinkage-estimators-1",
    "page": "Methods",
    "title": "Nonlinear shrinkage estimators",
    "category": "section",
    "text": "Read more on nonlinear shrinkage estimators..."
},

{
    "location": "man/methods/#Comparing-estimators-1",
    "page": "Methods",
    "title": "Comparing estimators",
    "category": "section",
    "text": "You may want to look at our simple comparison of covariance estimators which compares the MSE of the various estimators in a range of situations. Long story short, the LinearShrinkageEstimator with DiagonalUnequalVariance target performs well in the case np though most other estimators don\'t fare too badly in comparison. In the case np, the nonlinear shrinkage method does very well (though it is more expensive to compute)."
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
    "text": "Linear shrinkage estimators correspond to covariance estimators of the formhatSigma = (1-lambda)S + lambda Fwhere F is a target matrix of appropriate dimensions, lambdain01 is a shrinkage intensity and S is the sample covariance estimator (corrected or uncorrected depending on the corrected keyword)."
},

{
    "location": "man/lshrink/#Targets-and-intensities-1",
    "page": "Linear shrinkage estimators",
    "title": "Targets and intensities",
    "category": "section",
    "text": "There are several standard target matrices (denoted by F) that can be used (we follow here the notations and naming conventions of Schaffer & Strimmer 2005):Target name F_ii F_ij (ineq j) Comment\nDiagonalUnitVariance 1 0 F = mathbf I\nDiagonalCommonVariance v 0 F = vmathbf I\nDiagonalUnequalVariance S_ii 0 F = mathrmdiag(S), very common\nCommonCovariance v c \nPerfectPositiveCorrelation S_ii sqrtS_iiS_jj \nConstantCorrelation S_ii overlinersqrtS_iiS_jj used in Ledoit & Wolf 2004where $ v = \\mathrm{tr}(S)/p $ is the average variance, c = sum_ineq j S_ij(p(p-1)) is the average of off-diagonal terms of S and overliner is the average of sample correlations (see Schaffer & Strimmer 2005).For each of these targets, an optimal shrinkage intensity lambda^star can be computed. A standard approach is to apply the Ledoit-Wolf formula (shrinkage=:lw, see Ledoit & Wolf 2004) though there are some variants that can be applied too. Notably, Schaffer & Strimmer\'s variant (shrinkage=:ss) will ensure that the lambda^star computed is the same for X_c (the centered data matrix) as for X_s (the standardised data matrix). See Schaffer & Strimmer 2005.Chen\'s variant includes a Rao-Blackwellised estimator (shrinkage=:rblw) and an Oracle-Approximating one (shrinkage=:oas) for the DiagonalCommonVariance target. See Chen, Wiesel, Eldar & Hero 2010."
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
    "text": "Nonlinear shrinkage estimators correspond to covariance estimators based on the eigendecomposition of the sample matrix:F = eigen(X)\n# ... (transformation of eigenvalues)\nF.U*(d̃ .* F.U\') # d̃ is a vector of transformed eigenvaluesCurrently, only the analytical nonlinear shrinkage (AnalyticalNonlinearShrinkage) method is implemented."
},

{
    "location": "man/msecomp/#",
    "page": "MSE comparisons",
    "title": "MSE comparisons",
    "category": "page",
    "text": ""
},

{
    "location": "man/msecomp/#msecomp-1",
    "page": "MSE comparisons",
    "title": "MSE Comparison",
    "category": "section",
    "text": "Below are results obtained with a variety of data matrices of dimensions ntimes p. For each pair of dimension, 50 covariance matrices are generated with associated sample data matrices. The covariance obtained with the different estimators are then compared to the ground-truth and the MSE is reported.Abbreviation Method\nanshrink analytical nonlinear shrinkage\nccor LSE with constant correlation target\nccov LSE with constant covariance target\nd1v LSE with identity target\ndcv LSE with diagonal common variance target\nduv LSE with diagonal unequal variance target\nppc LSE with perfect positive correlation target\ns Simple estimator (baseline)\n_lw uses ledoit-wolf shrinkage\n_ss uses schaffer-strimmer shrinkage\n_oas uses oracle approximating shrinkage\n_rblw uses rao-blackwellised ledoit-wolf shrinkage"
},

{
    "location": "man/msecomp/#Fat-matrices-1",
    "page": "MSE comparisons",
    "title": "Fat matrices",
    "category": "section",
    "text": "(Image: ) (Image: ) (Image: ) (Image: )"
},

{
    "location": "man/msecomp/#Tall-matrices-1",
    "page": "MSE comparisons",
    "title": "Tall matrices",
    "category": "section",
    "text": "(Image: ) (Image: ) (Image: ) (Image: )"
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
    "location": "lib/public/#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": "cov\r\nStatsBase.CovarianceEstimator\r\nStatsBase.SimpleCovariance\r\nLinearShrinkage\r\nDiagonalUnitVariance\r\nDiagonalCommonVariance\r\nDiagonalUnequalVariance\r\nCommonCovariance\r\nPerfectPositiveCorrelation\r\nConstantCorrelation\r\nAnalyticalNonlinearShrinkage"
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
    "location": "lib/internals/#CovarianceEstimation.epanechnikov",
    "page": "Internals",
    "title": "CovarianceEstimation.epanechnikov",
    "category": "function",
    "text": "epanechnikov(x)\n\nReturn the Epanechnikov kernel evaluated at x.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.epanechnikov_HT1",
    "page": "Internals",
    "title": "CovarianceEstimation.epanechnikov_HT1",
    "category": "function",
    "text": "epnanechnikov_HT(x)\n\nReturn the Hilbert Transform of the Epanechnikov kernel evaluated at x if |x|≂̸√5.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.rescale",
    "page": "Internals",
    "title": "CovarianceEstimation.rescale",
    "category": "function",
    "text": "rescale(M, d)\n\nInternal function to scale the rows and the columns of a square matrix M according to the elements in d. This is useful when dealing with the standardised data matrix which can be written Xs=Xc*D where Xc is the centered data matrix and D=Diagonal(d) so that its simple covariance is D*S*D where S is the simple covariance of Xc. Such D*M*D terms appear often in the computations of optimal shrinkage λ.\n\nSpace complexity: O(p^2)\nTime complexity: O(p^2)\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.rescale!",
    "page": "Internals",
    "title": "CovarianceEstimation.rescale!",
    "category": "function",
    "text": "rescale!(M, d)\n\nSame as rescale but in place (no allocation).\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.uccov",
    "page": "Internals",
    "title": "CovarianceEstimation.uccov",
    "category": "function",
    "text": "uccov(X)\n\nInternal function to compute X*X\'/n where n is the number of rows of X. This corresponds to the uncorrected covariance of X if X is centered. This operation appears often in the computations of optimal shrinkage λ.\n\nSpace complexity: O(p^2)\nTime complexity: O(2np^2)\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.sumij",
    "page": "Internals",
    "title": "CovarianceEstimation.sumij",
    "category": "function",
    "text": "sumij(S)\n\nInternal function to compute the sum of elements of a square matrix S. A keyword with_diag can be passed to indicate whether to include or not the diagonal of S in the sum. Both cases happen often in the computations of optimal shrinkage λ.\n\nSpace complexity: O(1)\nTime complexity: O(p^2)\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.sumij2",
    "page": "Internals",
    "title": "CovarianceEstimation.sumij2",
    "category": "function",
    "text": "sumij2(S)\n\nInternal function identical to sumij except that it passes the function abs2 to the sum so that it is the sum of the elements of S squared which is computed. This is significantly more efficient than using sumij(S.^2) for large matrices as it allocates very little.\n\nSpace complexity: O(1)\nTime complexity: O(2p^2)\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#CovarianceEstimation.sum_fij",
    "page": "Internals",
    "title": "CovarianceEstimation.sum_fij",
    "category": "function",
    "text": "sum_fij(Xc, S, n, κ)\n\nInternal function corresponding to _ijf_ij that appears in http://strimmerlab.org/publications/journals/shrinkcov2005.pdf p.11.\n\nSpace complexity: O(np + 2p^2)\nTime complexity: O(2np^2)\n\n\n\n\n\n"
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
    "text": "analytical_nonlinear_shrinkage(S, n, p; decomp)\n\nInternal implementation of the analytical nonlinear shrinkage. The implementation is inspired from the Matlab code given in section C of Olivier Ledoit and Michael Wolf\'s paper \"Analytical Nonlinear Shrinkage of Large-Dimensional Covariance Matrices\". (Nov 2018) http://www.econ.uzh.ch/static/wp/econwp264.pdf\n\n\n\n\n\n"
},

{
    "location": "lib/internals/#Internals-1",
    "page": "Internals",
    "title": "Internals",
    "category": "section",
    "text": "CovarianceEstimation.epanechnikov\r\nCovarianceEstimation.epanechnikov_HT1\r\nCovarianceEstimation.rescale\r\nCovarianceEstimation.rescale!\r\nCovarianceEstimation.uccov\r\nCovarianceEstimation.sumij\r\nCovarianceEstimation.sumij2\r\nCovarianceEstimation.sum_fij\r\nCovarianceEstimation.linear_shrinkage\r\nCovarianceEstimation.analytical_nonlinear_shrinkage"
},

]}
