var documenterSearchIndex = {"docs":
[{"location":"man/lshrink/#lshrink","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"","category":"section"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"Linear shrinkage estimators correspond to covariance estimators of the form","category":"page"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"hatSigma = (1-lambda)S + lambda F","category":"page"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"where F is a target matrix of appropriate dimensions, lambdain01 is a shrinkage intensity and S is the sample covariance estimator (corrected or uncorrected depending on the corrected keyword).","category":"page"},{"location":"man/lshrink/#Targets-and-intensities","page":"Linear shrinkage estimators","title":"Targets and intensities","text":"","category":"section"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"There are several standard target matrices (denoted by F) that can be used (we follow here the notations and naming conventions of Schaffer & Strimmer 2005):","category":"page"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"Target name F_ii F_ij (ineq j) Comment\nDiagonalUnitVariance 1 0 F = mathbf I\nDiagonalCommonVariance v 0 F = vmathbf I\nDiagonalUnequalVariance S_ii 0 F = mathrmdiag(S), very common\nCommonCovariance v c \nPerfectPositiveCorrelation S_ii sqrtS_iiS_jj \nConstantCorrelation S_ii overlinersqrtS_iiS_jj used in Ledoit & Wolf 2004","category":"page"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"where $ v = \\mathrm{tr}(S)/p $ is the average variance, c = sum_ineq j S_ij(p(p-1)) is the average of off-diagonal terms of S and overliner is the average of sample correlations (see Schaffer & Strimmer 2005).","category":"page"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"For each of these targets, an optimal shrinkage intensity lambda^star can be computed. A standard approach is to apply the Ledoit-Wolf formula (shrinkage=:lw, see Ledoit & Wolf 2004) though there are some variants that can be applied too. Notably, Schaffer & Strimmer's variant (shrinkage=:ss) will ensure that the lambda^star computed is the same for X_c (the centered data matrix) as for X_s (the standardised data matrix). See Schaffer & Strimmer 2005.","category":"page"},{"location":"man/lshrink/","page":"Linear shrinkage estimators","title":"Linear shrinkage estimators","text":"Chen's variant includes a Rao-Blackwellised estimator (shrinkage=:rblw) and an Oracle-Approximating one (shrinkage=:oas) for the DiagonalCommonVariance target. See Chen, Wiesel, Eldar & Hero 2010.","category":"page"},{"location":"lib/internals/#Internal-Documentation","page":"Internals","title":"Internal Documentation","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"Documentation for CovarianceEstimation.jl's internal methods.","category":"page"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"See Public Documentation for public package docs.","category":"page"},{"location":"lib/internals/#Contents","page":"Internals","title":"Contents","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"lib/internals/#Index","page":"Internals","title":"Index","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"Pages = [\"internals.md\"]","category":"page"},{"location":"lib/internals/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"CovarianceEstimation.epanechnikov\r\nCovarianceEstimation.epanechnikov_HT1\r\nCovarianceEstimation.rescale\r\nCovarianceEstimation.rescale!\r\nCovarianceEstimation.uccov\r\nCovarianceEstimation.sumij\r\nCovarianceEstimation.sumij2\r\nCovarianceEstimation.sum_fij\r\nCovarianceEstimation.linear_shrinkage\r\nCovarianceEstimation.analytical_nonlinear_shrinkage","category":"page"},{"location":"lib/internals/#CovarianceEstimation.epanechnikov","page":"Internals","title":"CovarianceEstimation.epanechnikov","text":"epanechnikov(x)\n\nReturn the Epanechnikov kernel evaluated at x.\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.epanechnikov_HT1","page":"Internals","title":"CovarianceEstimation.epanechnikov_HT1","text":"epnanechnikov_HT(x)\n\nReturn the Hilbert Transform of the Epanechnikov kernel evaluated at x if |x|≂̸√5.\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.rescale","page":"Internals","title":"CovarianceEstimation.rescale","text":"rescale(M, d)\n\nInternal function to scale the rows and the columns of a square matrix M according to the elements in d. This is useful when dealing with the standardised data matrix which can be written Xs=Xc*D where Xc is the centered data matrix and D=Diagonal(d) so that its simple covariance is D*S*D where S is the simple covariance of Xc. Such D*M*D terms appear often in the computations of optimal shrinkage λ.\n\nSpace complexity: O(p^2)\nTime complexity: O(p^2)\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.rescale!","page":"Internals","title":"CovarianceEstimation.rescale!","text":"rescale!(M, d)\n\nSame as rescale but in place (no allocation).\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.uccov","page":"Internals","title":"CovarianceEstimation.uccov","text":"uccov(X)\n\nInternal function to compute X*X'/n where n is the number of rows of X. This corresponds to the uncorrected covariance of X if X is centered. This operation appears often in the computations of optimal shrinkage λ.\n\nSpace complexity: O(p^2)\nTime complexity: O(2np^2)\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.sumij","page":"Internals","title":"CovarianceEstimation.sumij","text":"sumij(S)\n\nInternal function to compute the sum of elements of a square matrix S. A keyword with_diag can be passed to indicate whether to include or not the diagonal of S in the sum. Both cases happen often in the computations of optimal shrinkage λ.\n\nSpace complexity: O(1)\nTime complexity: O(p^2)\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.sumij2","page":"Internals","title":"CovarianceEstimation.sumij2","text":"sumij2(S)\n\nInternal function identical to sumij except that it passes the function abs2 to the sum so that it is the sum of the elements of S squared which is computed. This is significantly more efficient than using sumij(S.^2) for large matrices as it allocates very little.\n\nSpace complexity: O(1)\nTime complexity: O(2p^2)\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.sum_fij","page":"Internals","title":"CovarianceEstimation.sum_fij","text":"sum_fij(Xc, S, n, κ)\n\nInternal function corresponding to _ijf_ij that appears in https://strimmerlab.github.io/publications/journals/shrinkcov2005.pdf p.11.\n\nSpace complexity: O(np + 2p^2)\nTime complexity: O(2np^2)\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.linear_shrinkage","page":"Internals","title":"CovarianceEstimation.linear_shrinkage","text":"linear_shrinkage(target, Xc, S, λ, n, p, corrected)\n\nPerforms linear shrinkage with target of type target for data matrix Xc of size n by p with covariance matrix S and shrinkage parameter λ. Calculates corrected covariance if corrected is true.\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/#CovarianceEstimation.analytical_nonlinear_shrinkage","page":"Internals","title":"CovarianceEstimation.analytical_nonlinear_shrinkage","text":"analytical_nonlinear_shrinkage(S, n, p; decomp)\n\nInternal implementation of the analytical nonlinear shrinkage. The implementation is inspired from the Matlab code given in section C of Olivier Ledoit and Michael Wolf's paper \"Analytical Nonlinear Shrinkage of Large-Dimensional Covariance Matrices\". (Nov 2018) http://www.econ.uzh.ch/static/wp/econwp264.pdf\n\n\n\n\n\n","category":"function"},{"location":"man/methods/#Methods","page":"Methods","title":"Methods","text":"","category":"section"},{"location":"man/methods/#Notations","page":"Methods","title":"Notations","text":"","category":"section"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"In these docs, X denotes an ntimes p data matrix describing n observations with p features (variables) and possibly p  n. For most estimators, the matrix is assumed to have entries in mathbb R.","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"note: Note\nWhile in the docs we assume that the rows of X correspond to observations and the columns to features, in the code you can specify this using the dims keyword with dims=1 being the default whereas dims=2 will take rows as features and columns as observations.","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"We will write X_c the centered matrix i.e. where each column sums up to zero. We will also write X_s the standardised matrix i.e. where each column not only sums up to zero but is scaled to have sample variance one. Finally, we will write S the standard sample covariance estimator (see below) of size ptimes p and D the diagonal matrix matching the diagonal of S.","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"Note that the sample variance and covariance can be corrected or uncorrected (and this can be specified in the code). In order to avoid having to specify this everywhere in the document, it is useful to introduce a last symbol: kappa which is either set to n (uncorrected case) or (n-1) (corrected case).","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"With these notations we can write:","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"beginaligned\n    S = kappa^-1X_c^T X_c\n    D_ij = S_ijmathbf 1_i=j\n    X_s = X_c D^-12\nendaligned","category":"page"},{"location":"man/methods/#Simple-estimator","page":"Methods","title":"Simple estimator","text":"","category":"section"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"The standard covariance estimator is easily obtained via the equation above. It can be specified with the constructor SimpleCovariance which can take a named argument corrected (either false (default) or true).","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 5, 7\nX = randn(n, p)\n# corrected covariance\nS = cov(SimpleCovariance(corrected=true), X)\n# we can also manually compute it and compare\nXc = (X .- sum(X, dims=1)/n) # centering\nκ = n-1 # correction factor\nS ≈ (Xc'*Xc)/κ","category":"page"},{"location":"man/methods/#Linear-shrinkage-estimators","page":"Methods","title":"Linear shrinkage estimators","text":"","category":"section"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"Linear shrinkage estimators correspond to covariance estimators of the form","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"hatSigma = (1-lambda)S + lambda F","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"where F is a target matrix of appropriate dimensions, lambdain01 is a shrinkage intensity and S is the sample covariance estimator. There are several standard targets that can be used, a simple example being the identity matrix.","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"The shrinkage intensity lambda can be specified manually or computed automatically. Depending on the target, different approaches are implemented to compute a good intensity such as, for example, the Ledoit-Wolf optimal intensity (which is the default intensity if you don't specify it).","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"You can read more on the targets that can be used and the corresponding automatic intensities here.","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"Here is an example using the identity matrix as a target and automatic shrinkage intensity (Ledoit-Wolfe):","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 2, 3\nX = randn(n, p)\ntarget = DiagonalUnitVariance()\nshrinkage = :lw # Ledoit-Wolf optimal shrinkage\nmethod = LinearShrinkage(target, shrinkage)\ncov(method, X)","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"You can also specify the intensity manually:","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"using CovarianceEstimation # hide\nusing Random # hide\nRandom.seed!(1) # hide\nn, p = 2, 3 # hide\nX = randn(n, p) # hide\ntarget = DiagonalUnitVariance() # hide\nshrinkage = 0.8\nmethod2 = LinearShrinkage(target, shrinkage)\ncov(method2, X)","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"Read more on linear shrinkage estimators...","category":"page"},{"location":"man/methods/#Nonlinear-shrinkage-estimators","page":"Methods","title":"Nonlinear shrinkage estimators","text":"","category":"section"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"Read more on nonlinear shrinkage estimators...","category":"page"},{"location":"man/methods/#Biweight-midcovariance","page":"Methods","title":"Biweight midcovariance","text":"","category":"section"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"The biweight midcovariance is a covariance estimator that is resilient to outliers. A full description of the technique is included on BiweightMidcovariance.","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"using CovarianceEstimation # hide\nusing Distributions # hide\nusing Random # hide\nRandom.seed!(1)\nn, p = 10, 3\nX = rand(TDist(3), (n, p))  # Moderately heavy-tailed data\ncov(BiweightMidcovariance(), X)","category":"page"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"The two controllable parameters are the keyword arguments c and modify_sample_size, whose purpose is described in the docstring.","category":"page"},{"location":"man/methods/#Comparing-estimators","page":"Methods","title":"Comparing estimators","text":"","category":"section"},{"location":"man/methods/","page":"Methods","title":"Methods","text":"You may want to look at our simple comparison of covariance estimators which compares the MSE of the various estimators in a range of situations. Long story short, the LinearShrinkageEstimator with DiagonalUnequalVariance target performs well in the case np though most other estimators don't fare too badly in comparison. In the case np, the nonlinear shrinkage method does very well (though it is more expensive to compute).","category":"page"},{"location":"lib/public/#Public-Documentation","page":"Public","title":"Public Documentation","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"Documentation for CovarianceEstimation.jl's public interface.","category":"page"},{"location":"lib/public/","page":"Public","title":"Public","text":"See Internal Documentation for internal package docs.","category":"page"},{"location":"lib/public/#Contents","page":"Public","title":"Contents","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#Index","page":"Public","title":"Index","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#Public-Interface","page":"Public","title":"Public Interface","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"cov\r\nStatsBase.CovarianceEstimator\r\nStatsBase.SimpleCovariance\r\nLinearShrinkage\r\nDiagonalUnitVariance\r\nDiagonalCommonVariance\r\nDiagonalUnequalVariance\r\nCommonCovariance\r\nPerfectPositiveCorrelation\r\nConstantCorrelation\r\nAnalyticalNonlinearShrinkage\r\nBiweightMidcovariance","category":"page"},{"location":"lib/public/#Statistics.cov","page":"Public","title":"Statistics.cov","text":"cov(lse::LinearShrinkage, X; dims=1)\n\nLinear shrinkage covariance estimator for matrix X along dimension dims. Computed using the method described by lse.\n\n\n\n\n\ncov(ans::AnalyticalNonlinearShrinkage, X; dims=1, mean=nothing)\n\nNonlinear covariance estimator derived from the sample covariance estimator S and its eigenvalue decomposition (which can be given through decomp). See Ledoit and Wolf's paper http://www.econ.uzh.ch/static/wp/econwp264.pdf The keyword mean can be nothing (centering via estimated mean), zero (no centering) or a provided vector. In the first case, a rank-1 modification is applied and therefore the effective sample size is decreased by one (see analytical_nonlinear_shrinkage). In the latter two case the mean cannot have been estimated on the data (otherwise the effective sample size will be 1 larger than it should be resulting in numerical instabilities). If you are unsure, use either nothing or provide an explicit (non-estimated) vector (possibly a zero vector) and avoid the use of mean=0.\n\nTime complexity (including formation of S)\n(p<n): O(np^2 + n^2) with moderate constant\n(p>n): O(p^3) with low constant (dominated by eigendecomposition of S)\n\n\n\n\n\n","category":"function"},{"location":"lib/public/#CovarianceEstimation.LinearShrinkage","page":"Public","title":"CovarianceEstimation.LinearShrinkage","text":"LinearShrinkage(target, shrinkage; corrected=false)\n\nLinear shrinkage estimator described by equation (1 - lambda) S + lambda F where S is standard covariance matrix, F is shrinkage target described by argument target and lambda is a shrinkage parameter, either given explicitly in shrinkage or automatically determined according to one of the supported methods.\n\nThe corrected estimator is used if corrected is true.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.DiagonalUnitVariance","page":"Public","title":"CovarianceEstimation.DiagonalUnitVariance","text":"DiagonalUnitVariance\n\nTarget for linear shrinkage: unit matrix. A subtype of LinearShrinkageTarget where\n\nF_ij=1 if i=j and\nF_ij=0 otherwise\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.DiagonalCommonVariance","page":"Public","title":"CovarianceEstimation.DiagonalCommonVariance","text":"DiagonalCommonVariance\n\nTarget for linear shrinkage: unit matrix multiplied by average variance of variables. A subtype of LinearShrinkageTarget where\n\nF_ij=v if i=j with v=mathrmtr(S)p and\nF_ij=0 otherwise\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.DiagonalUnequalVariance","page":"Public","title":"CovarianceEstimation.DiagonalUnequalVariance","text":"DiagonalUnequalVariance\n\nTarget for linear shrinkage: diagonal of covariance matrix. A subtype of LinearShrinkageTarget where\n\nF_ij=s_ij if i=j and\nF_ij=0 otherwise\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.CommonCovariance","page":"Public","title":"CovarianceEstimation.CommonCovariance","text":"CommonCovariance\n\nTarget for linear shrinkage: see target_C. A subtype of LinearShrinkageTarget where\n\nF_ij=v if i=j with v=mathrmtr(S)p and\nF_ij=c with c=sum_ineq j S_ij(p(p-1)) otherwise\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.PerfectPositiveCorrelation","page":"Public","title":"CovarianceEstimation.PerfectPositiveCorrelation","text":"PerfectPositiveCorrelation\n\nTarget for linear shrinkage: see target_E. A subtype of LinearShrinkageTarget where\n\nF_ij=S_ij if i=j and\nF_ij=sqrtS_iiS_jj otherwise\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.ConstantCorrelation","page":"Public","title":"CovarianceEstimation.ConstantCorrelation","text":"ConstantCorrelation\n\nTarget for linear shrinkage: see target_F. A subtype of LinearShrinkageTarget where\n\nF_ij=S_ij if i=j and\nF_ij=overlinersqrtS_iiS_jj otherwise where overliner is the average sample correlation\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.AnalyticalNonlinearShrinkage","page":"Public","title":"CovarianceEstimation.AnalyticalNonlinearShrinkage","text":"AnalyticalNonlinearShrinkage\n\nAnalytical nonlinear shrinkage estimator. See docs for analytical_nonlinear_shrinkage for details.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CovarianceEstimation.BiweightMidcovariance","page":"Public","title":"CovarianceEstimation.BiweightMidcovariance","text":"BiweightMidcovariance(; c=9.0, modify_sample_size=false)\n\nThe biweight midcovariance is a covariance estimator that is resilient to outliers.\n\nThe technique derives originally from astrophysics [1], and is implemented in the Python module Astropy [2], as well as in NIST's Dataplot [3].\n\nConsider two random variables x and y, for which we have n observations (x_i y_i). The biweight midcovariance is then defined to be:\n\nn_scdotfrac\n    sum_u_i1v_i1(x_i - M_x)(1 - u_i^2)^2(y_i - M_y)(1 - v_i^2)^2\n\n    left(sum_u_i1(1 - u_i^2)(1-5u_i^2)right)\n    left(sum_v_i1(1 - v_i^2)(1-5v_i^2)right)\n\n\nwhere n_s is the sample size, M_x and M_y are the medians of x_i and y_i respectively, and\n\nbeginaligned\nu_i = fracx_i - M_xc cdot mathrmMAD_x \nv_i = fracy_i - M_yc cdot mathrmMAD_y\nendaligned\n\nwhere mathrmMAD represents the median absolute deviation,\n\nbeginaligned\nmathrmMAD_x = mathrmmedian(leftleftx_i - M_xrightright) \nmathrmMAD_y = mathrmmedian(leftlefty_i - M_yrightright)\nendaligned\n\nIf either mathrmMAD_x = 0 or mathrmMAD_y = 0, the pairwise covariance is defined to be zero.\n\nThe parameter c is a tuning constant, for which the default is 90. Larger values will reduce the number of outliers that are removed — i.e. reducing robustness, but increasing sample efficiency.\n\nFields\n\nc::Float64: The tuning constant corresponding to c above.\nmodify_sample_size::Bool: If false, then we use a sample size n_s equal to the   total number of observations n. This is consistent with the standard definition of   biweight midcovariance in the literature. Otherwise, we count only those elements which   are not rejected as outliers in the numerator, i.e. those for which u_i1   and v_i1.   This follows the implementation in astropy [2].\n\nComplexity\n\nSpace: O(p^2)\nTime: O(np^2)\n\nReferences\n\n[1] Beers, Flynn, and Gebhardt (1990; AJ 100, 32) \"Measures of Location and Scale for Velocities in Clusters of Galaxies – A Robust Approach\"\n\n[2] Astropy biweight_midcovariance\n\n[3] NIST Dataplot biweight midcovariance\n\n\n\n\n\n","category":"type"},{"location":"#CovarianceEstimation.jl","page":"Home","title":"CovarianceEstimation.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Lightweight robust covariance estimation in Julia.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package for robustly estimating covariance matrices of real-valued data.","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Standard corrected and uncorrected covariance estimators,\nLinear and Nonlinear shrinkage estimators\nFocus on speed and lightweight dependencies","category":"page"},{"location":"#Manual-outline","page":"Home","title":"Manual outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"man/methods.md\", \"man/lshrink.md\", \"man/nlshrink.md\", \"man/msecomp.md\"]","category":"page"},{"location":"#Library-Outline","page":"Home","title":"Library Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"lib/public.md\", \"lib/internals.md\"]","category":"page"},{"location":"#main-index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"lib/public.md\"]","category":"page"},{"location":"man/nlshrink/#nlshrink","page":"Nonlinear shrinkage estimators","title":"Nonlinear shrinkage estimators","text":"","category":"section"},{"location":"man/nlshrink/","page":"Nonlinear shrinkage estimators","title":"Nonlinear shrinkage estimators","text":"Nonlinear shrinkage estimators correspond to covariance estimators based on the eigendecomposition of the sample matrix:","category":"page"},{"location":"man/nlshrink/","page":"Nonlinear shrinkage estimators","title":"Nonlinear shrinkage estimators","text":"F = eigen(X)\n# ... (transformation of eigenvalues)\nF.U*(d̃ .* F.U') # d̃ is a vector of transformed eigenvalues","category":"page"},{"location":"man/nlshrink/","page":"Nonlinear shrinkage estimators","title":"Nonlinear shrinkage estimators","text":"Currently, only the analytical nonlinear shrinkage (AnalyticalNonlinearShrinkage) method is implemented.","category":"page"},{"location":"man/msecomp/#msecomp","page":"MSE comparisons","title":"MSE Comparison","text":"","category":"section"},{"location":"man/msecomp/","page":"MSE comparisons","title":"MSE comparisons","text":"Below are results obtained with a variety of data matrices of dimensions ntimes p. For each pair of dimension, 50 covariance matrices are generated with associated sample data matrices. The covariance obtained with the different estimators are then compared to the ground-truth and the MSE is reported.","category":"page"},{"location":"man/msecomp/","page":"MSE comparisons","title":"MSE comparisons","text":"Abbreviation Method\nanshrink analytical nonlinear shrinkage\nccor LSE with constant correlation target\nccov LSE with constant covariance target\nd1v LSE with identity target\ndcv LSE with diagonal common variance target\nduv LSE with diagonal unequal variance target\nppc LSE with perfect positive correlation target\ns Simple estimator (baseline)\n_lw uses ledoit-wolf shrinkage\n_ss uses schaffer-strimmer shrinkage\n_oas uses oracle approximating shrinkage\n_rblw uses rao-blackwellised ledoit-wolf shrinkage","category":"page"},{"location":"man/msecomp/#Fat-matrices","page":"MSE comparisons","title":"Fat matrices","text":"","category":"section"},{"location":"man/msecomp/","page":"MSE comparisons","title":"MSE comparisons","text":"(Image: ) (Image: ) (Image: ) (Image: )","category":"page"},{"location":"man/msecomp/#Tall-matrices","page":"MSE comparisons","title":"Tall matrices","text":"","category":"section"},{"location":"man/msecomp/","page":"MSE comparisons","title":"MSE comparisons","text":"(Image: ) (Image: ) (Image: ) (Image: )","category":"page"}]
}
