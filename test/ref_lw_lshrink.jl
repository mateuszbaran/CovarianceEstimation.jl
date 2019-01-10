# This is a Julia adaptation of the Matlab code from O. Ledoit and M. Wolf
# solely used to test our own implementation of their estimator.
# Their code was released under the BSD-2 license with the copyright notice
# below.
# https://www.econ.uzh.ch/en/people/faculty/wolf/publications.html#9

###########################################################################
# This file is released under the BSD 2-clause license.

# Copyright (c) 2014, Olivier Ledoit and Michael Wolf
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###########################################################################

using LinearAlgebra
using Statistics

# comments marked with a `#` come from their code
# comments marked with a `#<#>` are ours

function matlab_ledoitwolf_covcor(z)
    x = copy(z)

    # de-mean returns
    t, n  = size(x)
    meanx = mean(x, dims=1)
    x -= repeat(meanx, outer=(t, 1))

    # compute sample covariance matrix
    sample = (1/t) * (x'*x)

    # compute prior
    diagsample = diag(sample)
    sqrtvar = sqrt.(diagsample)

    sqrtvar_mat = repeat(sqrtvar, outer=(1, n))

    rBar = (sum(sample./(sqrtvar_mat .* sqrtvar_mat'))-n)/(n*(n-1))
    prior = rBar * sqrtvar_mat .* sqrtvar_mat'

    #<#> replace the diagonal
    prior -= Diagonal(diag(prior))
    prior += Diagonal(diagsample)

    #<#> compute optimal shrinkage

    # what we call pi-hat
    y = x.^2;
    phiMat = y'*y/t - 2*(x'*x) .* sample/t + sample.^2
    phi = sum(phiMat)

    # what we call rho-hat
    term1    = ((x.^3)'*x)/t
    help     = x'*x/t
    helpDiag = diag(help)
    term2    = repeat(helpDiag, outer=(1, n)) .* sample
    term3    = help .* repeat(diagsample, outer=(1, n))
    term4    = repeat(diagsample, outer=(1, n)) .* sample
    thetaMat = term1 - term2 - term3 + term4
    #<#> remove the diagonal
    thetaMat -= Diagonal(diag(thetaMat))
    rho       = sum(diag(phiMat)) + rBar*sum(sqrtvar'./sqrtvar.*thetaMat)

    # what we call gamma hat
    gamma = sum((sample - prior).^2)
    kappa = (phi - rho) / gamma
    shrinkage = max(0, min(1, kappa/t))

    part_results = Dict(
        "r̄" => rBar,
        "F" => prior,
        "shrinkage" => shrinkage,
        "kappa" => kappa,
        "gamma" => gamma,
        "lwcov" => shrinkage*prior + (1-shrinkage)*sample)

    return part_results
end
