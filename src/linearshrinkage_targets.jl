# Taxonomy from http://strimmerlab.org/publications/journals/shrinkcov2005.pdf
# page 13

struct DiagonalUnitVariance <: LinearShrinkageTarget end

struct DiagonalCommonVariance <: LinearShrinkageTarget end

struct CommonCovariance <: LinearShrinkageTarget end

struct DiagonalUnequalVariance <: LinearShrinkageTarget end

struct PerfectPositiveCorrelation <: LinearShrinkageTarget end

struct ConstantCorrelation <: LinearShrinkageTarget end
