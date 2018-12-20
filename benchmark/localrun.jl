using PkgBenchmark

results = benchmarkpkg("CovarianceEstimation")
export_markdown("benchmark/results.md", results)
