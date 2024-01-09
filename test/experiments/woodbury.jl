using CairoMakie
using CairoMakie.Colors

# Attempt to match the colors in Fig 3. of Donoho et al. (2018)
const colors = [colorant"green1", colorant"cyan", colorant"steelblue1", colorant"blue", colorant"purple3", colorant"magenta", colorant"violetred"]

shrink(loss, λ) = CovarianceEstimation.shrink(loss, λ, 1, 1) + 1   # +1 because that's incorporated into the Woodbury σ² term
function ℓ(λ, γ=1)
    λ′ = λ + 1 - γ
    return (λ′ + sqrt(λ′^2 - 4λ)) / 2
end

λ = range(4, 9, length=101)

fig = Figure(size = (800, 800))#, font = "CMU Serif")

# L2 norm
ax = Axis(fig[1, 1], xlabel="λ", ylabel="η(λ)", title="norm=:L2, γ=1")
for i = 1:7
    i ∈ (5, 7) && continue
    loss = NormLossCov(:L2, i)
    lines!(ax, λ, shrink.(Ref(loss), λ); color=colors[i], label="F+$i")
end
lines!(ax, λ, λ, color=:red, linestyle=:dash, label="η(λ) = λ")
lines!(ax, λ, ℓ.(λ), color=:black, linestyle=:dash, label="ℓ(λ)")
ylims!(ax, (1, 9))
axislegend(ax, position=:lt)

# L1 norm
ax = Axis(fig[2, 1], xlabel="λ", ylabel="η(λ)", title="norm=:L1, γ=1")
for i = 1:7
    i ∈ (5, 7) && continue
    loss = NormLossCov(:L1, i)
    lines!(ax, λ, shrink.(Ref(loss), λ); color=colors[i], label="N+$i")
end
lines!(ax, λ, λ, color=:red, linestyle=:dash, label="η(λ) = λ")
lines!(ax, λ, ℓ.(λ), color=:black, linestyle=:dash, label="ℓ(λ)")
ylims!(ax, (1, 9))
axislegend(ax, position=:lt)

# Linf norm
ax = Axis(fig[1, 2], xlabel="λ", ylabel="η(λ)", title="norm=:Linf, γ=1")
for i = 1:7
    i ∈ (3, 4, 5, 7) && continue
    loss = NormLossCov(:Linf, i)
    lines!(ax, λ, shrink.(Ref(loss), λ); color=colors[i], label="O+$i")
end
lines!(ax, λ, λ, color=:red, linestyle=:dash, label="η(λ) = λ")
lines!(ax, λ, ℓ.(λ), color=:black, linestyle=:dash, label="ℓ(λ)")
ylims!(ax, (1, 9))
axislegend(ax, position=:lt)

# Statistical losses
ax = Axis(fig[2, 2], xlabel="λ", ylabel="η(λ)", title="Statistical losses (γ=1)")
lines!(ax, λ, shrink.(Ref(StatLossCov(:st)), λ); color=colors[1], label="Stein")
lines!(ax, λ, shrink.(Ref(StatLossCov(:ent)), λ); color=colors[2], label="Entropy")
lines!(ax, λ, shrink.(Ref(StatLossCov(:div)), λ); color=colors[3], label="Divergence")
lines!(ax, λ, shrink.(Ref(StatLossCov(:aff)), λ); color=colors[6], label="Affine")
lines!(ax, λ, shrink.(Ref(StatLossCov(:fre)), λ); color=colors[4], label="Frechet")
lines!(ax, λ, λ, color=:red, linestyle=:dash, label="η(λ) = λ")
lines!(ax, λ, ℓ.(λ), color=:black, linestyle=:dash, label="ℓ(λ)")
ylims!(ax, (1, 9))
axislegend(ax, position=:lt)

save(joinpath(pkgdir(CovarianceEstimation), "docs", "src", "assets", "donoho_fig3.png"), fig)
