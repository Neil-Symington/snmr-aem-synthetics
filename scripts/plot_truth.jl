#!/usr/bin/env julia
# Plot a built column: lithology with unit properties, clay, free water and
# conductivity on a common depth axis.
#
#   julia --project=. scripts/plot_truth.jl output/confined_aquifer_truth.jld2
using JLD2, Statistics, Printf, PyPlot
const mpl = PyPlot.matplotlib
length(ARGS) >= 1 || error("usage: plot_truth.jl <truth.jld2> [zmax]")
T = load(ARGS[1]); ZB = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 110.0
zf = T["z"]; units = T["units"]; tr = T["trans"]
out = replace(ARGS[1], "_truth.jld2" => "_truth.png")

# colour by lithology from the unit name, not by position in the list
function litho(nm)
    n = lowercase(nm)
    occursin("unweathered", n) || occursin("basement", n) ? ("#C9A0A8", "++") :
    occursin("saprock", n)                                ? ("#B0A0A8", "xxx") :
    occursin("weathered", n)                              ? ("#9AA86B", "xx") :
    occursin("clay", n) || occursin("aquitard", n)        ? ("#6E5B47", "---") :
    occursin("aquifer", n)                                ? ("#FFC94D", "..") :
    occursin("unsaturated", n)                            ? ("#E3D6B8", "..") :
    occursin("silt", n)                                   ? ("#B8A98C", "-.") :
                                                            ("#D8C9A8", "//")
end
um(x, a, b) = begin
    k = findall(zi -> a + 3 <= zi <= b - 3, zf); isempty(k) && (k = findall(zi -> a <= zi <= b, zf))
    isempty(k) ? NaN : mean(x[k])
end

fig = figure(figsize = (13.5, 8.0)); gs = fig.add_gridspec(1, 4, width_ratios = [2.4, 1, 1, 1], wspace = 0.28)
ax = fig.add_subplot(get(gs, (0, 0)))
for u in units
    t, b = u[1], min(u[2], ZB); t >= ZB && continue
    fc, hh = litho(u[3])
    ax.add_patch(mpl.patches.Rectangle((0, t), 1, b - t, facecolor = fc, edgecolor = "k", lw = .8, hatch = hh))
end
for (i, t) in enumerate(tr)
    b = units[i+1][1]; b - t/2 >= ZB && continue
    ax.add_patch(mpl.patches.Rectangle((0, b - t/2), 1, t, facecolor = "white", alpha = .45, edgecolor = "none", zorder = 2))
end
for u in units
    t, b = u[1], min(u[2], ZB); t >= ZB && continue; b - t < 6 && continue
    ax.text(0.5, (t + b)/2, @sprintf("%s\nφ %.2f  clay %.2f\nθ %.3f   σ_w %.2f",
            u[3], um(T["φ"], u[1], u[2]), um(T["clay"], u[1], u[2]),
            um(T["θ_SMR"], u[1], u[2]), um(T["σw"], u[1], u[2])),
        ha = "center", va = "center", fontsize = 7.6, zorder = 6, linespacing = 1.3,
        bbox = Dict("fc" => "white", "alpha" => 0.88, "ec" => "0.5", "lw" => 0.4, "pad" => 2.2))
end
ax.set_xlim(0, 1); ax.set_ylim(ZB, 0); ax.set_xticks([]); ax.set_ylabel("depth (m)")
ax.set_title("lithology and unit properties", fontsize = 11)
kk = findall(zi -> zi <= ZB, zf)
for (j, (x, lab, xl, c)) in enumerate(((T["clay"], "clay fraction", (0, 1.0), "tab:brown"),
                                       (T["θ_SMR"], "free water θ (v/v)", (0, 0.42), "tab:green")))
    a = fig.add_subplot(get(gs, (0, j)), sharey = ax)
    a.plot(x[kk], zf[kk], lw = 1.7, color = c); a.set_xlim(xl)
    a.set_title(lab, fontsize = 11); a.grid(alpha = .25); a.tick_params(labelleft = false)
    lab == "free water θ (v/v)" && a.axvline(0.15, color = "grey", ls = ":", lw = 1.0)
end
a = fig.add_subplot(get(gs, (0, 3)), sharey = ax)
a.semilogx(T["σ_bulk"][kk], zf[kk], color = "crimson", lw = 2.0, label = "σ_bulk (AEM sees)")
a.semilogx(T["σw"][kk], zf[kk], color = "tab:purple", lw = 1.5, ls = "--", label = "σ_w (pore water)")
a.set_xlim(5e-3, 3); a.legend(fontsize = 7.6, loc = "lower left"); a.grid(alpha = .25, which = "both")
a.set_title("conductivity (S/m)", fontsize = 11); a.tick_params(labelleft = false)
fig.suptitle(get(get(T, "config", Dict()), "name", "synthetic column"), fontsize = 12.5)
fig.subplots_adjust(top = .93, bottom = .06, left = .055, right = .985)
fig.savefig(out, dpi = 165); close(fig)
println("wrote ", out)
