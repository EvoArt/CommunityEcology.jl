module CommunityEcologyMakieExt

import Makie
import CommunityEcology as CE
using CommunityEcology: Community, Ordination, PhyloTree, BetaDisper, NMDSResult
import CommunityEcology: plot_ordination, plot_ordination!, plot_composition, plot_composition!,
    plot_alpha, plot_alpha!, plot_rarefaction, plot_rarefaction!, plot_heatmap,
    plot_heatmap!, plot_dispersion, plot_dispersion!, plot_shepard, plot_shepard!,
    plot_tree, plot_tree!

const PALETTE = (:teal, :darkorange, :mediumpurple, :firebrick,
                 :goldenrod, :dodgerblue, :forestgreen, :sienna)
_color(i) = PALETTE[mod1(i, length(PALETTE))]

function _limits(values)
    lo, hi = extrema(values)
    lo == hi && return (lo - 0.5, hi + 0.5)
    return (lo, hi)
end

function _groups(values)
    labels = string.(values)
    levels = sort!(unique(labels))
    lookup = Dict(level => i for (i, level) in enumerate(levels))
    return levels, [lookup[label] for label in labels]
end

function _ordinal_colors(o, colorby)
    colorby === nothing && return :teal
    colorby in propertynames(o.metadata) ||
        throw(ArgumentError("metadata has no column $colorby"))
    values = o.metadata[!, colorby]
    if all(v -> !ismissing(v) && v isa Real, values)
        return Float64.(values)
    end
    _, ids = _groups(values)
    return [_color(i) for i in ids]
end

Makie.@recipe OrdinationPlot (ordination,) begin
    colorby = nothing
    labels = false
    markersize = 13
end

function Makie.plot!(plot::OrdinationPlot)
    o = plot[1][]
    size(o.points, 2) >= 2 || throw(ArgumentError("ordination needs two axes to plot"))
    Makie.scatter!(plot, o.points[:, 1], o.points[:, 2];
                   color=_ordinal_colors(o, plot.colorby[]),
                   colormap=:viridis, markersize=plot.markersize[])
    if plot.labels[]
        Makie.text!(plot, o.points[:, 1], o.points[:, 2]; text=o.samples,
                    offset=(7, 7), align=(:left, :bottom), fontsize=11)
    end
    return plot
end

plot_ordination!(axis, o::Ordination; kwargs...) = ordinationplot!(axis, o; kwargs...)

function plot_ordination(o::Ordination; colorby=nothing, labels::Bool=false,
                         size=(760, 500), kwargs...)
    length(o.axislabels) >= 2 || throw(ArgumentError("ordination needs two axes to plot"))
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel=o.axislabels[1], ylabel=o.axislabels[2],
                    title="$(uppercase(string(o.method))) · $(o.distance)")
    plot_ordination!(ax, o; colorby, labels, kwargs...)
    if colorby !== nothing
        values = o.metadata[!, colorby]
        if all(v -> !ismissing(v) && v isa Real, values)
            Makie.Colorbar(fig[1, 2]; colormap=:viridis, limits=_limits(Float64.(values)),
                           label=string(colorby))
        else
            levels, _ = _groups(values)
            elements = [Makie.MarkerElement(; color=_color(i), marker=:circle,
                                            markersize=13) for i in eachindex(levels)]
            Makie.Legend(fig[1, 2], elements, levels, string(colorby))
        end
    end
    return fig
end

function _composition(c, top, proportions)
    top >= 1 || throw(ArgumentError("top must be positive"))
    reduced = CE.lump_rare(c, top)
    X = proportions ? CE.relative(reduced).counts : Float64.(reduced.counts)
    return reduced, X
end

Makie.@recipe CompositionPlot (community,) begin
    top = 10
    proportions = true
end

function Makie.plot!(plot::CompositionPlot)
    c, X = _composition(plot[1][], plot.top[], plot.proportions[])
    rectangles = Makie.Rect2f[]
    colors = Symbol[]
    for i in 1:CE.nsamples(c)
        bottom = 0.0
        for j in 1:CE.nfeatures(c)
            height = X[i, j]
            if height > 0
                push!(rectangles, Makie.Rect2f(i - 0.43, bottom, 0.86, height))
                push!(colors, _color(j))
            end
            bottom += height
        end
    end
    isempty(rectangles) || Makie.poly!(plot, rectangles; color=colors, strokewidth=0)
    return plot
end

plot_composition!(axis, c::Community; kwargs...) = compositionplot!(axis, c; kwargs...)

function plot_composition(c::Community; top::Integer=10, proportions::Bool=true,
                          size=(850, 500), kwargs...)
    CE.nsamples(c) > 0 || throw(ArgumentError("community has no samples"))
    reduced, _ = _composition(c, top, proportions)
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel="Sample", ylabel=proportions ? "Relative abundance" : "Count",
                    xticks=(collect(1:CE.nsamples(c)), c.samples), xticklabelrotation=pi / 4)
    plot_composition!(ax, c; top, proportions, kwargs...)
    elements = [Makie.PolyElement(; color=_color(j)) for j in 1:CE.nfeatures(reduced)]
    Makie.Legend(fig[1, 2], elements, reduced.features, "Feature")
    return fig
end

function _alpha_data(c, measure, group)
    values = CE.alpha_diversity(c; measures=[measure], metadata=false)[!, measure]
    if group === nothing
        return Float64.(collect(1:CE.nsamples(c))), Float64.(values), c.samples
    end
    group in propertynames(c.metadata) || throw(ArgumentError("metadata has no column $group"))
    levels, ids = _groups(c.metadata[!, group])
    x = [ids[i] + 0.13 * sin(2.4 * i) for i in eachindex(ids)]
    return x, Float64.(values), levels
end

Makie.@recipe AlphaDiversityPlot (community,) begin
    measure = :shannon
    group = nothing
    markersize = 11
end

function Makie.plot!(plot::AlphaDiversityPlot)
    x, y, _ = _alpha_data(plot[1][], plot.measure[], plot.group[])
    Makie.scatter!(plot, x, y; color=:teal, markersize=plot.markersize[])
    if plot.group[] !== nothing
        ids = round.(Int, x)
        for g in unique(ids)
            vals = y[ids .== g]
            isempty(vals) || Makie.lines!(plot, [g - 0.28, g + 0.28],
                                          fill(sum(vals) / length(vals), 2);
                                          color=:darkorange, linewidth=3)
        end
    end
    return plot
end

plot_alpha!(axis, c::Community; kwargs...) = alphadiversityplot!(axis, c; kwargs...)

function plot_alpha(c::Community; measure::Symbol=:shannon, group=nothing,
                    size=(760, 470), kwargs...)
    _, _, labels = _alpha_data(c, measure, group)
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel=group === nothing ? "Sample" : string(group),
                    ylabel=string(measure), xticks=(collect(1:length(labels)), labels),
                    xticklabelrotation=pi / 4)
    plot_alpha!(ax, c; measure, group, kwargs...)
    return fig
end

Makie.@recipe RarefactionPlot (community,) begin
    step = nothing
    linewidth = 2
end

function Makie.plot!(plot::RarefactionPlot)
    c = plot[1][]
    curves = CE.rarecurve(c; step=plot.step[])
    for (i, sample) in enumerate(c.samples)
        rows = curves.sample .== sample
        Makie.lines!(plot, curves.depth[rows], curves.richness[rows];
                     color=_color(i), linewidth=plot.linewidth[])
    end
    return plot
end

plot_rarefaction!(axis, c::Community; kwargs...) = rarefactionplot!(axis, c; kwargs...)

function plot_rarefaction(c::Community; size=(760, 500), kwargs...)
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel="Reads sampled", ylabel="Expected richness")
    plot_rarefaction!(ax, c; kwargs...)
    if CE.nsamples(c) <= 12
        elements = [Makie.LineElement(; color=_color(i), linewidth=2) for i in 1:CE.nsamples(c)]
        Makie.Legend(fig[1, 2], elements, c.samples, "Sample")
    end
    return fig
end

Makie.@recipe CommunityHeatmapPlot (community,) begin
    proportions = true
    colormap = :viridis
end

function Makie.plot!(plot::CommunityHeatmapPlot)
    c = plot[1][]
    X = plot.proportions[] ? CE.relative(c).counts : Float64.(c.counts)
    Makie.heatmap!(plot, 1:CE.nsamples(c), 1:CE.nfeatures(c), X;
                   colormap=plot.colormap[])
    return plot
end

plot_heatmap!(axis, c::Community; kwargs...) = communityheatmapplot!(axis, c; kwargs...)

function plot_heatmap(c::Community; proportions::Bool=true, size=(850, 550), kwargs...)
    CE.nsamples(c) > 0 && CE.nfeatures(c) > 0 ||
        throw(ArgumentError("heatmap needs samples and features"))
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel="Sample", ylabel="Feature",
                    xticks=(collect(1:CE.nsamples(c)), c.samples),
                    yticks=(collect(1:CE.nfeatures(c)), c.features),
                    xticklabelrotation=pi / 4)
    plot_heatmap!(ax, c; proportions, kwargs...)
    X = proportions ? CE.relative(c).counts : c.counts
    Makie.Colorbar(fig[1, 2]; colormap=:viridis, limits=_limits(X),
                   label=proportions ? "Relative abundance" : "Count")
    return fig
end

Makie.@recipe DispersionPlot (dispersion,) begin
    markersize = 11
end

function Makie.plot!(plot::DispersionPlot)
    bd = plot[1][]
    levels, ids = _groups(bd.group)
    x = [ids[i] + 0.13 * sin(2.4 * i) for i in eachindex(ids)]
    Makie.scatter!(plot, x, bd.distances;
                   color=[_color(i) for i in ids], markersize=plot.markersize[])
    for (j, level) in enumerate(levels)
        values = bd.distances[string.(bd.group) .== level]
        Makie.lines!(plot, [j - 0.3, j + 0.3], fill(sum(values) / length(values), 2);
                     color=:black, linewidth=3)
    end
    return plot
end

plot_dispersion!(axis, bd::BetaDisper; kwargs...) = dispersionplot!(axis, bd; kwargs...)

function plot_dispersion(bd::BetaDisper; size=(700, 450), kwargs...)
    levels, _ = _groups(bd.group)
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel="Group", ylabel="Distance to group centre",
                    xticks=(collect(1:length(levels)), levels))
    plot_dispersion!(ax, bd; kwargs...)
    return fig
end

Makie.@recipe ShepardPlot (result,) begin
    markersize = 5
end

function Makie.plot!(plot::ShepardPlot)
    r = plot[1][]
    Makie.scatter!(plot, r.dissimilarities, r.distances;
                   color=(:teal, 0.35), markersize=plot.markersize[])
    order = sortperm(r.dissimilarities)
    Makie.lines!(plot, r.dissimilarities[order], r.fitted[order];
                 color=:darkorange, linewidth=2)
    return plot
end

plot_shepard!(axis, r::NMDSResult; kwargs...) = shepardplot!(axis, r; kwargs...)
plot_shepard!(axis, o::Ordination; kwargs...) = plot_shepard!(axis, _nmds_result(o); kwargs...)

function _nmds_result(o::Ordination)
    o.result isa NMDSResult || throw(ArgumentError("Shepard plot needs an NMDS ordination"))
    return o.result
end

function plot_shepard(r::NMDSResult; size=(700, 480), kwargs...)
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel="Observed dissimilarity",
                    ylabel="Ordination distance", title="Shepard diagram · stress $(round(r.stress; digits=3))")
    plot_shepard!(ax, r; kwargs...)
    return fig
end
plot_shepard(o::Ordination; kwargs...) = plot_shepard(_nmds_result(o); kwargs...)

function _tree_coordinates(t::PhyloTree)
    n = length(t.parent)
    x = zeros(n)
    y = zeros(n)
    for (i, tip) in enumerate(1:t.ntips)
        y[tip] = t.ntips - i + 1
    end
    for v in reverse(CE.postorder(t))
        for child in t.children[v]
            x[child] = x[v] + t.edge[child]
        end
    end
    for v in CE.postorder(t)
        isempty(t.children[v]) || (y[v] = sum(y[child] for child in t.children[v]) / length(t.children[v]))
    end
    return x, y
end

Makie.@recipe CommunityTreePlot (tree,) begin
    linecolor = :teal
    linewidth = 2
    labels = true
end

function Makie.plot!(plot::CommunityTreePlot)
    t = plot[1][]
    x, y = _tree_coordinates(t)
    xs = Float64[]
    ys = Float64[]
    for v in eachindex(t.parent)
        children = t.children[v]
        isempty(children) && continue
        append!(xs, (x[v], x[v], NaN))
        append!(ys, (minimum(y[children]), maximum(y[children]), NaN))
        for child in children
            append!(xs, (x[v], x[child], NaN))
            append!(ys, (y[child], y[child], NaN))
        end
    end
    Makie.lines!(plot, xs, ys; color=plot.linecolor[], linewidth=plot.linewidth[])
    if plot.labels[]
        Makie.text!(plot, x[1:t.ntips], y[1:t.ntips]; text=CE.tiplabels(t),
                    offset=(7, 0), align=(:left, :center), fontsize=11)
    end
    return plot
end

plot_tree!(axis, t::PhyloTree; kwargs...) = communitytreeplot!(axis, t; kwargs...)
plot_tree!(axis, c::Community; kwargs...) = plot_tree!(axis, _tree(c); kwargs...)

function _tree(c::Community)
    c.tree === nothing && throw(ArgumentError("community has no tree"))
    return c.tree
end

function plot_tree(t::PhyloTree; size=(800, 550), kwargs...)
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; xlabel="Branch length", yticks=([], []))
    plot_tree!(ax, t; kwargs...)
    Makie.hidespines!(ax, :l, :r, :t)
    return fig
end
plot_tree(c::Community; kwargs...) = plot_tree(_tree(c); kwargs...)

end
