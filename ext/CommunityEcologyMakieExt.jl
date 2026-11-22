module CommunityEcologyMakieExt

import Makie
import CommunityEcology as CE
using CommunityEcology: Community, Ordination, RDAResult, PhyloTree, BetaDisper, NMDSResult
import CommunityEcology: plot_ordination, plot_ordination!, plot_rda, plot_rda!,
    plot_composition, plot_composition!,
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

const MARKERS = (:circle, :utriangle, :diamond, :rect, :dtriangle, :cross, :star5, :hexagon)

# level => style, from `nothing` (cycle `default`), a vector (matched to the sorted levels),
# or a Dict keyed by the metadata values (or their strings)
function _styles(spec, levels, default)
    spec === nothing && return Dict(l => default[mod1(i, length(default))] for (i, l) in enumerate(levels))
    if spec isa AbstractDict
        s = Dict(string(k) => v for (k, v) in spec)
        return Dict(l => haskey(s, l) ? s[l] : throw(ArgumentError("no style given for \"$l\"")) for l in levels)
    end
    return Dict(l => spec[mod1(i, length(spec))] for (i, l) in enumerate(levels))
end

function _column(metadata, col)
    col in propertynames(metadata) || throw(ArgumentError("metadata has no column $col"))
    return string.(metadata[!, col])
end

_levels(values) = sort!(unique(values))

# numeric columns go on a colormap, categorical ones get a level => colour Dict
function _point_colors(o, colorby, palette)
    colorby === nothing && return :teal, nothing
    v = _column(o.metadata, colorby)
    raw = o.metadata[!, colorby]
    palette === nothing && all(x -> !ismissing(x) && x isa Real, raw) && return Float64.(raw), nothing
    styles = _styles(palette, _levels(v), PALETTE)
    return [styles[x] for x in v], styles
end

"""
    ordinationplot(o::Ordination; colorby, markerby, palette, markers, alpha, centroids, ...)

Sample scores of `o`, coloured by the metadata column `colorby` and shaped by `markerby`.
`palette` and `markers` are vectors (matched to the sorted levels) or Dicts from level to
colour or marker; `alpha` applies to the sample points. With `centroids=true` (needs a
categorical `colorby`) each `colorby` group gets a centroid (`centroidmarker`), joined by
dashed lines to the centroid of each of its `markerby` subgroups or, without `markerby`, to
each of its points.
"""
Makie.@recipe OrdinationPlot (ordination,) begin
    colorby = nothing
    markerby = nothing
    palette = nothing
    markers = nothing
    alpha = 1.0
    labels = false
    markersize = 13
    centroids = false
    centroidsize = 22
    centroidmarker = :rect
    strokewidth = 1.5
end

function Makie.plot!(plot::OrdinationPlot)
    o = plot[1][]
    size(o.points, 2) >= 2 || throw(ArgumentError("ordination needs two axes to plot"))
    x, y = o.points[:, 1], o.points[:, 2]
    colors, cstyles = _point_colors(o, plot.colorby[], plot.palette[])
    a = plot.alpha[]
    pointcolors = colors isa AbstractVector{<:Real} ? colors :
                  colors isa AbstractVector ? [(c, a) for c in colors] : (colors, a)
    mb = plot.markerby[]
    mv = mb === nothing ? nothing : _column(o.metadata, mb)
    mstyles = mb === nothing ? nothing : _styles(plot.markers[], _levels(mv), MARKERS)
    Makie.scatter!(plot, x, y; color=pointcolors, colormap=:viridis, markersize=plot.markersize[],
                   marker=mb === nothing ? :circle : [mstyles[m] for m in mv])
    if plot.centroids[]
        cstyles === nothing && throw(ArgumentError("centroids need a categorical colorby"))
        cv = _column(o.metadata, plot.colorby[])
        centre(sel) = Makie.Point2f(sum(x[sel]) / count(sel), sum(y[sel]) / count(sel))
        cs, sw = plot.centroidsize[], plot.strokewidth[]
        for g in _levels(cv)
            sel = cv .== g
            c = centre(sel)
            if mv === nothing
                for i in findall(sel)
                    Makie.lines!(plot, [c, Makie.Point2f(x[i], y[i])]; linestyle=:dash, color=(cstyles[g], a))
                end
            else
                for m in _levels(mv[sel])
                    cc = centre(sel .& (mv .== m))
                    Makie.lines!(plot, [c, cc]; linestyle=:dash, color=cstyles[g])
                    Makie.scatter!(plot, cc; color=cstyles[g], strokewidth=sw, marker=mstyles[m], markersize=cs)
                end
            end
            Makie.scatter!(plot, c; color=cstyles[g], strokewidth=sw, marker=plot.centroidmarker[], markersize=cs)
        end
    end
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

Makie.@recipe RDAPlot (result,) begin
    labels = false
    showfeatures = true
    markersize = 11
end

function Makie.plot!(plot::RDAPlot)
    r = plot[1][]
    size(r.sites, 2) >= 2 || throw(ArgumentError("RDA needs two constrained axes to plot"))
    Makie.scatter!(plot, r.sites[:, 1], r.sites[:, 2];
                   color=:teal, markersize=plot.markersize[])
    if plot.labels[]
        Makie.text!(plot, r.sites[:, 1], r.sites[:, 2]; text=r.samples,
                    offset=(7, 7), align=(:left, :bottom), fontsize=10)
    end
    if plot.showfeatures[]
        site_extent = maximum(abs, r.sites[:, 1:2])
        feature_extent = maximum(abs, r.feature_scores[:, 1:2])
        scale = iszero(feature_extent) ? 0.0 : 0.8 * site_extent / feature_extent
        for (i, feature) in enumerate(r.features)
            x = scale * r.feature_scores[i, 1]
            y = scale * r.feature_scores[i, 2]
            Makie.lines!(plot, [0.0, x], [0.0, y]; color=:darkorange, linewidth=1.5)
            Makie.text!(plot, x, y; text=feature, color=:darkorange,
                        offset=(4, 4), align=(:left, :bottom), fontsize=10)
        end
    end
    return plot
end

plot_rda!(axis, r::RDAResult; kwargs...) = rdaplot!(axis, r; kwargs...)

function plot_rda(r::RDAResult; size=(800, 550), kwargs...)
    Base.size(r.sites, 2) >= 2 || throw(ArgumentError("RDA needs two constrained axes to plot"))
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1];
                    xlabel="RDA1 ($(round(100 * r.explained[1]; digits=1))%)",
                    ylabel="RDA2 ($(round(100 * r.explained[2]; digits=1))%)",
                    title="RDA · adjusted R² $(round(r.adjusted_r2; digits=3))")
    plot_rda!(ax, r; kwargs...)
    return fig
end

function _composition(c, top, proportions)
    top >= 1 || throw(ArgumentError("top must be positive"))
    reduced = CE.lump_rare(c, top)
    X = proportions ? CE.relative(reduced).counts : Float64.(reduced.counts)
    return reduced, X
end

# sample order from `sortby`: a column, `column => level order`, or a vector of these (outermost first)
function _sample_order(metadata, sortby)
    sortby === nothing && return collect(1:size(metadata, 1))
    specs = sortby isa AbstractVector ? sortby : [sortby]
    key(s::Symbol) = _column(metadata, s)
    function key(s::Pair)
        v, levels = _column(metadata, first(s)), string.(last(s))
        pos = indexin(v, levels)
        any(isnothing, pos) && throw(ArgumentError("$(first(s)) has values not in the given level order"))
        return pos
    end
    keys = [key(s) for s in specs]
    return sortperm(collect(zip(keys...)))         # stable: keeps the original order within groups
end

# runs of equal consecutive values, as (value, first, last)
function _runs(v)
    runs = Tuple{String,Int,Int}[]
    for (i, x) in enumerate(v)
        if isempty(runs) || runs[end][1] != x
            push!(runs, (x, i, i))
        else
            runs[end] = (x, runs[end][2], i)
        end
    end
    return runs
end

"""
    compositionplot(c::Community; top, features, colors, sortby, above, below, ...)

Stacked bars of the feature abundances in each sample. `features` picks and orders the
features (otherwise the `top` most abundant, with the rest lumped), coloured by `colors`.
`sortby` orders the samples by metadata columns (`column` or `column => level order`,
outermost first). `above` and `below` name metadata columns whose consecutive groups are
marked by labelled brackets above or below the bars; `grouplabels` and `groupcolors` are
Dicts from group value to label text and bracket colour.
"""
Makie.@recipe CompositionPlot (community,) begin
    top = 10
    proportions = true
    features = nothing
    colors = nothing
    sortby = nothing
    above = nothing
    below = nothing
    grouplabels = nothing
    groupcolors = nothing
    barwidth = 0.86
    bracketwidth = 1.5
    fontsize = 14
end

function Makie.plot!(plot::CompositionPlot)
    c0 = plot[1][]
    order = _sample_order(c0.metadata, plot.sortby[])
    fs = plot.features[]
    if fs === nothing
        c, X = _composition(c0[order, :], plot.top[], plot.proportions[])
        names = c.features
    else
        base = plot.proportions[] ? CE.relative(c0) : c0
        cols = indexin(string.(fs), base.features)
        any(isnothing, cols) && throw(ArgumentError("features not in the community: $(fs[isnothing.(cols)])"))
        X, names = Float64.(base.counts[order, cols]), string.(fs)
    end
    meta = c0.metadata[order, :]
    palette = plot.colors[] === nothing ? [_color(j) for j in eachindex(names)] : Makie.to_color.(plot.colors[])
    w = plot.barwidth[]
    rectangles, colors = Makie.Rect2f[], eltype(palette)[]
    for i in axes(X, 1)
        bottom = 0.0
        for j in axes(X, 2)
            if X[i, j] > 0
                push!(rectangles, Makie.Rect2f(i - w / 2, bottom, w, X[i, j]))
                push!(colors, palette[mod1(j, length(palette))])
            end
            bottom += X[i, j]
        end
    end
    isempty(rectangles) || Makie.poly!(plot, rectangles; color=colors, strokewidth=0)
    h = plot.proportions[] ? 1.0 : maximum(sum(X; dims=2))
    labels = plot.grouplabels[] === nothing ? Dict{String,String}() :
             Dict(string(k) => string(v) for (k, v) in plot.grouplabels[])
    bcolors = plot.groupcolors[] === nothing ? Dict{String,Any}() :
              Dict{String,Any}(string(k) => v for (k, v) in plot.groupcolors[])
    for (col, y0, y1, ylab) in ((plot.above[], 1.02, 1.05, 1.09), (plot.below[], -0.02, -0.06, -0.14))
        col === nothing && continue
        for (g, s, e) in _runs(_column(meta, col))
            Makie.lines!(plot, [s - 0.4, s - 0.4, e + 0.4, e + 0.4], h .* [y0, y1, y1, y0];
                         color=get(bcolors, g, :black), linewidth=plot.bracketwidth[])
            Makie.text!(plot, (s + e) / 2, ylab * h; text=get(labels, g, g), align=(:center, :center),
                        fontsize=plot.fontsize[])
        end
    end
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
