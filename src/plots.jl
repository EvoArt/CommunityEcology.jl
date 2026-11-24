# Makie recipes live in the optional CommunityEcologyMakieExt extension.
# Load CairoMakie, GLMakie, or WGLMakie before calling these functions.

"""
    plot_ordination(o::Ordination; colorby, labels, size, kwargs...) -> Figure
    plot_ordination!(axis, o::Ordination; kwargs...)

Plot the sample scores of `o`. Load a Makie backend first. Keywords:

- `colorby`: metadata column for the colour (a colormap if numeric).
- `markerby`: metadata column for the marker.
- `palette`, `markers`: colours and markers, as vectors (matched to the sorted levels) or
  Dicts from level to colour or marker.
- `alpha`: opacity of the sample points; `markersize`; `labels` (sample names).
- `centroids=true` (needs a categorical `colorby`): a centroid for each `colorby` group
  (`centroidmarker`, `centroidsize`), joined by dashed lines to the centroid of each of its
  `markerby` subgroups or, without `markerby`, to each of its points.

To plot saved coordinates, wrap them in `Ordination(points, metadata)`.
"""
function plot_ordination end
function plot_ordination! end

"""
    plot_rda(r::RDAResult; size=(800, 550), labels=false, showfeatures=true, markersize=11) -> Figure
    plot_rda!(axis, r::RDAResult; kwargs...)

Plot RDA site scores on the first two constrained axes, with feature arrows if
`showfeatures`. `labels` adds sample names. Needs two constrained axes.
"""
function plot_rda end
function plot_rda! end

"""
    plot_composition(c::Community; top, proportions, size, kwargs...) -> Figure
    plot_composition!(axis, c::Community; kwargs...)

Plot the composition of each sample as stacked bars. Load a Makie backend first. Keywords:

- `top=10`: plot the `top` most abundant features and lump the rest, unless `features`
  names the features to plot, in order; `colors` colours them.
- `proportions=true`: relative abundances (over all features) rather than counts.
- `sortby`: order the samples by metadata columns, each `column` or `column => level order`,
  outermost first (stable, so the original order is kept within groups).
- `above`, `below`: metadata columns whose runs of consecutive samples get labelled brackets
  above or below the bars; `grouplabels` and `groupcolors` are Dicts from a group value to
  its label and bracket colour.
- `barwidth`, `bracketwidth`, `fontsize`.
"""
function plot_composition end
function plot_composition! end

"""
    plot_alpha(c::Community; measure=:shannon, group=nothing, size=(760, 470), markersize=11) -> Figure
    plot_alpha!(axis, c::Community; kwargs...)

Plot one [`alpha_diversity`](@ref) measure per sample, or jittered by metadata
column `group` with group means.
"""
function plot_alpha end
function plot_alpha! end

"""
    plot_rarefaction(c::Community; step=nothing, linewidth=2, size=(760, 500)) -> Figure
    plot_rarefaction!(axis, c::Community; kwargs...)

Plot [`rarecurve`](@ref) expected richness against depth, one line per sample.
Legend shown for ≤ 12 samples.
"""
function plot_rarefaction end
function plot_rarefaction! end

"""
    plot_heatmap(c::Community; proportions=true, colormap=:viridis, size=(850, 550)) -> Figure
    plot_heatmap!(axis, c::Community; kwargs...)

Plot a samples × features heatmap of relative abundances (`proportions=true`)
or counts.
"""
function plot_heatmap end
function plot_heatmap! end

"""
    plot_dispersion(bd::BetaDisper; size=(700, 450), markersize=11) -> Figure
    plot_dispersion!(axis, bd::BetaDisper; kwargs...)

Plot each sample's distance to its group centre, with group means.
"""
function plot_dispersion end
function plot_dispersion! end

"""
    plot_shepard(r; size=(700, 480), markersize=5) -> Figure
    plot_shepard!(axis, r; kwargs...)

Shepard diagram (ordination distance vs observed dissimilarity, with the fitted
monotone curve). `r` is an `NMDSResult` or an NMDS [`Ordination`](@ref).
"""
function plot_shepard end
function plot_shepard! end

"""
    plot_tree(t; size=(800, 550), linecolor=:teal, linewidth=2, labels=true) -> Figure
    plot_tree!(axis, t; kwargs...)

Plot a rooted `PhyloTree` as a rectangular cladogram with branch lengths.
`t` may be a [`Community`](@ref) with a tree. `labels` shows tip names.
"""
function plot_tree end
function plot_tree! end
