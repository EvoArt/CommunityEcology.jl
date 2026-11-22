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

"""Plot site scores and feature directions from a two-axis `RDAResult`."""
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

"""Plot one alpha diversity measure by sample or metadata group."""
function plot_alpha end
function plot_alpha! end

"""Plot exact expected rarefaction curves for the samples of a community."""
function plot_rarefaction end
function plot_rarefaction! end

"""Plot a sample-by-feature abundance heatmap."""
function plot_heatmap end
function plot_heatmap! end

"""Plot distances to group centres from a `BetaDisper` result."""
function plot_dispersion end
function plot_dispersion! end

"""Plot a Shepard diagram for an `NMDSResult` or an NMDS `Ordination`."""
function plot_shepard end
function plot_shepard! end

"""Plot a rooted `PhyloTree` as a rectangular tree."""
function plot_tree end
function plot_tree! end
