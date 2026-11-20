# Makie recipes live in the optional CommunityEcologyMakieExt extension.
# Load CairoMakie, GLMakie, or WGLMakie before calling these functions.

"""Plot community sample scores. Load a Makie backend first. Returns a `Figure`."""
function plot_ordination end
function plot_ordination! end

"""Plot relative composition as stacked bars. Load a Makie backend first."""
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
