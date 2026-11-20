"""
    CommunityEcology

Community data for ecologists: counts, diversity, beta diversity, ordination,
phylogenetic analyses, PERMANOVA, and tests of dispersion.
"""
module CommunityEcology

using DataFrames, Distances, LinearAlgebra, Random, Statistics, Tables
import CSV, Downloads, p7zip_jll
import DataFrames: metadata   # extend DataAPI.metadata, so `using DataFrames` does not clash
using Printf: @sprintf
using SpecialFunctions: loggamma

import PERMANOVA
import BetaDispersion
using PERMANOVA: permanova, pairwise_permanova, PermanovaResult, PairwisePermanova, adjust_p, @formula
using BetaDispersion: betadisper, permutest, dispersion_anova, BetaDisper, DispersionTest

import NonmetricMDS
using NonmetricMDS: nmds, pcoa, NMDSResult, PCoAResult, procrustes, shepard

export Community, nsamples, nfeatures, samples, features, counts, metadata, taxonomy, sequences,
       tree, ranks, sample_depths, add_metadata, add_taxonomy, add_tree, read_metadata
export filter_samples, filter_features, prune, rarefy, relative, hellinger, chord, wisconsin, clr, presence_absence,
       collapse, top_features, lump_rare
export richness, shannon, simpson, invsimpson, chao1, pielou, hill, alpha_diversity, rarecurve, specaccum
export DistanceMatrix, beta_diversity, unifrac
export PhyloTree, nj_tree, build_tree, sequence_distances, midpoint_root, write_newick, read_newick,
       faith_pd, ntips, tiplabels, cophenetic
export Ordination, ordinate, RDAResult, rda
export MantelResult, mantel, envfit, factorfit
export simper
export write_community, read_community
export dune, varespec, miseq_sop, miseq_sop_reads
export permanova, pairwise_permanova, PermanovaResult, PairwisePermanova, adjust_p, @formula
export betadisper, permutest, dispersion_anova, BetaDisper, DispersionTest
export nmds, pcoa, NMDSResult, PCoAResult, procrustes, shepard
export plot_ordination, plot_ordination!, plot_composition, plot_composition!,
       plot_rda, plot_rda!,
       plot_alpha, plot_alpha!, plot_rarefaction, plot_rarefaction!,
       plot_heatmap, plot_heatmap!, plot_dispersion, plot_dispersion!,
       plot_shepard, plot_shepard!, plot_tree, plot_tree!

include("style.jl")
include("tree.jl")
include("community.jl")
include("transform.jl")
include("alpha.jl")
include("beta.jl")
include("align.jl")
include("ordination.jl")
include("rda.jl")
include("permutation.jl")
include("simper.jl")
include("io.jl")
include("show.jl")
include("plots.jl")
include("datasets.jl")

end
