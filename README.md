# CommunityEcology.jl

Community ecology analyses for Julia. A `Community` holds a samples × features
count table, sample metadata, and optional taxonomy, sequences, and a
phylogenetic tree. The package covers a useful subset of the workflows commonly
done with R's [vegan](https://vegandevs.github.io/vegan/): diversity, distance
matrices, ordination, permutation tests, and transformations. It does not aim
for API or numerical equivalence with every vegan function.

## Quick start

```julia
using CommunityEcology, DataFrames, Random

counts = [8 1 0; 7 2 1; 0 1 9; 1 0 8]
meta = DataFrame(sample=["A1", "A2", "B1", "B2"],
                 site=["A", "A", "B", "B"],
                 temperature=[10.0, 12.0, 20.0, 23.0])
c = Community(counts; samples=meta.sample,
              features=["oak", "pine", "birch"], metadata=meta)

alpha_diversity(c)                  # richness, Shannon, Simpson, Chao1, etc.
d = beta_diversity(c, :braycurtis)
o = ordinate(d, c, :pcoa)           # or :nmds
rda(c, [:site, :temperature])       # constrained ordination
permanova(c, @formula(1 ~ site); permutations=999)
bd = betadisper(c, :site)
permutest(bd; permutations=999)
envfit(o; variables=[:temperature], permutations=999)
mantel(d, beta_diversity(c, :jaccard); permutations=999)
simper(c, :site)                  # descriptive feature contributions
```

The examples above use four samples to show the API. Use larger datasets for
permutation tests; the number of distinct group arrangements limits attainable
p-values.

## Main operations

| Task | CommunityEcology.jl | Related vegan operation |
| --- | --- | --- |
| Alpha diversity | `alpha_diversity`, `richness`, `shannon`, `simpson`, `chao1`, `hill` | `diversity`, `specnumber`, `estimateR` |
| Rarefaction | `rarefy`, `rarecurve` | `rrarefy`, `rarecurve` |
| Species accumulation | `specaccum` | `specaccum(method="exact")` |
| Transformations | `relative`, `hellinger`, `chord`, `wisconsin`, `clr`, `presence_absence` | `decostand`, `wisconsin` |
| Dissimilarities | `beta_diversity`, `unifrac` | a subset of `vegdist` |
| Ordination | `ordinate`, `nmds`, `pcoa`, `rda` | `metaMDS`, `wcmdscale`, basic `rda` |
| Group comparison | `permanova`, `pairwise_permanova`, `betadisper`, `permutest` | `adonis2`, `betadisper`, `permutest` |
| Association | `mantel`, `envfit`, `factorfit` | `mantel`, `envfit` |
| Feature contributions | `simper` | `simper` (descriptive breakdown) |
| Phylogeny | `read_newick`, `nj_tree`, `faith_pd`, `unifrac` | related phylogenetic workflows |

Samples are rows and features are columns. Metadata is matched by a `sample`
column; see `?Community` for taxonomy, sequence-table, and tree inputs. All
transformations return new `Community` values. `write_community` and
`read_community` round-trip counts and annotations through CSV, FASTA, and
Newick files.

`nmds` and `pcoa` come from the separate **NonmetricMDS.jl** package;
PERMANOVA and dispersion methods likewise use their separate packages.
Load `CairoMakie`, `GLMakie`, or `WGLMakie` to enable the optional Makie
recipes. For example:

```julia
using CairoMakie
fig = plot_ordination(o; colorby=:site, labels=true)
save("ordination.png", fig)
plot_composition(c; top=10)
plot_alpha(c; measure=:shannon, group=:site)
plot_rarefaction(c)
plot_heatmap(c)
```

Each plotting function also has an `!` form for an existing Makie `Axis`.
`plot_dispersion`, `plot_shepard`, and `plot_tree` cover their respective
analysis results.

`rda` currently supports numeric and categorical metadata predictors, with
Hellinger, relative-abundance, or raw response data. `envfit` fits numeric
vectors and categorical centroids. CCA, conditioning terms, and ordination
permutation tests are not yet provided.

`beta_diversity(c, :jaccard)` uses presence/absence. Vegan's default
`vegdist(method="jaccard")` uses quantitative abundances; choose `:ruzicka`
for that measure here. `simper` reports descriptive contributions and does
not provide vegan's per-feature permutation p-values.

## Scope relative to vegan

The supported workflows cover common count-table preparation, alpha diversity,
rarefaction, several ecological dissimilarities, unconstrained ordination,
basic RDA, PERMANOVA, dispersion, environmental fitting, Mantel testing,
SIMPER, species accumulation, and phylogenetic diversity.

Substantial vegan areas remain: CCA and correspondence analysis, distance-based
RDA, partial ordination, constrained-ordination permutation tests and model
selection, null models, variation partitioning, additional dissimilarity
indices, and advanced ordination plotting. Results from the Julia implementations
should be checked against vegan for any workflow that requires numerical
parity, especially across different scaling or permutation conventions.

## Development

This checkout uses locally developed dependencies. After those packages are
available in your Julia environment, run:

```julia
using Pkg
Pkg.develop(path="path/to/CommunityEcology")
Pkg.test("CommunityEcology")
```
