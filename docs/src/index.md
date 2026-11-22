```@meta
CurrentModule = CommunityEcology
```

# CommunityEcology.jl

CommunityEcology.jl analyses community count tables in Julia. Samples are rows and features are columns. A `Community` also holds sample metadata and can include taxonomy, sequences, and a phylogenetic tree.

## Create and analyse a community

Install with Julia's package manager (`] add CommunityEcology`) after registration. Until then, add the repository URL and the [NonmetricMDS.jl](https://github.com/EvoArt/NonmetricMDS.jl) dependency.

```julia
using CommunityEcology, DataFrames

counts = [8 1 0; 7 2 1; 0 1 9; 1 0 8]
meta = DataFrame(sample=["A1", "A2", "B1", "B2"], site=["A", "A", "B", "B"])
c = Community(counts; samples=meta.sample,
              features=["oak", "pine", "birch"], metadata=meta)

alpha_diversity(c)                 # one row per sample
d = beta_diversity(c, :braycurtis) # pairwise dissimilarities
ordinate(d, c, :pcoa)              # principal coordinates
```

For a multisite Jaccard decomposition, use `beta_partition(c)` or group by a metadata column:

```julia
beta_partition(c; family=:jaccard)
beta_partition(c, :site; family=:sorensen)
```

The decomposition returns total dissimilarity, turnover, and nestedness. The grouped result also reports within-group dispersion. Each group needs at least two samples.

## Main workflows

- `alpha_diversity`, `rarefy`, and `specaccum` describe within-sample diversity and accumulation.
- `beta_diversity`, `ordinate`, and `rda` analyse differences among samples.
- `permanova`, `betadisper`, `mantel`, and `envfit` test associations and group differences.
- `relative`, `hellinger`, and `presence_absence` transform counts into new `Community` values.
- `write_community` and `read_community` round-trip a community through files.

`nmds` and `pcoa` are provided by NonmetricMDS.jl; PERMANOVA and dispersion methods likewise use separate packages. Load a Makie backend to enable the optional plotting recipes. The [README](https://github.com/EvoArt/CommunityEcology.jl#readme) has a longer example and describes differences from R's vegan package.

## API reference

```@autodocs
Modules = [CommunityEcology]
```
