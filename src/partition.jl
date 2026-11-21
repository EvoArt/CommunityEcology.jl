# Multi-site beta diversity and its partition into turnover and nestedness (Baselga 2010, 2012)

const _PARTITION_NAMES = (jaccard=("βJAC", "βJTU", "βJNE"), sorensen=("βSOR", "βSIM", "βSNE"))

"""
    beta_partition(c; family=:jaccard)

Multi-site β-diversity of the samples of `c`, from presence/absence, and its partition into
species turnover and nestedness (Baselga 2010, 2012; as `betapart::beta.multi`). With
`Sᵢ` the richness of sample `i`, `S_T` the total richness, and `bᵢⱼ` the number of species in
`i` but not `j`, summed over pairs of samples as `min = Σ min(bᵢⱼ, bⱼᵢ)` and
`max = Σ max(bᵢⱼ, bⱼᵢ)`, and `a = Σ Sᵢ - S_T`:

- `family=:jaccard`: total `βJAC = (min + max) / (a + min + max)`, turnover
  `βJTU = 2min / (a + 2min)`, nestedness `βJNE = βJAC - βJTU`;
- `family=:sorensen`: `βSOR = (min + max) / (2a + min + max)`, `βSIM = min / (a + min)`,
  `βSNE = βSOR - βSIM`.

Returns `(total, turnover, nestedness, relative_turnover)`, the last being turnover / total.
For one table per group of samples, see `beta_partition(c, groups...)`.

```julia
beta_partition(dune())          # (total = 0.88…, turnover = 0.81…, nestedness = 0.07…, …)
```
"""
function beta_partition(c::Community; family::Symbol=:jaccard)
    return _beta_multi(c.counts .> 0, family)
end

function _beta_multi(P::AbstractMatrix{Bool}, family::Symbol; strict::Bool=true)
    family in keys(_PARTITION_NAMES) || throw(ArgumentError("family must be :jaccard or :sorensen, got :$family"))
    n = size(P, 1)
    if n < 2         # undefined for one sample: an error on its own, NaN in a table of groups
        strict && throw(ArgumentError("multi-site β-diversity needs at least two samples"))
        return (total=NaN, turnover=NaN, nestedness=NaN, relative_turnover=NaN)
    end
    a = count(P) - count(vec(any(P; dims=1)))          # Σ Sᵢ - S_T
    smin = smax = 0
    for i in 1:n-1, j in i+1:n
        bij = count(k -> P[i, k] & !P[j, k], axes(P, 2))
        bji = count(k -> P[j, k] & !P[i, k], axes(P, 2))
        smin += min(bij, bji)
        smax += max(bij, bji)
    end
    if family === :jaccard
        total = (smin + smax) / (a + smin + smax)
        turnover = 2smin / (a + 2smin)
    else
        total = (smin + smax) / (2a + smin + smax)
        turnover = smin / (a + smin)
    end
    return (total=total, turnover=turnover, nestedness=total - turnover, relative_turnover=turnover / total)
end

"""
    BetaPartition

The table returned by `beta_partition(c, groups...)`: one row per group of samples, with the
multi-site β-diversity partition and, optionally, the dispersion of the group. A Tables.jl
table (`DataFrame(bp)` gives the values); it displays rounded to `bp.digits` places, as plain
text, HTML or Markdown (so it renders as a table in Quarto documents, Word included).
"""
struct BetaPartition
    table::DataFrame
    family::Symbol
    dispersion::Union{Nothing,Symbol}
    digits::Int
end

"""
    beta_partition(c, groups::Symbol...; family=:jaccard, dispersion=:braycurtis, digits=2)

The multi-site β-diversity partition (see `beta_partition(c)`) for each group of samples,
the groups being the combinations of the metadata columns `groups`, in order of appearance.
With `dispersion` set to a dissimilarity (see [`beta_diversity`](@ref)), each group's
dispersion is added: the mean distance of its samples to the group centroid in principal
coordinate space ([`betadisper`](@ref) with `type=:centroid`), computed on `c` as given; pass
`relative(c)` for dispersion of proportions. `dispersion=nothing` leaves it out.

```julia
beta_partition(relative(c), :transfer, :carbon)   # one row per transfer × carbon source
```
"""
function beta_partition(c::Community, group::Symbol, groups::Symbol...; family::Symbol=:jaccard,
                        dispersion::Union{Nothing,Symbol}=:braycurtis, digits::Integer=2)
    cols = (group, groups...)
    c.metadata === nothing && throw(ArgumentError("the community has no metadata"))
    for g in cols
        g in propertynames(c.metadata) || throw(ArgumentError("no metadata column :$g"))
    end
    key = [Tuple(c.metadata[i, g] for g in cols) for i in 1:nsamples(c)]
    levels = unique(key)
    P = c.counts .> 0
    disp = dispersion === nothing ? nothing :
           betadisper(beta_diversity(c, dispersion).D, string.(key); type=:centroid)
    tot, tu, ne = _PARTITION_NAMES[family]
    t = DataFrame([String(g) => [l[k] for l in levels] for (k, g) in enumerate(cols)])
    parts = [_beta_multi(P[key .== Ref(l), :], family; strict=false) for l in levels]
    t[!, tot] = getfield.(parts, :total)
    t[!, tu] = getfield.(parts, :turnover)
    t[!, ne] = getfield.(parts, :nestedness)
    t[!, "relative turnover"] = getfield.(parts, :relative_turnover)
    if disp !== nothing
        t[!, "dispersion"] = [mean(disp.distances[key .== Ref(l)]) for l in levels]
    end
    return BetaPartition(t, family, dispersion, digits)
end

Tables.istable(::Type{BetaPartition}) = true
Tables.columnaccess(::Type{BetaPartition}) = true
Tables.columns(bp::BetaPartition) = Tables.columns(bp.table)
DataFrames.DataFrame(bp::BetaPartition) = copy(bp.table)
Base.size(bp::BetaPartition) = size(bp.table)

function _rounded(bp::BetaPartition)
    t = copy(bp.table)
    for n in names(t)
        eltype(t[!, n]) <: AbstractFloat && (t[!, n] = round.(t[!, n]; digits=bp.digits))
    end
    return t
end

function _partition_title(bp::BetaPartition)
    fam = bp.family === :jaccard ? "Jaccard" : "Sørensen"
    d = bp.dispersion === nothing ? "" : "; dispersion: mean distance to centroid ($(bp.dispersion))"
    return "Multi-site β-diversity partition ($fam family, Baselga 2010)$d"
end

function Base.show(io::IO, ::MIME"text/plain", bp::BetaPartition)
    println(io, _partition_title(bp))
    show(io, _rounded(bp); summary=false, eltypes=false, show_row_number=false)
end
Base.show(io::IO, bp::BetaPartition) = show(io, MIME"text/plain"(), bp)
Base.show(io::IO, m::MIME"text/html", bp::BetaPartition) = show(io, m, _rounded(bp); summary=false, eltypes=false)

function Base.show(io::IO, ::MIME"text/markdown", bp::BetaPartition)
    t = _rounded(bp)
    cells = [string.(collect(row)) for row in eachrow(t)]
    # separator widths follow the contents, which pandoc uses as relative column widths
    widths = [maximum(length, vcat(n, [c[j] for c in cells])) for (j, n) in enumerate(names(t))]
    println(io, "| ", join(names(t), " | "), " |")
    println(io, "|", join([repeat("-", w + 2) for w in widths], "|"), "|")
    for c in cells
        println(io, "| ", join(c, " | "), " |")
    end
end
