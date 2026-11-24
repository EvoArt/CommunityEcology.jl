# Alpha diversity of one sample (a vector of counts or abundances), and of
# every sample of a community. Definitions follow vegan's `diversity`,
# `specnumber`, `estimateR` and `rarefy`.

_props(x) = (s = sum(x); s > 0 ? [v / s for v in x if v > 0] : Float64[])

"""
    richness(x) -> Int

Number of entries of `x` greater than zero.
"""
richness(x::AbstractVector) = count(>(0), x)

"""
    shannon(x) -> Float64

Shannon index `-Σ pᵢ log pᵢ` (natural log), as `vegan::diversity`.
"""
shannon(x::AbstractVector) = -sum((p * log(p) for p in _props(x)); init=0.0)

"""
    simpson(x) -> Float64

Gini-Simpson index `1 - Σ pᵢ²`, as `vegan::diversity(index="simpson")`.
"""
simpson(x::AbstractVector) = 1 - sum(abs2, _props(x); init=0.0)

"""
    invsimpson(x) -> Float64

Inverse Simpson index `1 / Σ pᵢ²`; `0.0` for an empty sample.
"""
invsimpson(x::AbstractVector) = (s = sum(abs2, _props(x); init=0.0); s > 0 ? 1 / s : 0.0)

"""
    pielou(x) -> Float64

Pielou's evenness `shannon(x) / log(richness(x))`; `NaN` if richness ≤ 1.
"""
pielou(x::AbstractVector) = (S = richness(x); S > 1 ? shannon(x) / log(S) : NaN)

"""
    chao1(x)

Bias-corrected Chao1 estimate of richness, `S + F₁(F₁ - 1) / (2(F₂ + 1))`,
from integer counts (F₁ singletons, F₂ doubletons), as `vegan::estimateR`.
Denoised ASV tables contain few singletons, which makes Chao1 close to the
observed richness there.
"""
function chao1(x::AbstractVector)
    S = richness(x)
    F1 = count(==(1), x)
    F2 = count(==(2), x)
    return S + F1 * (F1 - 1) / (2 * (F2 + 1))
end

"""
    hill(x, q)

Hill number of order `q`, the effective number of features: `q = 0`
richness, `q = 1` exp(Shannon), `q = 2` inverse Simpson.
"""
function hill(x::AbstractVector, q::Real)
    p = _props(x)
    isempty(p) && return 0.0
    q == 0 && return float(length(p))
    q == 1 && return exp(-sum(v * log(v) for v in p))
    return sum(v^q for v in p)^(1 / (1 - q))
end

const _ALPHA = Dict(:richness => richness, :shannon => shannon, :simpson => simpson,
                    :invsimpson => invsimpson, :chao1 => chao1, :pielou => pielou)

"""
    alpha_diversity(c; measures=[:richness, :shannon, :simpson, :invsimpson, :chao1, :pielou],
                    metadata=true)

A `DataFrame` with one row per sample: `sample`, `depth`, the requested
measures and (with `metadata=true`) the sample metadata, ready for plotting
or modelling. With a tree, `:faith_pd` is also available.

Richness depends on sequencing depth: compare it across samples after
[`rarefy`](@ref), or look at [`rarecurve`](@ref).
"""
function alpha_diversity(c::Community; measures=[:richness, :shannon, :simpson, :invsimpson, :chao1, :pielou],
                         metadata::Bool=true)
    df = DataFrame(sample=c.samples, depth=sample_depths(c))
    for m in measures
        if m === :faith_pd
            c.tree === nothing && throw(ArgumentError(":faith_pd needs a tree (see build_tree)"))
            tips = Dict(n => i for (i, n) in enumerate(tiplabels(c.tree)))
            pos = [tips[f] for f in c.features]
            df[!, m] = [begin
                present = falses(ntips(c.tree))
                present[pos[view(c.counts, i, :) .> 0]] .= true
                faith_pd(c.tree, present)
            end for i in 1:nsamples(c)]
        else
            f = get(_ALPHA, m) do
                throw(ArgumentError("unknown measure :$m; choose from " * join(sort!(collect(keys(_ALPHA))), ", ") * ", faith_pd"))
            end
            df[!, m] = [f(view(c.counts, i, :)) for i in 1:nsamples(c)]
        end
    end
    if metadata
        for n in names(c.metadata)
            n == "sample" && continue
            df[!, n] = c.metadata[!, n]
        end
    end
    return df
end

"""
    rarecurve(c; step=nothing, depths=nothing)

Rarefaction curves: the expected number of features in a random subsample
of each depth (Hurlbert 1971; `vegan::rarefy`), computed exactly rather
than by resampling. Returns a long `DataFrame` with `sample`, `depth` and
`richness`, plus the sample metadata. `step` defaults to 1/50 of the
largest depth.
"""
function rarecurve(c::Community; step=nothing, depths=nothing)
    all(isinteger, c.counts) || throw(ArgumentError("rarecurve needs integer counts"))
    N = sample_depths(c)
    maxN = maximum(N)
    st = step === nothing ? max(1, round(Int, maxN / 50)) : step
    rows = DataFrame(sample=String[], depth=Int[], richness=Float64[])
    for i in 1:nsamples(c)
        x = round.(Int, filter(>(0), c.counts[i, :]))
        n = Int(N[i])
        ds = depths === nothing ? unique(vcat(1, st:st:n, n)) : filter(<=(n), collect(depths))
        for d in ds
            push!(rows, (c.samples[i], d, _expected_richness(x, n, d)))
        end
    end
    meta = select(c.metadata, :sample, Not(:sample))
    return leftjoin(rows, meta; on=:sample, order=:left)
end

# E[S_d] = Σ 1 - C(N - Nᵢ, d) / C(N, d), with log-binomials
_lchoose(n, k) = loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1)
function _expected_richness(x, N, d)
    d >= N && return float(length(x))
    lc = _lchoose(N, d)
    return sum(xi -> N - xi < d ? 1.0 : 1.0 - exp(_lchoose(N - xi, d) - lc), x)
end

"""
    specaccum(c::Community)

Expected number of features observed after pooling 1 through `nsamples(c)`
samples in random order. This is the exact, combinatorial expectation used by
`vegan::specaccum(method="exact")`; no random subsampling is required. Returns
a `DataFrame` with `samples` and `richness` columns. Empty samples contribute
no new features.
"""
function specaccum(c::Community)
    n = nsamples(c)
    n == 0 && return DataFrame(samples=Int[], richness=Float64[])
    prevalence = vec(sum(c.counts .> 0; dims=1))
    out = DataFrame(samples=collect(1:n), richness=zeros(n))
    for k in 1:n
        denom = _lchoose(n, k)
        out.richness[k] = sum(prevalence) do f
            f == 0 ? 0.0 : n - f < k ? 1.0 : -expm1(_lchoose(n - f, k) - denom)
        end
    end
    return out
end
