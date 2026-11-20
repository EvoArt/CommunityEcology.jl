# Filtering and transforming communities. Every function returns a new
# Community; the input is never modified.

_with_counts(c::Community, X) = Community(X, c.samples, c.features, c.sequences, c.taxonomy, c.metadata, c.tree)

"""
    filter_samples(f, c)
    filter_samples(c, column => predicate)

Keep the samples whose metadata row satisfies `f(row)`, e.g.
`filter_samples(r -> r.day > 10, c)` or `filter_samples(c, :site => ==("A"))`.
"""
filter_samples(f, c::Community) = c[[f(r) === true for r in eachrow(c.metadata)], :]
filter_samples(c::Community, p::Pair{Symbol}) = c[[p.second(x) === true for x in c.metadata[!, p.first]], :]

"""
    filter_features(f, c)
    filter_features(c, rank => predicate)

Keep features whose taxonomy row satisfies `f(row)`, e.g.
`filter_features(c, :Kingdom => ==("Bacteria"))`. Missing ranks never match.
"""
function filter_features(f, c::Community)
    c.taxonomy === nothing && throw(ArgumentError("the community has no taxonomy"))
    return c[:, [f(r) === true for r in eachrow(c.taxonomy)]]
end
function filter_features(c::Community, p::Pair{Symbol})
    c.taxonomy === nothing && throw(ArgumentError("the community has no taxonomy"))
    col = c.taxonomy[!, p.first]
    return c[:, [!ismissing(x) && p.second(x) === true for x in col]]
end

"""
    prune(c; min_count=1, min_samples=1, min_depth=1)

Drop features seen fewer than `min_count` times in total or in fewer than
`min_samples` samples, then samples with fewer than `min_depth` reads.
"""
function prune(c::Community; min_count::Real=1, min_samples::Integer=1, min_depth::Real=1)
    tot = vec(sum(c.counts; dims=1))
    prev = vec(sum(c.counts .> 0; dims=1))
    c2 = c[:, (tot .>= min_count) .& (prev .>= min_samples)]
    return c2[sample_depths(c2) .>= min_depth, :]
end

"""
    relative(c)

Counts as proportions of each sample's total.
"""
relative(c::Community) = _with_counts(c, c.counts ./ max.(sum(c.counts; dims=2), eps()))

"""
    hellinger(c)

The Hellinger transformation: square roots of relative abundances
(Legendre & Gallagher 2001). Euclidean distances on it are Hellinger distances.
"""
hellinger(c::Community) = _with_counts(c, sqrt.(relative(c).counts))

"""
    chord(c)

Scale each sample to Euclidean length one. Euclidean distances on this
transformation are chord distances. Empty samples remain zero.
"""
function chord(c::Community)
    X = Float64.(c.counts)
    norms = sqrt.(sum(abs2, X; dims=2))
    return _with_counts(c, X ./ max.(norms, eps(Float64)))
end

"""
    wisconsin(c)

Wisconsin double standardization: divide each feature by its maximum across
samples, then divide each sample by its resulting total. All-zero features
and samples remain zero.
"""
function wisconsin(c::Community)
    X = Float64.(c.counts)
    isempty(X) && return _with_counts(c, X)
    maxima = maximum(X; dims=1)
    X ./= max.(maxima, eps(Float64))
    X ./= max.(sum(X; dims=2), eps(Float64))
    return _with_counts(c, X)
end

"""
    clr(c; pseudocount=0.5)

Centred log-ratio transformation, `log(x + pseudocount)` minus the sample's
mean log. Euclidean distances on it are Aitchison distances.
"""
function clr(c::Community; pseudocount::Real=0.5)
    L = log.(c.counts .+ pseudocount)
    return _with_counts(c, L .- mean(L; dims=2))
end

"""
    presence_absence(c)

1 where a feature was seen in a sample, 0 otherwise.
"""
presence_absence(c::Community) = _with_counts(c, Int.(c.counts .> 0))

"""
    rarefy(c, depth=minimum(sample_depths(c)); rng=Random.default_rng())

Subsample every sample to `depth` reads without replacement. Samples with
fewer reads are dropped (with a message), and features left with no reads
are removed. Rarefying discards data; use it for richness comparisons or
when a method requires equal depths.
"""
function rarefy(c::Community, depth::Integer=minimum(sample_depths(c)); rng::AbstractRNG=Random.default_rng())
    all(x -> isinteger(x), c.counts) || throw(ArgumentError("rarefy needs integer counts"))
    d = sample_depths(c)
    keep = findall(>=(depth), d)
    length(keep) < nsamples(c) &&
        @info "$(nsamples(c) - length(keep)) sample(s) with fewer than $depth reads dropped: " *
              join(c.samples[setdiff(1:nsamples(c), keep)], ", ")
    X = zeros(Int, length(keep), nfeatures(c))
    for (r, i) in enumerate(keep)
        _subsample!(view(X, r, :), round.(Int, view(c.counts, i, :)), depth, rng)
    end
    out = Community(X, c.samples[keep], c.features, c.sequences, c.taxonomy, c.metadata[keep, :], c.tree)
    return prune(out; min_count=1, min_depth=0)
end

# draw `depth` reads without replacement from the reads of `x`
function _subsample!(out, x, depth, rng)
    pool = Vector{Int}(undef, sum(x))
    k = 0
    for (j, n) in enumerate(x), _ in 1:n
        pool[k+=1] = j
    end
    fill!(out, 0)
    # partial Fisher-Yates: the first `depth` entries become a random sample
    for i in 1:depth
        r = rand(rng, i:length(pool))
        pool[i], pool[r] = pool[r], pool[i]
        out[pool[i]] += 1
    end
    return out
end

"""
    collapse(c, rank::Symbol)

Sum the features that share a taxon at `rank` (e.g. `:Genus`). A feature
not assigned at that rank is grouped as "unclassified X", where X is its
deepest assigned rank above. The result's features are the taxon names and
its taxonomy keeps the ranks down to `rank`. Sequences and tree are dropped.
"""
function collapse(c::Community, rank::Symbol)
    c.taxonomy === nothing && throw(ArgumentError("the community has no taxonomy"))
    rs = ranks(c)
    k = findfirst(==(rank), rs)
    k === nothing && throw(ArgumentError("no rank `$rank`; available: " * join(rs, ", ")))
    tax = c.taxonomy
    labels = [_taxon_label(tax, i, k) for i in 1:nfeatures(c)]
    groups = unique(labels)
    gi = Dict(g => i for (i, g) in enumerate(groups))
    X = zeros(eltype(c.counts), nsamples(c), length(groups))
    for (j, l) in enumerate(labels)
        X[:, gi[l]] .+= view(c.counts, :, j)
    end
    newtax = DataFrame([n => Vector{Union{Missing,String}}(missing, length(groups)) for n in names(tax)[1:k]])
    for (j, l) in enumerate(labels)
        g = gi[l]
        for n in names(tax)[1:k]
            ismissing(newtax[g, n]) && (newtax[g, n] = tax[j, n])
        end
    end
    newtax[!, names(tax)[k]] = groups
    o = sortperm(vec(sum(X; dims=1)); rev=true)
    return Community(X[:, o], c.samples, groups[o], String[], newtax[o, :], c.metadata, nothing)
end

function _taxon_label(tax, i, k)
    v = tax[i, k]
    ismissing(v) || return v
    for r in k-1:-1:1
        u = tax[i, r]
        ismissing(u) || return "unclassified " * u
    end
    return "unclassified"
end

"""
    top_features(c, n=10; by=:mean)

Names of the `n` features with the highest mean relative abundance across
samples (`by=:mean`), or the highest total count (`by=:total`).
"""
function top_features(c::Community, n::Integer=10; by::Symbol=:mean)
    score = by === :total ? vec(sum(c.counts; dims=1)) : vec(mean(relative(c).counts; dims=1))
    o = sortperm(score; rev=true)
    return c.features[o[1:min(n, length(o))]]
end

"""
    lump_rare(c, n=10; other="Other")

Keep the `n` most abundant features (see [`top_features`](@ref)) and sum
the rest into one feature named `other`. Handy before a composition plot.
"""
function lump_rare(c::Community, n::Integer=10; other::AbstractString="Other")
    top = top_features(c, n)
    nfeatures(c) <= n && return c
    keep = in(Set(top)).(c.features)
    ti = _index(top, c.features)
    X = hcat(c.counts[:, ti], sum(c.counts[:, .!keep]; dims=2))
    tax = nothing
    if c.taxonomy !== nothing
        tax = vcat(c.taxonomy[ti, :], DataFrame([n => Union{Missing,String}[other] for n in names(c.taxonomy)]))
    end
    return Community(X, c.samples, vcat(top, other), String[], tax, c.metadata, nothing)
end
