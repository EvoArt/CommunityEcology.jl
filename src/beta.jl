# Beta diversity: labelled distance matrices between samples.

"""
    DistanceMatrix

A symmetric matrix of dissimilarities between samples, with their `labels`
and the `method` used. It is an `AbstractMatrix`, so it can be passed to
`permanova`, `betadisper`, `nmds` or `pcoa` directly.
"""
struct DistanceMatrix <: AbstractMatrix{Float64}
    D::Matrix{Float64}
    labels::Vector{String}
    method::Symbol
end
Base.size(d::DistanceMatrix) = size(d.D)
Base.getindex(d::DistanceMatrix, i::Int, j::Int) = d.D[i, j]
Base.getindex(d::DistanceMatrix, s::AbstractString, t::AbstractString) =
    d.D[findfirst(==(s), d.labels), findfirst(==(t), d.labels)]
Base.Matrix(d::DistanceMatrix) = copy(d.D)

const _BETA = (:braycurtis, :jaccard, :sorensen, :ruzicka, :euclidean,
               :manhattan, :canberra, :gower, :chord, :hellinger, :aitchison,
               :unifrac, :wunifrac)

"""
    beta_diversity(c, method=:braycurtis; pseudocount=0.5)

Dissimilarities between the samples of `c`:

- `:braycurtis`: Bray-Curtis on counts, Σ|xᵢ - yᵢ| / Σ(xᵢ + yᵢ)
  (use `relative(c)` first to weight samples equally)
- `:jaccard`, `:sorensen`: on presence/absence
- `:ruzicka`: quantitative Jaccard (as `vegan::vegdist(method="jaccard")`)
- `:manhattan`, `:canberra` (mean over non-double-zero features), `:gower`
  (feature-range-scaled mean absolute difference)
- `:euclidean`, `:chord`, `:hellinger` (Euclidean on transformed data),
  `:aitchison` (Euclidean on the centred log-ratio, with `pseudocount`)
- `:unifrac`, `:wunifrac`: unweighted and weighted (normalised) UniFrac
  (Lozupone & Knight 2005; Lozupone et al. 2007), which need a tree

Samples with no reads are rejected: drop them with [`prune`](@ref).
"""
function beta_diversity(c::Community, method::Symbol=:braycurtis; pseudocount::Real=0.5)
    method in _BETA || throw(ArgumentError("unknown method :$method; choose from " * join(_BETA, ", ")))
    nsamples(c) > 0 || throw(ArgumentError("beta diversity needs at least one sample"))
    any(==(0), sample_depths(c)) && throw(ArgumentError("some samples have no reads; drop them first with prune(c)"))
    X = float.(c.counts)
    D = if method === :braycurtis
        pairwise(BrayCurtis(), X; dims=1)
    elseif method === :jaccard
        _binary_dist(X .> 0, :jaccard)
    elseif method === :sorensen
        _binary_dist(X .> 0, :sorensen)
    elseif method === :ruzicka
        _abundance_dist(X, :ruzicka)
    elseif method === :euclidean
        pairwise(Euclidean(), X; dims=1)
    elseif method === :manhattan
        pairwise(Cityblock(), X; dims=1)
    elseif method === :canberra
        _abundance_dist(X, :canberra)
    elseif method === :gower
        _abundance_dist(X, :gower)
    elseif method === :chord
        pairwise(Euclidean(), chord(c).counts; dims=1)
    elseif method === :hellinger
        pairwise(Euclidean(), hellinger(c).counts; dims=1)
    elseif method === :aitchison
        pairwise(Euclidean(), clr(c; pseudocount=pseudocount).counts; dims=1)
    else
        c.tree === nothing && throw(ArgumentError(":$method needs a tree: c = build_tree(c)"))
        unifrac(c; weighted=method === :wunifrac)
    end
    for i in axes(D, 1)
        D[i, i] = 0.0
    end
    return DistanceMatrix(D, copy(c.samples), method)
end

function _abundance_dist(X, method)
    n, p = size(X)
    D = zeros(n, n)
    ranges = method === :gower ? vec(maximum(X; dims=1) .- minimum(X; dims=1)) : Float64[]
    for i in 1:n, j in i+1:n
        a = view(X, i, :)
        b = view(X, j, :)
        d = if method === :ruzicka
            den = sum(max.(a, b))
            den == 0 ? 0.0 : sum(abs.(a .- b)) / den
        elseif method === :canberra
            nz = 0
            s = 0.0
            for f in 1:p
                den = a[f] + b[f]
                if den > 0
                    s += abs(a[f] - b[f]) / den
                    nz += 1
                end
            end
            nz == 0 ? 0.0 : s / nz
        else
            p == 0 ? 0.0 : sum(f -> ranges[f] == 0 ? 0.0 :
                               abs(a[f] - b[f]) / ranges[f], 1:p) / p
        end
        D[i, j] = D[j, i] = d
    end
    return D
end

function _binary_dist(B, kind)
    n = size(B, 1)
    D = zeros(n, n)
    for i in 1:n, j in i+1:n
        a = sum(B[i, :] .& B[j, :])
        bc = sum(B[i, :] .⊻ B[j, :])
        d = kind === :jaccard ? (a + bc == 0 ? 0.0 : bc / (a + bc)) : (bc / (2a + bc))
        D[i, j] = D[j, i] = d
    end
    return D
end

"""
    unifrac(c; weighted=false)

UniFrac distances between the samples of `c`, using its tree. Unweighted:
the fraction of branch length leading to features in one sample but not
both. Weighted (normalised, as `phyloseq::UniFrac(weighted=TRUE)` and
`GUniFrac`'s `d_1`): branch lengths weighted by differences in relative
abundance below them.
"""
function unifrac(c::Community; weighted::Bool=false)
    t = c.tree
    t === nothing && throw(ArgumentError("unifrac needs a tree"))
    tips = Dict(n => i for (i, n) in enumerate(tiplabels(t)))
    ns = nsamples(c)
    # A[v, s]: share of sample s's reads below node v
    A = zeros(length(t.parent), ns)
    P = relative(c).counts
    for (j, f) in enumerate(c.features)
        v = tips[f]
        for s in 1:ns
            A[v, s] += P[s, j]
        end
    end
    order = postorder(t)
    for v in order
        p = t.parent[v]
        p == 0 && continue
        for s in 1:ns
            A[p, s] += A[v, s]
        end
    end
    edges = [v for v in order if t.parent[v] != 0]
    b = t.edge
    D = zeros(ns, ns)
    if weighted
        dep = _depths(t)
        for i in 1:ns, j in i+1:ns
            num = sum(b[v] * abs(A[v, i] - A[v, j]) for v in edges)
            den = sum(dep[k] * (A[k, i] + A[k, j]) for k in 1:t.ntips)
            D[i, j] = D[j, i] = den > 0 ? num / den : 0.0
        end
    else
        for i in 1:ns, j in i+1:ns
            num = 0.0; den = 0.0
            for v in edges
                hi, hj = A[v, i] > 0, A[v, j] > 0
                (hi || hj) || continue
                den += b[v]
                hi != hj && (num += b[v])
            end
            D[i, j] = D[j, i] = den > 0 ? num / den : 0.0
        end
    end
    return D
end

function Base.show(io::IO, ::MIME"text/plain", d::DistanceMatrix)
    n = size(d, 1)
    _title(io, "DistanceMatrix", "$n samples, $(d.method)")
    n == 0 && return
    off = [d.D[i, j] for j in 1:n for i in j+1:n]
    isempty(off) || _kv(io, "range", _g(minimum(off)) * " – " * _g(maximum(off)) * "  (mean " * _g(mean(off)) * ")")
    k = min(n, 6)
    w = max(6, maximum(length, d.labels[1:k]))
    _section(io, "Dissimilarities", min(_termwidth(io), 100))
    println(io, "  ", " "^w, join(lpad.(first.(d.labels[1:k], 7), 8)))
    for i in 1:k
        print(io, "  ", rpad(d.labels[i], w), join(lpad.(_g.(d.D[i, 1:k]), 8)))
        i < k && println(io)
    end
    n > k && print(io, "\n", _note(io, "  … $(n - k) more samples"))
end
Base.show(io::IO, d::DistanceMatrix) = print(io, "DistanceMatrix($(size(d, 1)) samples, $(d.method))")
