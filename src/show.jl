# Compact REPL summaries use the same colour-aware helpers as DenoiseDNA.jl.
# `show(io, MIME"text/plain"(), x)` is deliberately capped even without :limit.

Base.show(io::IO, c::Community) =
    print(io, "Community($(nsamples(c)) samples × $(nfeatures(c)) features)")

function _composition_strip(io, counts, order, depth)
    b = IOBuffer()
    for j in order
        x = counts[j]
        glyph = x == 0 ? "·" : x / depth >= 0.1 ? "█" : x / depth >= 0.01 ? "▓" : "░"
        print(b, x == 0 ? _sty(io, glyph; color=_PAL.tail) :
              _sty(io, glyph; color=_PAL.spark))
    end
    return String(take!(b))
end

function Base.show(io::IO, ::MIME"text/plain", c::Community)
    ns, nf = size(c)
    total = sum(c.counts)
    _title(io, "Community", "$(_fmt(ns)) samples × $(_fmt(nf)) features" *
           _note(io, " · $(_fmt(total)) counts"))
    _checks(io, [(c.taxonomy !== nothing, "taxonomy"),
                 (ncol(c.metadata) > 1, "metadata"),
                 (!isempty(c.sequences), "sequences"),
                 (c.tree !== nothing, "tree")])
    (ns == 0 || nf == 0) && return
    depths = sample_depths(c)
    present = vec(sum(c.counts .> 0; dims=2))
    _kv(io, "depth", "$(_fmt(minimum(depths)))–$(_fmt(maximum(depths))) per sample" *
                     _note(io, " · median $(_fmt(median(depths)))"))
    _kv(io, "richness", "$(_fmt(minimum(present)))–$(_fmt(maximum(present))) features per sample")
    if c.taxonomy !== nothing
        _kv(io, "ranks", join(names(c.taxonomy), ", "))
    end
    W = _termwidth(io)
    order = sortperm(vec(sum(c.counts; dims=1)); rev=true)
    nshow = _nshow(io, ns, 10)
    nw = min(18, max(6, maximum(length, c.samples[1:nshow])))
    nstrip = min(nf, max(0, W - nw - 30))
    rows = Vector{String}[]
    for i in 1:nshow
        name = length(c.samples[i]) > nw ? first(c.samples[i], nw - 1) * "…" : c.samples[i]
        strip = nstrip == 0 || depths[i] == 0 ? "" :
                _composition_strip(io, view(c.counts, i, :), view(order, 1:nstrip), depths[i])
        push!(rows, [name, _fmt(depths[i]), _fmt(present[i]), strip])
    end
    _section(io, "Samples", min(W, 100))
    _table(io, ["sample", "counts", "features", "abundance →"], rows;
           align=[:l, :r, :r, :l])
    _more(io, ns, nshow, "samples")
    totals = vec(sum(c.counts; dims=1))
    prev = vec(sum(c.counts .> 0; dims=1))
    k = _nshow(io, nf, 8)
    _section(io, "Most abundant features", min(W, 100))
    frows = Vector{String}[]
    for j in order[1:k]
        name = c.features[j]
        name = length(name) > 28 ? first(name, 27) * "…" : name
        push!(frows, [name, _fmt(totals[j]), total == 0 ? "—" : _pct(totals[j] / total),
                      _fmt(prev[j])])
    end
    _table(io, ["feature", "counts", "share", "samples"], frows;
           align=[:l, :r, :r, :r])
    _more(io, nf, k, "features")
end

Base.show(io::IO, m::MantelResult) =
    print(io, "MantelResult(r=$(_g(m.r)), p=$(_g(m.pvalue)))")

function Base.show(io::IO, ::MIME"text/plain", m::MantelResult)
    _title(io, "MantelResult", "$(m.method) correlation of distance matrices")
    _kv(io, "r", _g(m.r))
    _kv(io, "p", _g(m.pvalue) * _note(io, " · $(m.permutations) permutations"))
    _kv(io, "alternative", string(m.alternative))
end

Base.show(io::IO, r::RDAResult) =
    print(io, "RDAResult($(length(r.samples)) samples, $(length(r.features)) features, R²=$(_g(r.r2)))")

function Base.show(io::IO, ::MIME"text/plain", r::RDAResult)
    _title(io, "RDAResult", "$(length(r.samples)) samples × $(length(r.features)) features" *
           _note(io, " · $(r.transform) response"))
    _kv(io, "predictors", join(string.(r.predictors), ", "))
    _kv(io, "R²", _pct(r.r2) * _note(io, " · adjusted " * _pct(r.adjusted_r2)))
    _kv(io, "axes", join((_pct(x) for x in r.explained), ", ") * " of total inertia")
    n = length(r.samples)
    k = _nshow(io, n, 10)
    k == 0 && return
    dims = min(size(r.sites, 2), 3)
    rows = [[r.samples[i], [_g(r.sites[i, j]) for j in 1:dims]...] for i in 1:k]
    _section(io, "Constrained site scores", min(_termwidth(io), 100))
    _table(io, ["sample", ["RDA$j" for j in 1:dims]...], rows;
           align=[:l, fill(:r, dims)...])
    _more(io, n, k, "samples")
end
