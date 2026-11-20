# Ordination of communities, and the statistics packages applied to them.

"""
    Ordination

Sample coordinates from [`ordinate`](@ref), with the sample metadata. Fields:
`points` (samples × axes), `samples`,
`metadata`, `method` (`:nmds` or `:pcoa`), `distance`, `axislabels`, and the
underlying `result` (an `NMDSResult` or `PCoAResult`).
"""
struct Ordination
    points::Matrix{Float64}
    samples::Vector{String}
    metadata::DataFrame
    method::Symbol
    distance::Symbol
    axislabels::Vector{String}
    result::Any
end

"""
    ordinate(c, method=:nmds; distance=:braycurtis, k=2, kw...)
    ordinate(D::DistanceMatrix, c, method=:nmds; k=2, kw...)

Ordinate the samples of `c`: non-metric MDS (`:nmds`, several random
starts; keywords go to `nmds`) or principal coordinates (`:pcoa`), on the
dissimilarity `distance` (see [`beta_diversity`](@ref)) or a given
`DistanceMatrix`.
"""
function ordinate(c::Community, method::Symbol=:nmds; distance::Symbol=:braycurtis, k::Integer=2, kw...)
    return ordinate(beta_diversity(c, distance), c, method; k=k, kw...)
end

function ordinate(D::DistanceMatrix, c::Community, method::Symbol=:nmds; k::Integer=2, kw...)
    D.labels == c.samples || throw(ArgumentError("the distance matrix and the community have different samples"))
    if method === :nmds
        r = nmds(D.D; k=k, labels=D.labels, kw...)
        labs = ["NMDS$i" for i in 1:size(r.points, 2)]
    elseif method === :pcoa
        r = pcoa(D.D; k=k, labels=D.labels)
        labs = [@sprintf("PCo%d (%.1f%%)", i, 100e) for (i, e) in enumerate(r.explained)]
    else
        throw(ArgumentError("method must be :nmds or :pcoa"))
    end
    return Ordination(r.points, copy(c.samples), c.metadata, method, D.method, labs, r)
end

function Base.show(io::IO, ::MIME"text/plain", o::Ordination)
    _title(io, "Ordination", "$(uppercase(string(o.method))) of $(length(o.samples)) samples, $(o.distance)")
    if o.result isa NMDSResult
        _kv(io, "stress", @sprintf("%.4f", o.result.stress) * "  (best of $(o.result.tries) starts, reached $(o.result.repeats)×)")
    else
        _kv(io, "axes", join(o.axislabels, ", "))
    end
    n = length(o.samples)
    k = _nshow(io, n, 10)
    k == 0 && return
    dims = min(size(o.points, 2), 3)
    rows = [[o.samples[i], [_g(o.points[i, j]) for j in 1:dims]...] for i in 1:k]
    _section(io, "Sample coordinates", min(_termwidth(io), 100))
    _table(io, ["sample", o.axislabels[1:dims]...], rows;
           align=[:l, fill(:r, dims)...])
    _more(io, n, k, "samples")
end
Base.show(io::IO, o::Ordination) = print(io, "Ordination($(o.method), $(length(o.samples)) samples)")

# ------------------------------------------------- statistics on communities

"""
    permanova(c::Community, formula; distance=:braycurtis, kw...)

PERMANOVA of the samples of `c` (dissimilarity `distance`, see
[`beta_diversity`](@ref)) on its metadata. Keywords go to
`PERMANOVA.permanova`: `permutations`, `by`, `strata`, `rng`.
"""
PERMANOVA.permanova(c::Community, formula::PERMANOVA.StatsModels.FormulaTerm; distance::Symbol=:braycurtis, kw...) =
    permanova(c.metadata, beta_diversity(c, distance).D, formula; kw...)

"""
    pairwise_permanova(c::Community, variable::Symbol; distance=:braycurtis, kw...)

PERMANOVA between every pair of levels of the metadata column `variable`.
"""
PERMANOVA.pairwise_permanova(c::Community, variable::Symbol; distance::Symbol=:braycurtis, kw...) =
    pairwise_permanova(c.metadata, beta_diversity(c, distance).D, variable; kw...)

"""
    betadisper(c::Community, variable::Symbol; distance=:braycurtis, kw...)

Distances of samples to their group's spatial median, grouping by the
metadata column `variable`. Test with `permutest`.
"""
BetaDispersion.betadisper(c::Community, variable::Symbol; distance::Symbol=:braycurtis, kw...) =
    betadisper(beta_diversity(c, distance).D, c.metadata[!, variable]; labels=c.samples, kw...)

# The package's `DistanceMatrix` overloads
NonmetricMDS.nmds(D::DistanceMatrix; kw...) = nmds(D.D; labels=D.labels, kw...)
NonmetricMDS.pcoa(D::DistanceMatrix; kw...) = pcoa(D.D; labels=D.labels, kw...)
