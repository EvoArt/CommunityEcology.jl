module CommunityEcologyBIOMFormatExt

using CommunityEcology, DataFrames
import BIOMFormat
using BIOMFormat: BIOMTable

const RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

_isdna(s) = length(s) >= 20 && all(in("ACGTNacgtn"), s)

# "k__Bacteria" -> "Bacteria"; "g__" or "" -> missing
function _rank(x)
    s = strip(replace(x, r"^[a-zA-Z]__" => ""))
    return isempty(s) ? missing : String(s)
end

function _taxonomy(lineages)
    all(isempty, lineages) && return nothing
    n = maximum(length, lineages)
    cols = n <= length(RANKS) ? RANKS[1:n] : vcat(RANKS, ["Rank$i" for i in length(RANKS)+1:n])
    return DataFrame([c => [k <= length(l) ? _rank(l[k]) : missing for l in lineages]
                      for (k, c) in enumerate(cols)])
end

function _sample_metadata(t)
    t.sample_metadata === nothing && return nothing
    keys_ = sort!(unique!([k for d in t.sample_metadata for k in keys(d)]))
    isempty(keys_) && return nothing
    df = DataFrame(sample=t.samples)
    for k in keys_
        df[!, k] = [get(d, k, missing) for d in t.sample_metadata]
    end
    return df
end

"""
    Community(t::BIOMTable; top=nothing, prefix="ASV", metadata=nothing, kw...)

A community from a BIOM table, transposed to samples × features. When the
observation ids are DNA sequences (as in Deblur or DADA2 output) they are
kept as the feature sequences and the features are named `ASV1`, `ASV2`, …
in table order. Lineages become a taxonomy table, and sample metadata in the
file is used unless `metadata` is given.

`top = n` keeps only the `n` most abundant observations, which keeps a large
table (such as the Earth Microbiome Project's) small in memory.
"""
function CommunityEcology.Community(t::BIOMTable; top=nothing, prefix="ASV", metadata=nothing, kw...)
    if top !== nothing
        totals = vec(sum(BIOMFormat.counts(t); dims=2))
        t = t[sort(partialsortperm(totals, 1:min(top, length(totals)); rev=true)), :]
    end
    ids = BIOMFormat.observations(t)
    dna = all(_isdna, ids)
    return Community(Matrix(permutedims(BIOMFormat.counts(t)));
                     samples=BIOMFormat.samples(t),
                     features=dna ? string.(prefix, eachindex(ids)) : ids,
                     sequences=dna ? uppercase.(ids) : String[],
                     taxonomy=_taxonomy(BIOMFormat.taxonomy(t)),
                     metadata=metadata === nothing ? _sample_metadata(t) : metadata,
                     kw...)
end

end
