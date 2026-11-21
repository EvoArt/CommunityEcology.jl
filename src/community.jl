"""
    Community

A community data set: a table of counts with one row per sample and one
column per feature (ASV, OTU or taxon), with optional feature sequences,
a taxonomy (one row per feature, one column per rank), sample metadata
(one row per sample, first column `sample`) and a phylogenetic tree.

Construct one from a sequence table (DenoiseDNA's `SequenceTable`), from a
count matrix, or read one saved with [`write_community`](@ref):

```julia
c = Community(seqtab; taxonomy=taxa, metadata="samples.csv")
c = Community(counts; samples=names, features=ids)
```

Index it like a matrix, by position, `Bool` vector or name:
`c[1:5, :]`, `c[:, ["ASV1", "ASV2"]]`, `c[c.metadata.site .== "A", :]`.
"""
struct Community{T<:Real}
    counts::Matrix{T}
    samples::Vector{String}
    features::Vector{String}
    sequences::Vector{String}            # empty, or one per feature
    taxonomy::Union{Nothing,DataFrame}   # one row per feature: ranks only
    metadata::DataFrame                  # one row per sample, column :sample first
    tree::Union{Nothing,PhyloTree}       # tips named by feature
end

"""
    Community(counts; samples, features, sequences=String[], taxonomy=nothing,
              metadata=nothing, tree=nothing, by=:sample)
    Community(seqtab; taxonomy=nothing, metadata=nothing, tree=nothing, prefix="ASV", by=:sample)

`counts` is a samples × features matrix. A sequence table from DenoiseDNA
(anything with `samples`, `sequences` and `table` fields) gives the counts,
sample names and sequences; features are named `ASV1`, `ASV2`, … in table
order.

- `taxonomy`: any table with one row per feature (e.g. the `TaxonomyTable`
  from AssignTaxonomy). If it has a `Sequence` column, rows are matched to
  features by sequence; otherwise they must be in feature order.
- `metadata`: a table or a CSV/TSV file with one row per sample, matched by
  the column `by` (default `:sample`). Samples without metadata get `missing`.
- `tree`: a [`PhyloTree`](@ref) whose tips are the features.
"""
function Community(counts::AbstractMatrix{<:Real}; samples=string.("S", 1:size(counts, 1)),
                   features=string.("F", 1:size(counts, 2)), sequences=String[], taxonomy=nothing,
                   metadata=nothing, tree=nothing, by::Symbol=:sample)
    ns, nf = size(counts)
    samples = String.(string.(collect(samples)))      # String, not e.g. CSV's String7
    features = String.(string.(collect(features)))
    sequences = String.(string.(collect(sequences)))
    length(samples) == ns || throw(DimensionMismatch("$(length(samples)) sample names for $ns rows"))
    length(features) == nf || throw(DimensionMismatch("$(length(features)) feature names for $nf columns"))
    isempty(sequences) || length(sequences) == nf ||
        throw(DimensionMismatch("$(length(sequences)) sequences for $nf features"))
    allunique(samples) || throw(ArgumentError("sample names must be unique"))
    allunique(features) || throw(ArgumentError("feature names must be unique"))
    any(<(0), counts) && throw(ArgumentError("counts must be non-negative"))
    tax = taxonomy === nothing ? nothing : _match_taxonomy(taxonomy, features, sequences)
    meta = _match_metadata(metadata, samples, by)
    tree === nothing || _check_tree(tree, features)
    return Community(Matrix(counts), samples, features, sequences, tax, meta, tree)
end

"""
    Community(table; sample=:sample, features=nothing, metadata=nothing, kw...)

A community from a wide table (a `DataFrame`, CSV file, or any Tables.jl
table) with one row per sample: the column `sample` names the samples and
the `features` columns (default: every other column) hold the counts.

```julia
c = Community("counts.csv"; metadata="samples.csv")
```
"""
function Community(df::AbstractDataFrame; sample::Symbol=:sample, features=nothing, kw...)
    String(sample) in names(df) ||
        throw(ArgumentError("the table has no column `$sample` naming the samples; its columns are: " *
                            join(first(names(df), 8), ", ") * (ncol(df) > 8 ? ", …" : "") *
                            ". Pass `sample=:column`."))
    cols = features === nothing ? setdiff(names(df), [String(sample)]) : string.(collect(features))
    bad = [n for n in cols if !(eltype(df[!, n]) <: Union{Missing,Real})]
    isempty(bad) || throw(ArgumentError("count columns must be numeric; these are not: " * join(bad, ", ")))
    X = Matrix(coalesce.(df[!, cols], 0))
    X = all(isinteger, X) ? Int.(X) : float.(X)
    return Community(X; samples=string.(df[!, sample]), features=cols, kw...)
end

Community(file::AbstractString; delim=nothing, kw...) = Community(read_metadata(file; delim=delim); kw...)

function Community(st; prefix="ASV", kw...)
    if !all(p -> hasproperty(st, p), (:samples, :sequences, :table)) && Tables.istable(st)
        return Community(DataFrame(Tables.columntable(st)); kw...)
    end
    all(p -> hasproperty(st, p), (:samples, :sequences, :table)) ||
        throw(ArgumentError("expected a count matrix, a table of counts, or a sequence table (with samples, sequences and table)"))
    nf = length(st.sequences)
    return Community(st.table; samples=st.samples, features=string.(prefix, 1:nf), sequences=st.sequences, kw...)
end

function _match_taxonomy(tax, features, sequences)
    df = DataFrame(Tables.columntable(tax))
    nf = length(features)
    if "Sequence" in names(df)
        isempty(sequences) && throw(ArgumentError("the taxonomy has a Sequence column but the community has no sequences"))
        pos = Dict(s => i for (i, s) in enumerate(df.Sequence))
        idx = [get(pos, s, 0) for s in sequences]
        nmiss = count(==(0), idx)
        nmiss > 0 && @warn "$nmiss feature(s) have no row in the taxonomy; their ranks are missing"
        df = select(df, Not(:Sequence))
        out = DataFrame([n => Vector{Union{Missing,String}}(missing, nf) for n in names(df)])
        for (i, k) in enumerate(idx)
            k == 0 && continue
            for n in names(df)
                v = df[k, n]
                out[i, n] = ismissing(v) ? missing : string(v)
            end
        end
        return out
    end
    nrow(df) == nf || throw(DimensionMismatch("taxonomy has $(nrow(df)) rows for $nf features (and no Sequence column to match on)"))
    return DataFrame([n => Union{Missing,String}[ismissing(v) ? missing : string(v) for v in df[!, n]]
                      for n in names(df)])
end

function _match_metadata(meta, samples, by)
    meta === nothing && return DataFrame(sample=copy(samples))
    df = meta isa AbstractString ? read_metadata(meta) : DataFrame(Tables.columntable(meta))
    if !(String(by) in names(df))
        throw(ArgumentError("metadata has no column `$by` to match samples on; its columns are: " *
                            join(names(df), ", ") * ". Pass `by=:column`."))
    end
    key = string.(df[!, by])
    allunique(key) || throw(ArgumentError("metadata column `$by` has duplicated sample names"))
    pos = Dict(k => i for (i, k) in enumerate(key))
    idx = [get(pos, s, 0) for s in samples]
    missing_s = samples[idx .== 0]
    if !isempty(missing_s)
        @warn "$(length(missing_s)) sample(s) have no metadata row: " * join(first(missing_s, 5), ", ") *
              (length(missing_s) > 5 ? ", …" : "")
    end
    extra = length(key) - count(>(0), idx)
    extra > 0 && @info "$extra metadata row(s) match no sample and were ignored"
    out = DataFrame(sample=copy(samples))
    for n in names(df)
        n == String(by) && continue
        col = df[!, n]
        out[!, n] = [i == 0 ? missing : col[i] for i in idx]
        any(ismissing, out[!, n]) || (out[!, n] = disallowmissing(out[!, n]))
    end
    return out
end

function _check_tree(t, features)
    labs = Set(tiplabels(t))
    miss = count(f -> !(f in labs), features)
    miss == 0 || throw(ArgumentError("$miss feature(s) are not tips of the tree"))
    return t
end

# ---------------------------------------------------------------- accessors

nsamples(c::Community) = size(c.counts, 1)
nfeatures(c::Community) = size(c.counts, 2)
Base.size(c::Community) = size(c.counts)
Base.size(c::Community, d) = size(c.counts, d)
"Sample names."
samples(c::Community) = c.samples
"Feature (ASV/taxon) names."
features(c::Community) = c.features
"The samples × features count matrix."
counts(c::Community) = c.counts
"Sample metadata (a `DataFrame`, one row per sample)."
metadata(c::Community) = c.metadata
"Taxonomy (a `DataFrame`, one row per feature), or `nothing`."
taxonomy(c::Community) = c.taxonomy
"Feature sequences (empty if unknown)."
sequences(c::Community) = c.sequences
"The phylogenetic tree, or `nothing`."
tree(c::Community) = c.tree
"Taxonomic ranks available, as `Symbol`s."
ranks(c::Community) = c.taxonomy === nothing ? Symbol[] : Symbol.(names(c.taxonomy))
"Total count (sequencing depth) of each sample."
sample_depths(c::Community) = vec(sum(c.counts; dims=2))

"""
    add_metadata(c, meta; by=:sample)

A copy of `c` with sample metadata from `meta` (a table or a CSV/TSV file),
replacing any it had.
"""
add_metadata(c::Community, meta; by::Symbol=:sample) =
    Community(c.counts, c.samples, c.features, c.sequences, c.taxonomy, _match_metadata(meta, c.samples, by), c.tree)

"""
    read_metadata(file; delim=nothing)

Read a CSV or TSV sample sheet into a `DataFrame` (the delimiter is
detected unless given).
"""
read_metadata(file::AbstractString; delim=nothing) =
    delim === nothing ? CSV.read(file, DataFrame) : CSV.read(file, DataFrame; delim=delim)

# ---------------------------------------------------------------- indexing

_index(::Colon, names) = Colon()
_index(i::Union{Integer,AbstractVector{<:Integer},AbstractVector{Bool},AbstractRange}, names) = i
_index(s::AbstractString, names) = _index([s], names)
function _index(v::AbstractVector{<:AbstractString}, names)
    pos = Dict(n => i for (i, n) in enumerate(names))
    return [haskey(pos, s) ? pos[s] : throw(KeyError(s)) for s in v]
end
_rows(i::Integer) = [i]
_rows(i) = i

function Base.getindex(c::Community, i, j)
    ii = _rows(_index(i, c.samples))
    jj = _rows(_index(j, c.features))
    seqs = isempty(c.sequences) ? String[] : c.sequences[jj]
    tax = c.taxonomy === nothing ? nothing : c.taxonomy[jj, :]
    f = c.features[jj]
    tr = c.tree === nothing ? nothing : (jj isa Colon ? c.tree : keep_tips(c.tree, f))
    return Community(c.counts[ii, jj], c.samples[ii], f, seqs, tax, c.metadata[ii, :], tr)
end

"""
    keep_tips(t::PhyloTree, names)

The subtree spanned by the tips `names`, with unary nodes removed.
"""
function keep_tips(t::PhyloTree, names)
    want = Set(string.(names))
    labs = tiplabels(t)
    all(in(Set(labs)), want) || throw(ArgumentError("some names are not tips of the tree"))
    length(want) == t.ntips && return t
    keep = falses(nnodes(t))
    for v in 1:t.ntips
        labs[v] in want || continue
        u = v
        while u != 0 && !keep[u]
            keep[u] = true
            u = t.parent[u]
        end
    end
    # tips first (in the original order), then internal nodes
    order = vcat(findall(v -> keep[v], 1:t.ntips), findall(v -> keep[v], t.ntips+1:nnodes(t)) .+ t.ntips)
    newid = zeros(Int, nnodes(t))
    newid[order] = 1:length(order)
    parent = [t.parent[v] == 0 ? 0 : newid[t.parent[v]] for v in order]
    sub = _tree(parent, t.edge[order], t.names[order], count(v -> keep[v], 1:t.ntips))
    # drop unary nodes, including a unary root
    sub = _drop_unary(sub)
    while length(sub.children[sub.root]) == 1
        r = sub.root
        c = only(sub.children[r])
        p = copy(sub.parent); p[c] = 0; p[r] = -1
        keepn = findall(>=(0), p)
        nid = zeros(Int, length(p)); nid[keepn] = 1:length(keepn)
        sub = _tree([p[v] <= 0 ? 0 : nid[p[v]] for v in keepn], sub.edge[keepn], sub.names[keepn], sub.ntips)
    end
    return sub
end

Base.:(==)(a::Community, b::Community) =
    a.counts == b.counts && a.samples == b.samples && a.features == b.features &&
    a.sequences == b.sequences && isequal(a.taxonomy, b.taxonomy) && isequal(a.metadata, b.metadata)

"""
    cophenetic(t::PhyloTree)
    cophenetic(c::Community)

Patristic (tip-to-tip) distances along the branches of a tree. For a
community, rows and columns follow `features(c)`.
"""
function cophenetic(c::Community)
    c.tree === nothing && throw(ArgumentError("the community has no tree"))
    pos = Dict(n => i for (i, n) in enumerate(tiplabels(c.tree)))
    idx = [pos[f] for f in c.features]
    return cophenetic(c.tree)[idx, idx]
end

"""
    add_taxonomy(c, taxonomy)

A copy of `c` with a new taxonomy: any table with one row per feature, such
as AssignTaxonomy's result. Rows are matched to features by a `Sequence`
column when there is one, otherwise taken in feature order.
"""
add_taxonomy(c::Community, tax) =
    Community(c.counts, c.samples, c.features, c.sequences, _match_taxonomy(tax, c.features, c.sequences),
              c.metadata, c.tree)

"""
    add_tree(c, tree)

A copy of `c` with a phylogenetic tree: a [`PhyloTree`](@ref), Newick text,
or a Newick file. Every feature must be a tip; other tips are dropped.
"""
function add_tree(c::Community, t)
    t isa PhyloTree || (t = read_newick(t))
    _check_tree(t, c.features)
    return Community(c.counts, c.samples, c.features, c.sequences, c.taxonomy, c.metadata,
                     keep_tips(t, c.features))
end

"""
    add_taxonomy(c, reference::AbstractString; kw...)

Classify the community's sequences against a DADA2-formatted reference file with
AssignTaxonomy's naive Bayesian classifier and attach the ranks; keywords go to
`AssignTaxonomy.assignTaxonomy`. Needs `using AssignTaxonomy`. Pipes well:

```julia
@chain seqtab begin
    Community
    add_taxonomy("silva_nr99_v138.2_toGenus_trainset.fa.gz")
    collapse(:Family)
end
```
"""
add_taxonomy(c::Community, reference::AbstractString; kw...) = _classify(c, reference; kw...)
_classify(c, reference; kw...) =
    error("add_taxonomy(community, reference_file) classifies with AssignTaxonomy: run `using AssignTaxonomy` first")
