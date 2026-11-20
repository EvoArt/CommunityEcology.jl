# Saving and loading communities as plain files that R, Excel or a text
# editor can open: counts.csv, taxonomy.csv, metadata.csv, sequences.fasta
# and tree.nwk in one folder.

"""
    write_community(dir, c)

Save `c` as plain files in the folder `dir`: `counts.csv` (one row per
sample), `taxonomy.csv` and `sequences.fasta` (one entry per feature),
`metadata.csv` and `tree.nwk`, when present. In R,
`read.csv("counts.csv", row.names = 1)` gives the matrix phyloseq's
`otu_table` expects.
"""
function write_community(dir::AbstractString, c::Community)
    mkpath(dir)
    cnt = DataFrame(c.counts, c.features)
    insertcols!(cnt, 1, :sample => c.samples)
    CSV.write(joinpath(dir, "counts.csv"), cnt)
    CSV.write(joinpath(dir, "metadata.csv"), c.metadata)
    if c.taxonomy !== nothing
        tax = copy(c.taxonomy)
        insertcols!(tax, 1, :feature => c.features)
        CSV.write(joinpath(dir, "taxonomy.csv"), tax)
    end
    if !isempty(c.sequences)
        open(joinpath(dir, "sequences.fasta"), "w") do io
            for (f, s) in zip(c.features, c.sequences)
                println(io, ">", f, "\n", s)
            end
        end
    end
    c.tree === nothing || write_newick(joinpath(dir, "tree.nwk"), c.tree)
    return dir
end

"""
    read_community(dir)

Read a community saved with [`write_community`](@ref).
"""
function read_community(dir::AbstractString)
    f = joinpath(dir, "counts.csv")
    isfile(f) || throw(ArgumentError("no counts.csv in $dir"))
    cnt = CSV.read(f, DataFrame)
    samples = string.(cnt[!, 1])
    feats = names(cnt)[2:end]
    X = Matrix(cnt[!, 2:end])
    tax = isfile(joinpath(dir, "taxonomy.csv")) ? select(CSV.read(joinpath(dir, "taxonomy.csv"), DataFrame), Not(1)) : nothing
    meta = isfile(joinpath(dir, "metadata.csv")) ? CSV.read(joinpath(dir, "metadata.csv"), DataFrame) : nothing
    seqs = String[]
    if isfile(joinpath(dir, "sequences.fasta"))
        d = Dict{String,String}()
        name = ""
        for l in eachline(joinpath(dir, "sequences.fasta"))
            startswith(l, ">") ? (name = strip(l[2:end]); d[name] = "") : (d[name] *= strip(l))
        end
        seqs = [d[f] for f in feats]
    end
    tr = isfile(joinpath(dir, "tree.nwk")) ? read_newick(joinpath(dir, "tree.nwk")) : nothing
    if meta !== nothing
        meta[!, :sample] = string.(meta[!, :sample])
    end
    return Community(X; samples=samples, features=feats, sequences=seqs, taxonomy=tax, metadata=meta, tree=tr)
end
