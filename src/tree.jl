# Rooted phylogenetic trees: neighbour joining, midpoint rooting, Newick.
#
# Nodes 1:ntips are the tips (in the order of `names`); internal nodes
# follow. `parent[root] == 0`; `edge[v]` is the length of the branch above v.

"""
    PhyloTree

A rooted tree. Tips are nodes `1:ntips(t)`, named `tiplabels(t)`. Fields:
`parent` (0 at the root), `edge` (length of the branch above each node),
`children`, `names` (tip labels; internal nodes have `""`), `root`.

Build one with [`build_tree`](@ref) (from sequences) or [`nj_tree`](@ref)
(from a distance matrix), or read one with [`read_newick`](@ref).
"""
struct PhyloTree
    parent::Vector{Int}
    edge::Vector{Float64}
    children::Vector{Vector{Int}}
    names::Vector{String}
    ntips::Int
    root::Int
end

ntips(t::PhyloTree) = t.ntips
tiplabels(t::PhyloTree) = t.names[1:t.ntips]
nnodes(t::PhyloTree) = length(t.parent)

function _tree(parent, edge, names, ntips)
    n = length(parent)
    ch = [Int[] for _ in 1:n]
    root = 0
    for v in 1:n
        p = parent[v]
        p == 0 ? (root = v) : push!(ch[p], v)
    end
    return PhyloTree(parent, edge, ch, names, ntips, root)
end

# nodes with every child before its parent
function postorder(t::PhyloTree)
    out = Int[]
    stack = [(t.root, false)]
    while !isempty(stack)
        v, done = pop!(stack)
        if done || isempty(t.children[v])
            push!(out, v)
        else
            push!(stack, (v, true))
            for c in reverse(t.children[v])
                push!(stack, (c, false))
            end
        end
    end
    return out
end

"""
    nj_tree(D; labels=string.(1:n), nonnegative=true)

Neighbour-joining tree (Saitou & Nei 1987; the algorithm of Studier &
Keppler 1988, as `ape::nj`) of the symmetric distance matrix `D`, rooted at
its midpoint. Negative branch lengths, which NJ can produce, are set to
zero unless `nonnegative=false`.
"""
function nj_tree(D; labels=string.(1:size(D, 1)), nonnegative::Bool=true, root::Bool=true)
    n = size(D, 1)
    n >= 2 || throw(ArgumentError("need at least two sequences to build a tree"))
    labels = string.(collect(labels))
    length(labels) == n || throw(DimensionMismatch("$(length(labels)) labels for $n sequences"))
    if n == 2
        h = D[1, 2] / 2
        return _tree([3, 3, 0], [h, h, 0.0], vcat(labels, ""), 2)
    end
    d = Matrix{Float64}(D)
    active = collect(1:n)            # matrix rows still in play
    node = collect(1:n)              # tree node of each row
    parent = zeros(Int, 2n - 2)
    edge = zeros(2n - 2)
    next = n
    r = vec(sum(d; dims=2))
    while length(active) > 3
        m = length(active)
        # pair minimising Q(i,j) = (m-2) d(i,j) - r(i) - r(j)
        best, bi, bj = Inf, 0, 0
        for a in 1:m-1, b in a+1:m
            i, j = active[a], active[b]
            q = (m - 2) * d[i, j] - r[i] - r[j]
            if q < best
                best, bi, bj = q, a, b
            end
        end
        i, j = active[bi], active[bj]
        li = d[i, j] / 2 + (r[i] - r[j]) / (2 * (m - 2))
        lj = d[i, j] - li
        next += 1
        parent[node[i]], edge[node[i]] = next, li
        parent[node[j]], edge[node[j]] = next, lj
        # the new node takes row i; row j leaves
        for k in active
            k == i || k == j || (d[i, k] = d[k, i] = (d[i, k] + d[j, k] - d[i, j]) / 2)
        end
        d[i, i] = 0.0
        node[i] = next
        deleteat!(active, bj)
        for k in active
            r[k] = sum(d[k, l] for l in active)
        end
    end
    # the last three rows join at one centre node
    a, b, c = active
    next += 1
    la = (d[a, b] + d[a, c] - d[b, c]) / 2
    lb = d[a, b] - la
    lc = d[a, c] - la
    for (k, l) in ((a, la), (b, lb), (c, lc))
        parent[node[k]], edge[node[k]] = next, l
    end
    resize!(parent, next); resize!(edge, next)
    nonnegative && (edge .= max.(edge, 0.0))
    t = _tree(parent, edge, vcat(labels, fill("", next - n)), n)
    return root ? midpoint_root(t) : t
end

# undirected adjacency with branch lengths
function _adjacency(t::PhyloTree)
    adj = [Tuple{Int,Float64}[] for _ in 1:nnodes(t)]
    for v in 1:nnodes(t)
        p = t.parent[v]
        p == 0 && continue
        push!(adj[v], (p, t.edge[v]))
        push!(adj[p], (v, t.edge[v]))
    end
    return adj
end

function _distances_from(adj, s)
    dist = fill(Inf, length(adj))
    prev = zeros(Int, length(adj))
    dist[s] = 0.0
    stack = [s]
    while !isempty(stack)
        v = pop!(stack)
        for (w, l) in adj[v]
            if isinf(dist[w])
                dist[w] = dist[v] + l
                prev[w] = v
                push!(stack, w)
            end
        end
    end
    return dist, prev
end

"""
    midpoint_root(t::PhyloTree)

Re-root `t` at the midpoint of its longest tip-to-tip path.
"""
function midpoint_root(t::PhyloTree)
    adj = _adjacency(t)
    nt = t.ntips
    d1, _ = _distances_from(adj, 1)
    a = argmax(d1[1:nt])
    da, prev = _distances_from(adj, a)
    b = argmax(da[1:nt])
    half = da[b] / 2
    # walk from b towards a until the midpoint is passed
    v = b
    while prev[v] != 0 && da[prev[v]] > half
        v = prev[v]
    end
    u = prev[v]                                   # the midpoint lies on edge u–v
    u == 0 && return t
    lv = da[v] - half                             # part of the edge below the new root
    lu = half - da[u]
    # new root node splits edge u–v; orient every edge away from it
    n = nnodes(t)
    root = n + 1
    names = vcat(t.names, "")
    nadj = vcat(adj, [Tuple{Int,Float64}[]])
    filter!(x -> x[1] != v, nadj[u]); filter!(x -> x[1] != u, nadj[v])
    push!(nadj[u], (root, lu)); push!(nadj[v], (root, lv))
    push!(nadj[root], (u, lu)); push!(nadj[root], (v, lv))
    parent = zeros(Int, n + 1)
    edge = zeros(n + 1)
    seen = falses(n + 1)
    seen[root] = true
    stack = [root]
    while !isempty(stack)
        x = pop!(stack)
        for (w, l) in nadj[x]
            seen[w] && continue
            seen[w] = true
            parent[w], edge[w] = x, l
            push!(stack, w)
        end
    end
    return _drop_unary(_tree(parent, edge, names, nt))
end

# remove internal nodes with a single child (left by re-rooting at an old root)
function _drop_unary(t::PhyloTree)
    parent, edge = copy(t.parent), copy(t.edge)
    for v in t.ntips+1:nnodes(t)
        if length(t.children[v]) == 1 && parent[v] != 0
            c = only(t.children[v])
            parent[c] = parent[v]
            edge[c] += edge[v]
            parent[v] = -1
        end
    end
    keep = findall(>=(0), parent)
    newid = zeros(Int, length(parent))
    newid[keep] = 1:length(keep)
    p2 = [parent[v] == 0 ? 0 : newid[parent[v]] for v in keep]
    return _tree(p2, edge[keep], t.names[keep], t.ntips)
end

"""
    write_newick([io or file,] t::PhyloTree)

The tree in Newick format (a string, or written to `io` or a file), for
FigTree, iTOL, R's `ape::read.tree` or Julia's Phylo.jl.
"""
function write_newick(t::PhyloTree)
    io = IOBuffer()
    _newick(io, t, t.root)
    print(io, ";")
    return String(take!(io))
end
write_newick(io::IO, t::PhyloTree) = print(io, write_newick(t))
write_newick(file::AbstractString, t::PhyloTree) = open(io -> println(io, write_newick(t)), file, "w")

function _newick(io, t, v)
    if !isempty(t.children[v])
        print(io, "(")
        for (k, c) in enumerate(t.children[v])
            k > 1 && print(io, ",")
            _newick(io, t, c)
        end
        print(io, ")")
    end
    nm = t.names[v]
    print(io, occursin(r"[\s(),:;\[\]']", nm) ? "'" * replace(nm, "'" => "''") * "'" : nm)
    t.parent[v] != 0 && print(io, ":", repr(t.edge[v]))
end

"""
    read_newick(text_or_file)

Read a tree in Newick format (from a string that starts with `(`, or a file).
"""
function read_newick(s::AbstractString)
    txt = strip(startswith(lstrip(s), "(") ? s : read(s, String))
    endswith(txt, ";") && (txt = txt[1:end-1])
    parent = Int[]; edge = Float64[]; name = String[]; istip = Bool[]
    stack = Int[]
    cur = 0
    i = firstindex(txt)
    function newnode()
        push!(parent, isempty(stack) ? 0 : stack[end]); push!(edge, 0.0); push!(name, ""); push!(istip, true)
        isempty(stack) || (istip[stack[end]] = false)
        return length(parent)
    end
    readlabel() = begin
        if i <= lastindex(txt) && txt[i] == '\''
            j = i + 1
            buf = IOBuffer()
            while j <= lastindex(txt)
                if txt[j] == '\''
                    if j < lastindex(txt) && txt[nextind(txt, j)] == '\''
                        print(buf, '\''); j = nextind(txt, nextind(txt, j)); continue
                    end
                    break
                end
                print(buf, txt[j]); j = nextind(txt, j)
            end
            i = nextind(txt, j)
            String(take!(buf))
        else
            j = i
            while j <= lastindex(txt) && !(txt[j] in "(),:;")
                j = nextind(txt, j)
            end
            lab = strip(txt[i:prevind(txt, j)])
            i = j
            String(lab)
        end
    end
    readlength() = begin
        if i <= lastindex(txt) && txt[i] == ':'
            j = nextind(txt, i)
            k = j
            while k <= lastindex(txt) && !(txt[k] in "(),;")
                k = nextind(txt, k)
            end
            v = parse(Float64, strip(txt[j:prevind(txt, k)]))
            i = k
            v
        else
            0.0
        end
    end
    while i <= lastindex(txt)
        ch = txt[i]
        if ch == '('
            cur = newnode()
            push!(stack, cur)
            i = nextind(txt, i)
        elseif ch == ','
            i = nextind(txt, i)
        elseif ch == ')'
            cur = pop!(stack)
            i = nextind(txt, i)
            name[cur] = readlabel()
            edge[cur] = readlength()
        elseif isspace(ch)
            i = nextind(txt, i)
        else
            cur = newnode()
            name[cur] = readlabel()
            edge[cur] = readlength()
        end
    end
    # renumber: tips first, in reading order
    tips = findall(istip)
    ints = findall(!, istip)
    order = vcat(tips, ints)
    newid = zeros(Int, length(order))
    newid[order] = 1:length(order)
    p2 = [parent[v] == 0 ? 0 : newid[parent[v]] for v in order]
    nm = [istip[v] ? name[v] : "" for v in order]
    return _tree(p2, edge[order], nm, length(tips))
end

"""
    faith_pd(t::PhyloTree, present)

Faith's phylogenetic diversity: total length of the branches joining the
root to the tips in `present` (a `Bool` vector over the tips, or tip names).
"""
function faith_pd(t::PhyloTree, present::AbstractVector{Bool})
    length(present) == t.ntips || throw(DimensionMismatch("need one Bool per tip"))
    below = falses(nnodes(t))
    below[1:t.ntips] .= present
    pd = 0.0
    for v in postorder(t)
        v == t.root && continue
        if below[v]
            pd += t.edge[v]
            below[t.parent[v]] = true
        end
    end
    return pd
end
faith_pd(t::PhyloTree, names::AbstractVector{<:AbstractString}) = faith_pd(t, in(Set(names)).(tiplabels(t)))

# distance of every node from the root
function _depths(t::PhyloTree)
    dep = zeros(nnodes(t))
    for v in reverse(postorder(t))
        v == t.root || (dep[v] = dep[t.parent[v]] + t.edge[v])
    end
    return dep
end

# tip-to-tip (patristic) distances
function cophenetic(t::PhyloTree)
    adj = _adjacency(t)
    nt = t.ntips
    M = zeros(nt, nt)
    for i in 1:nt
        d, _ = _distances_from(adj, i)
        M[i, :] = d[1:nt]
    end
    return M
end

Base.show(io::IO, t::PhyloTree) = print(io, "PhyloTree($(t.ntips) tips)")
function Base.show(io::IO, ::MIME"text/plain", t::PhyloTree)
    _title(io, "PhyloTree", "$(_fmt(t.ntips)) tips, rooted")
    dep = _depths(t)
    _kv(io, "depth", "max root-to-tip " * _g(maximum(dep[1:t.ntips])))
    _kv(io, "length", "total branch length " * _g(sum(t.edge)))
    labs = tiplabels(t)
    k = min(length(labs), 4)
    print(io, _key(io, "tips"), join(labs[1:k], ", "), length(labs) > k ? ", …" : "")
end
