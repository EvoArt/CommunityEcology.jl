# Pairwise distances between amplicon sequences, for building a tree.
#
# Each pair is aligned globally with free end gaps (match 1, mismatch -1,
# gap -2) in a band around the diagonal: ASVs of one amplicon differ by
# few indels, so a band of ±32 columns (plus the length difference) holds
# the best alignment and keeps each pair to a few tens of microseconds.

const _MATCH, _MISMATCH, _GAP = Int32(1), Int32(-1), Int32(-2)

struct _AlignBuffer
    H::Matrix{Int32}
    T::Matrix{UInt8}
end
_AlignBuffer(m, n) = _AlignBuffer(Matrix{Int32}(undef, m + 1, n + 1), Matrix{UInt8}(undef, m + 1, n + 1))

# (differences, aligned columns) between a and b, end gaps not counted
function _align_diffs(a, b, band, buf)
    m, n = length(a), length(b)
    H, T = buf.H, buf.T
    NEG = typemin(Int32) ÷ Int32(4)
    kmin = min(0, n - m) - band
    kmax = max(0, n - m) + band
    @inbounds begin
        for j in 0:n
            H[1, j+1] = 0                     # free leading gaps
        end
        for i in 1:m
            H[i+1, 1] = 0
            jlo = max(1, i + kmin)
            jhi = min(n, i + kmax)
            jlo > 1 && (H[i+1, jlo] = NEG)    # left edge of the band
            for j in jlo:jhi
                up = j - i + 1 <= kmax ? H[i, j+1] : NEG   # above lies outside the band?
                d = H[i, j] + (a[i] == b[j] ? _MATCH : _MISMATCH)
                u = up + _GAP
                l = H[i+1, j] + _GAP
                if d >= u && d >= l
                    H[i+1, j+1] = d; T[i+1, j+1] = 0x01
                elseif u >= l
                    H[i+1, j+1] = u; T[i+1, j+1] = 0x02
                else
                    H[i+1, j+1] = l; T[i+1, j+1] = 0x03
                end
            end
            jhi < n && (H[i+1, jhi+2] = NEG)  # right edge of the band
        end
        # free trailing gaps: best cell in the last row or column (in band)
        bi, bj, best = m, min(n, m + kmax), NEG
        for j in max(1, m + kmin):min(n, m + kmax)
            H[m+1, j+1] > best && ((bi, bj, best) = (m, j, H[m+1, j+1]))
        end
        for i in max(1, n - kmax):m
            H[i+1, n+1] > best && ((bi, bj, best) = (i, n, H[i+1, n+1]))
        end
        diffs = 0
        cols = 0
        i, j = bi, bj
        while i > 0 && j > 0
            t = T[i+1, j+1]
            cols += 1
            if t == 0x01
                diffs += a[i] != b[j]
                i -= 1; j -= 1
            elseif t == 0x02
                diffs += 1; i -= 1
            else
                diffs += 1; j -= 1
            end
        end
    end
    return diffs, cols
end

"""
    sequence_distances(seqs; band=32, model=:jc69)

Pairwise evolutionary distances between DNA sequences (strings or anything
that converts to one), from banded global alignments with free end gaps.
`model=:jc69` applies the Jukes-Cantor correction to the proportion of
differing aligned sites (gaps count as differences); `:p` returns the raw
proportion. Runs on all Julia threads.
"""
function sequence_distances(seqs; band::Integer=32, model::Symbol=:jc69)
    model in (:jc69, :p) || throw(ArgumentError("model must be :jc69 or :p"))
    s = [Vector{UInt8}(uppercase(string(x))) for x in seqs]
    n = length(s)
    L = maximum(length, s; init=0)
    D = zeros(n, n)
    pairs = [(i, j) for i in 1:n for j in i+1:n]
    nchunk = max(1, min(Threads.nthreads(), length(pairs)))
    Threads.@threads for c in 1:nchunk
        buf = _AlignBuffer(L, L)
        for k in c:nchunk:length(pairs)
            i, j = pairs[k]
            d, cols = _align_diffs(s[i], s[j], band, buf)
            p = cols == 0 ? 1.0 : d / cols
            D[i, j] = D[j, i] = model === :p ? p : _jc69(p)
        end
    end
    return D
end

# Jukes-Cantor (1969); saturated pairs get a large finite distance
_jc69(p) = p < 0.7499 ? -0.75 * log(1 - 4p / 3) : 6.0

"""
    build_tree(c::Community; band=32)
    build_tree(seqs; labels, band=32)

A neighbour-joining tree (rooted at its midpoint) of the feature sequences
of `c`, from Jukes-Cantor distances of pairwise alignments; returns a copy
of `c` with the tree attached, ready for UniFrac and Faith's PD. Given
sequences, returns the [`PhyloTree`](@ref).

This is a quick distance tree, good enough for UniFrac and for a picture.
For a phylogeny to interpret on its own, align with a multiple aligner and
build a maximum-likelihood tree; `read_newick` reads the result.
"""
function build_tree(c::Community; band::Integer=32)
    isempty(c.sequences) && throw(ArgumentError("the community has no sequences to build a tree from"))
    t = build_tree(c.sequences; labels=c.features, band=band)
    return Community(c.counts, c.samples, c.features, c.sequences, c.taxonomy, c.metadata, t)
end
build_tree(seqs::AbstractVector; labels=string.(1:length(seqs)), band::Integer=32) =
    nj_tree(sequence_distances(seqs; band=band); labels=labels)
