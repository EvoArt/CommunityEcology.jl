"""
    simper(c::Community, variable::Symbol; relative=false)

Break down mean between-group Bray–Curtis dissimilarity into contributions
from individual features. Returns a `DataFrame` for every pair of levels in
the metadata column `variable`, ordered by decreasing mean contribution.
`average` is the mean contribution on the Bray–Curtis scale; `fraction` and
`cumulative` divide by the total mean dissimilarity. `mean_a` and `mean_b`
are group mean abundances. `relative=true` first converts each sample to
relative abundance. Missing group labels are omitted.

SIMPER is descriptive: a high contribution can reflect abundant, variable
features even without a consistent group difference. Use [`permanova`](@ref)
for a group-level permutation test.
"""
function simper(c::Community, variable::Symbol; relative::Bool=false)
    variable in propertynames(c.metadata) || throw(ArgumentError("metadata has no column $variable"))
    group = c.metadata[!, variable]
    keep = findall(!ismissing, group)
    length(keep) >= 2 || throw(ArgumentError("SIMPER needs at least two labelled samples"))
    any(==(0), sample_depths(c)[keep]) &&
        throw(ArgumentError("samples with no reads must be removed before SIMPER"))
    labels = string.(group)
    levels = sort!(unique(labels[keep]))
    length(levels) >= 2 || throw(ArgumentError("SIMPER needs at least two groups"))
    X = relative ? relative(c).counts : Float64.(c.counts)
    out = DataFrame(group_a=String[], group_b=String[], feature=String[],
                    average=Float64[], sd=Float64[], ratio=Float64[],
                    fraction=Float64[], cumulative=Float64[],
                    mean_a=Float64[], mean_b=Float64[])
    p = nfeatures(c)
    for ai in 1:length(levels)-1, bi in ai+1:length(levels)
        a, b = levels[ai], levels[bi]
        ia = [i for i in keep if labels[i] == a]
        ib = [i for i in keep if labels[i] == b]
        pair_contributions = Matrix{Float64}(undef, length(ia) * length(ib), p)
        row = 0
        for i in ia, j in ib
            row += 1
            denominator = sum(view(X, i, :)) + sum(view(X, j, :))
            pair_contributions[row, :] .= abs.(view(X, i, :) .- view(X, j, :)) ./ denominator
        end
        average = vec(mean(pair_contributions; dims=1))
        spread = size(pair_contributions, 1) > 1 ?
            vec(std(pair_contributions; dims=1)) : zeros(p)
        mean_a = vec(mean(X[ia, :]; dims=1))
        mean_b = vec(mean(X[ib, :]; dims=1))
        overall = sum(average)
        order = sortperm(average; rev=true)
        cumulative = 0.0
        for f in order
            fraction = overall == 0 ? 0.0 : average[f] / overall
            cumulative += fraction
            ratio = spread[f] == 0 ? NaN : average[f] / spread[f]
            push!(out, (a, b, c.features[f], average[f], spread[f], ratio,
                        fraction, cumulative, mean_a[f], mean_b[f]))
        end
    end
    return out
end
