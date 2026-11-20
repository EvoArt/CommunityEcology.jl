"""
    MantelResult

Result of [`mantel`](@ref): correlation `r`, permutation `pvalue`, number of
`permutations`, correlation `method`, and test `alternative`.
"""
struct MantelResult
    r::Float64
    pvalue::Float64
    permutations::Int
    method::Symbol
    alternative::Symbol
end

"""
    mantel(a, b; method=:pearson, alternative=:greater, permutations=999,
           rng=Random.default_rng())

Mantel test of association between two distance matrices with the same sample
order. `method=:spearman` correlates average ranks (including ties). The
permutation test relabels the rows and columns of `b` together. The p-value
includes the observed arrangement, so its minimum is `1/(permutations+1)`.
Use `alternative=:less` or `:two_sided` for other tails.
"""
function mantel(a::AbstractMatrix, b::AbstractMatrix; method::Symbol=:pearson,
                alternative::Symbol=:greater, permutations::Integer=999,
                rng::AbstractRNG=Random.default_rng())
    method in (:pearson, :spearman) || throw(ArgumentError("method must be :pearson or :spearman"))
    alternative in (:greater, :less, :two_sided) ||
        throw(ArgumentError("alternative must be :greater, :less or :two_sided"))
    permutations >= 0 || throw(ArgumentError("permutations must be non-negative"))
    n = size(a, 1)
    size(a) == (n, n) && size(b) == (n, n) && n >= 3 ||
        throw(DimensionMismatch("both distance matrices must be square, equally sized, and have at least 3 samples"))
    _validate_distance(a)
    _validate_distance(b)
    x = _lower_triangle(a)
    y = _lower_triangle(b)
    method === :spearman && (x = _average_ranks(x); y = _average_ranks(y))
    r = cor(x, y)
    isfinite(r) || throw(ArgumentError("distance vectors must have nonzero variance"))
    exceed = 0
    for _ in 1:permutations
        p = randperm(rng, n)
        yp = [b[p[i], p[j]] for j in 1:n for i in j+1:n]
        method === :spearman && (yp = _average_ranks(yp))
        rp = cor(x, yp)
        exceed += alternative === :greater ? rp >= r - 1e-12 :
                  alternative === :less ? rp <= r + 1e-12 : abs(rp) >= abs(r) - 1e-12
    end
    return MantelResult(r, (exceed + 1) / (permutations + 1), Int(permutations), method, alternative)
end

function mantel(a::DistanceMatrix, b::DistanceMatrix; kw...)
    a.labels == b.labels || throw(ArgumentError("distance matrices have different sample labels or order"))
    return mantel(a.D, b.D; kw...)
end

_lower_triangle(D) = Float64[D[i, j] for j in 1:size(D, 1) for i in j+1:size(D, 1)]

function _validate_distance(D)
    all(isfinite, D) && all(>=(0), D) ||
        throw(ArgumentError("distances must be finite and non-negative"))
    for i in axes(D, 1)
        isapprox(D[i, i], 0; atol=1e-10) || throw(ArgumentError("distance diagonal must be zero"))
        for j in i+1:size(D, 1)
            isapprox(D[i, j], D[j, i]; rtol=1e-8, atol=1e-10) ||
                throw(ArgumentError("distance matrix must be symmetric"))
        end
    end
end

function _average_ranks(x)
    order = sortperm(x)
    ranks = Vector{Float64}(undef, length(x))
    i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[order[j+1]] == x[order[i]]
            j += 1
        end
        ranks[order[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return ranks
end

"""
    envfit(o::Ordination; variables=nothing, permutations=999,
           rng=Random.default_rng())

Fit sample metadata to an ordination, in the style of `vegan::envfit`.
Numeric variables use vector fitting; categorical variables use factor
centroids. Returns a `DataFrame` with `variable`, `level` (`missing` for numeric
vectors), squared correlation `r2`, permutation `pvalue`, and `axis1`, `axis2`,
etc. Numeric coordinates are unit vector directions; factor coordinates are
level centroids. Missing values are omitted per variable.
`variables=nothing` uses all metadata columns other than `sample`.
"""
function envfit(o::Ordination; variables=nothing, permutations::Integer=999,
                rng::AbstractRNG=Random.default_rng())
    permutations >= 0 || throw(ArgumentError("permutations must be non-negative"))
    vars = variables === nothing ?
        [Symbol(n) for n in names(o.metadata) if n != "sample"] : Symbol.(collect(variables))
    isempty(vars) && throw(ArgumentError("no metadata variables to fit"))
    k = size(o.points, 2)
    out = DataFrame(variable=String[], level=Union{Missing,String}[],
                    r2=Float64[], pvalue=Float64[])
    for j in 1:k
        out[!, Symbol("axis$j")] = Float64[]
    end
    for variable in vars
        variable in propertynames(o.metadata) || throw(ArgumentError("metadata has no column $variable"))
        variable === :sample && throw(ArgumentError("sample names are identifiers, not environmental variables"))
        col = o.metadata[!, variable]
        if !all(v -> ismissing(v) || (v isa Real && !(v isa Bool)), col)
            append!(out, factorfit(o, variable; permutations=permutations, rng=rng))
            continue
        end
        keep = findall(!ismissing, col)
        length(keep) > k + 1 || throw(ArgumentError("$variable needs more than $(k + 1) non-missing samples"))
        X = o.points[keep, :]
        X = X .- mean(X; dims=1)
        rank(X) == k || throw(ArgumentError("ordination axes are linearly dependent for $variable"))
        y = Float64.(col[keep])
        y .-= mean(y)
        sum(abs2, y) > 0 || throw(ArgumentError("metadata column $variable is constant"))
        beta = X \ y
        fitted = X * beta
        r2 = clamp(sum(abs2, fitted) / sum(abs2, y), 0.0, 1.0)
        direction = iszero(norm(beta)) ? zeros(k) : beta / norm(beta)
        exceed = 0
        for _ in 1:permutations
            yp = y[randperm(rng, length(y))]
            fp = X * (X \ yp)
            exceed += sum(abs2, fp) / sum(abs2, yp) >= r2 - 1e-12
        end
        push!(out, (string(variable), missing, r2,
                    (exceed + 1) / (permutations + 1), direction...))
    end
    return out
end

"""
    factorfit(o::Ordination, variable::Symbol; permutations=999,
              rng=Random.default_rng())

Fit a categorical metadata variable to an ordination. Returns level centroids
with a shared between-group R² and permutation p-value. The result has the
same columns as [`envfit`](@ref), so results can be combined.
"""
function factorfit(o::Ordination, variable::Symbol; permutations::Integer=999,
                   rng::AbstractRNG=Random.default_rng())
    permutations >= 0 || throw(ArgumentError("permutations must be non-negative"))
    variable in propertynames(o.metadata) || throw(ArgumentError("metadata has no column $variable"))
    col = o.metadata[!, variable]
    all(v -> ismissing(v) || v isa Real && !(v isa Bool), col) &&
        throw(ArgumentError("metadata column $variable is numeric; use envfit for vector fitting"))
    keep = findall(!ismissing, col)
    groups = string.(col[keep])
    levels = sort!(unique(groups))
    length(levels) >= 2 || throw(ArgumentError("$variable needs at least two levels"))
    length(keep) > length(levels) ||
        throw(ArgumentError("$variable needs more samples than levels"))
    X = o.points[keep, :]
    r2, centroids = _factor_stats(X, groups, levels)
    exceed = 0
    for _ in 1:permutations
        permuted = groups[randperm(rng, length(groups))]
        rp, _ = _factor_stats(X, permuted, levels)
        exceed += rp >= r2 - 1e-12
    end
    k = size(X, 2)
    out = DataFrame(variable=String[], level=Union{Missing,String}[],
                    r2=Float64[], pvalue=Float64[])
    for j in 1:k
        out[!, Symbol("axis$j")] = Float64[]
    end
    pvalue = (exceed + 1) / (permutations + 1)
    for (level, centroid) in zip(levels, centroids)
        push!(out, (string(variable), level, r2, pvalue, centroid...))
    end
    return out
end

function _factor_stats(X, groups, levels)
    grand = vec(mean(X; dims=1))
    total = sum(abs2, X .- grand')
    total > 0 || throw(ArgumentError("ordination points have no variation"))
    centroids = [vec(mean(X[findall(==(level), groups), :]; dims=1)) for level in levels]
    between = sum(count(==(level), groups) * sum(abs2, centroid .- grand)
                  for (level, centroid) in zip(levels, centroids))
    return clamp(between / total, 0.0, 1.0), centroids
end
