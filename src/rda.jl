"""
    RDAResult

Constrained redundancy analysis from [`rda`](@ref). `sites` and `features`
are score matrices (rows are samples and features respectively), `eigenvalues`
are constrained eigenvalues, `explained` gives each axis' share of total
response inertia, and `r2`/`adjusted_r2` summarize the full fitted model.
"""
struct RDAResult
    sites::Matrix{Float64}
    feature_scores::Matrix{Float64}
    samples::Vector{String}
    features::Vector{String}
    predictors::Vector{Symbol}
    eigenvalues::Vector{Float64}
    explained::Vector{Float64}
    r2::Float64
    adjusted_r2::Float64
    transform::Symbol
end

"""
    rda(c::Community, predictors; transform=:hellinger, k=2)

Redundancy analysis of a community table constrained by numeric or categorical
sample metadata columns. The response is Hellinger transformed by default;
`:relative` and `:raw` are also available. Numeric predictors are centred;
categorical predictors use reference-level dummy coding. The result contains
fitted sample coordinates (`sites`), feature loadings (`feature_scores`),
constrained eigenvalues, and adjusted R². Permutation inference for a model
can be performed separately with [`permanova`](@ref).

This is a basic RDA with no conditioning terms, weights, or formula interface.
"""
function rda(c::Community, predictors::AbstractVector{Symbol}; transform::Symbol=:hellinger,
             k::Integer=2)
    n, p = size(c)
    n >= 3 && p >= 1 || throw(ArgumentError("RDA needs at least 3 samples and 1 feature"))
    k >= 1 || throw(ArgumentError("k must be positive"))
    isempty(predictors) && throw(ArgumentError("supply at least one metadata predictor"))
    allunique(predictors) || throw(ArgumentError("predictors must be unique"))
    transform in (:hellinger, :relative, :raw) ||
        throw(ArgumentError("transform must be :hellinger, :relative, or :raw"))
    Y = transform === :hellinger ? hellinger(c).counts :
        transform === :relative ? relative(c).counts : Float64.(c.counts)
    Y = Matrix{Float64}(Y)
    Y .-= mean(Y; dims=1)
    total = sum(abs2, Y)
    total > 0 || throw(ArgumentError("community responses have no variation"))
    cols = Vector{Float64}[]
    for predictor in predictors
        predictor in propertynames(c.metadata) ||
            throw(ArgumentError("metadata has no column $predictor"))
        col = c.metadata[!, predictor]
        any(ismissing, col) && throw(ArgumentError("predictor $predictor has missing values"))
        if all(v -> v isa Real && !(v isa Bool), col)
            v = Float64.(col)
            all(isfinite, v) || throw(ArgumentError("predictor $predictor has non-finite values"))
            v .-= mean(v)
            iszero(norm(v)) || push!(cols, v)
        else
            levels = sort!(unique(string.(col)))
            for level in levels[2:end]
                v = Float64.(string.(col) .== level)
                v .-= mean(v)
                push!(cols, v)
            end
        end
    end
    isempty(cols) && throw(ArgumentError("predictors have no variation"))
    X = hcat(cols...)
    q = size(X, 2)
    q < n - 1 || throw(ArgumentError("too many predictor columns for $n samples"))
    rank(X) == q || throw(ArgumentError("predictors are linearly dependent"))
    fitted = X * (X \ Y)
    decomposition = svd(fitted; full=false)
    singular = decomposition.S
    tol = sqrt(eps(Float64)) * maximum(singular)
    available = count(>(tol), singular)
    available > 0 || throw(ArgumentError("predictors explain no response variation"))
    axes = min(k, available)
    sites = Matrix(decomposition.U[:, 1:axes] * Diagonal(singular[1:axes]))
    feature_scores = Matrix(decomposition.V[:, 1:axes])
    for j in 1:axes
        i = argmax(abs.(view(sites, :, j)))
        if sites[i, j] < 0
            sites[:, j] .*= -1
            feature_scores[:, j] .*= -1
        end
    end
    eigenvalues = singular[1:axes] .^ 2 ./ (n - 1)
    explained = singular[1:axes] .^ 2 ./ total
    r2 = clamp(sum(abs2, fitted) / total, 0.0, 1.0)
    adjusted_r2 = 1 - (1 - r2) * (n - 1) / (n - q - 1)
    return RDAResult(sites, feature_scores, copy(c.samples), copy(c.features),
                     collect(predictors), eigenvalues, explained, r2,
                     adjusted_r2, transform)
end

rda(c::Community, predictor::Symbol; kw...) = rda(c, [predictor]; kw...)
