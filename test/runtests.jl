using Test, Random, DataFrames, CommunityEcology
using CairoMakie

const X = [8 1 0 0;
           7 2 0 1;
           6 1 1 0;
           0 1 7 2;
           1 0 6 3;
           0 2 8 1]
const IDS = ["S$i" for i in 1:6]
const FEATURES = ["F$i" for i in 1:4]
const META = DataFrame(sample=IDS, group=repeat(["A", "B"], inner=3),
                       temperature=[10.0, 12, 15, 20, 22, 25])
const C = Community(X; samples=IDS, features=FEATURES, metadata=META)

@testset "Community and diversity" begin
    @test size(C) == (6, 4)
    @test sample_depths(C) == vec(sum(X; dims=2))
    @test counts(C["S1", "F1"]) == reshape([8], 1, 1)
    @test richness([2, 0, 3]) == 2
    @test shannon([1, 1]) ≈ log(2)
    @test simpson([1, 1]) ≈ 0.5
    @test hill([1, 1], 2) ≈ 2
    @test alpha_diversity(C).group == META.group
    @test specaccum(C).richness[end] == count(>(0), vec(sum(X; dims=1)))
    @test specaccum(C).richness[1] ≈ sum(count(>(0), X[i, :]) for i in 1:6) / 6
    @test occursin("Most abundant features", sprint(show, MIME"text/plain"(), C))
    @test all(isapprox.(vec(sum(relative(C).counts; dims=2)), 1))
    @test all(isapprox.(vec(sum(hellinger(C).counts .^ 2; dims=2)), 1))
    @test all(isapprox.(vec(sum(chord(C).counts .^ 2; dims=2)), 1))
    @test all(isapprox.(vec(sum(wisconsin(C).counts; dims=2)), 1))
    @test size(prune(C; min_samples=2)) == (6, 4)
    @test sum(counts(rarefy(C, 5; rng=MersenneTwister(1))); dims=2) == fill(5, 6, 1)
end

@testset "Makie figures" begin
    o = ordinate(C, :pcoa)
    bd = betadisper(C, :group)
    t = read_newick("((F1:1,F2:1):1,(F3:1,F4:1):1);")
    figures = [plot_ordination(o; colorby=:group), plot_composition(C),
               plot_alpha(C; group=:group), plot_rarefaction(C),
               plot_heatmap(C), plot_dispersion(bd), plot_tree(t)]
    @test all(f -> f isa Figure, figures)
    mktempdir() do dir
        file = joinpath(dir, "ordination.png")
        save(file, first(figures))
        @test filesize(file) > 0
    end
end

@testset "Distance, ordination and inference" begin
    d = beta_diversity(C)
    @test d.labels == IDS
    @test d["S1", "S1"] == 0
    @test d["S1", "S4"] > d["S1", "S2"]
    @test Matrix(d) ≈ Matrix(d)'
    @test all(isfinite, beta_diversity(C, :jaccard))
    @test all(isfinite, beta_diversity(C, :hellinger))
    @test all(isfinite, beta_diversity(C, :canberra))
    @test all(isfinite, beta_diversity(C, :gower))
    @test beta_diversity(C, :ruzicka)[1, 4] > beta_diversity(C, :ruzicka)[1, 2]
    @test beta_diversity(C, :ruzicka)[1, 4] ≈ 17 / 18
    @test beta_diversity(C, :canberra)[1, 4] ≈ 0.75
    @test Matrix(beta_diversity(C, :chord)) ≈ Matrix(beta_diversity(chord(C), :euclidean))
    o = ordinate(d, C, :pcoa)
    @test size(o.points) == (6, 2)
    @test o.samples == IDS
    @test size(ordinate(d, C, :nmds; tries=2, maxit=80, threaded=false).points) == (6, 2)
    fit = envfit(o; variables=[:temperature], permutations=19, rng=MersenneTwister(2))
    @test fit.variable == ["temperature"]
    @test 0 <= fit.r2[1] <= 1
    @test 0 < fit.pvalue[1] <= 1
    @test isapprox(hypot(fit.axis1[1], fit.axis2[1]), 1; atol=1e-8)
    factors = envfit(o; variables=[:group], permutations=19, rng=MersenneTwister(2))
    @test factors.level == ["A", "B"]
    @test all((0 .<= factors.r2) .& (factors.r2 .<= 1))
    @test factorfit(o, :group; permutations=0).r2[1] ≈ factors.r2[1]
    constrained = rda(C, [:group, :temperature])
    @test size(constrained.sites, 1) == 6
    @test length(constrained.explained) <= 2
    @test 0 <= constrained.r2 <= 1
    exact = Community([1 0; 1 0; 0 1; 0 1];
                      metadata=DataFrame(sample=["S1", "S2", "S3", "S4"],
                                         group=["A", "A", "B", "B"]))
    @test rda(exact, :group; transform=:raw).r2 ≈ 1
    @test occursin("RDAResult", sprint(show, MIME"text/plain"(), constrained))
    parts = simper(C, :group)
    @test nrow(parts) == nfeatures(C)
    @test sum(parts.average) ≈ sum(d[i, j] for i in 1:3 for j in 4:6) / 9
    @test parts.cumulative[end] ≈ 1
    m = mantel(d, beta_diversity(C, :jaccard); permutations=19, rng=MersenneTwister(3))
    @test -1 <= m.r <= 1
    @test 0 < m.pvalue <= 1
    @test mantel(d, d; permutations=0).r ≈ 1
    @test mantel(d, d; method=:spearman, permutations=0).r ≈ 1
    @test occursin("MantelResult", sprint(show, MIME"text/plain"(), m))
    @test_throws ArgumentError mantel(d, beta_diversity(C[reverse(IDS), :]))
    p = permanova(C, @formula(1 ~ group); permutations=19, rng=MersenneTwister(4))
    @test p.n == 6
    @test p.p[1] <= 1
    bd = betadisper(C, :group)
    @test length(bd.distances) == 6
end

@testset "Phylogeny and file round trip" begin
    t = read_newick("((F1:1,F2:1):1,(F3:1,F4:1):1);")
    c = Community(X; samples=IDS, features=FEATURES, metadata=META, tree=t)
    @test ntips(tree(c)) == 4
    @test faith_pd(t, ["F1", "F2"]) ≈ 3
    @test beta_diversity(c, :unifrac)[1, 4] > 0
    @test alpha_diversity(c; measures=[:faith_pd]).faith_pd[1] > 0
    mktempdir() do dir
        write_community(dir, c)
        loaded = read_community(dir)
        @test counts(loaded) == X
        @test samples(loaded) == IDS
        @test metadata(loaded).group == META.group
        @test ntips(tree(loaded)) == 4
    end
end
