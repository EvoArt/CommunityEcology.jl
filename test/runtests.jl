using Test, Random, DataFrames, CommunityEcology
using CommunityEcology: nsamples, ntips, PhyloTree, DistanceMatrix
using CairoMakie, CSV
import BIOMFormat, Microbiome

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
    figures = [plot_ordination(o; colorby=:group),
               plot_rda(rda(C, [:group, :temperature])), plot_composition(C),
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

@testset "Example data and table constructors" begin
    d = dune()
    @test size(d) == (20, 30)
    @test sum(counts(d)) == 685
    @test metadata(d).Management[1:2] == ["SF", "BF"]
    @test taxonomy(d).Family[1:2] == ["Asteraceae", "Poaceae"]
    @test size(varespec()) == (24, 44)
    @test "pH" in names(metadata(varespec()))

    wide = DataFrame(sample=IDS, F1=X[:, 1], F2=X[:, 2], F3=X[:, 3], F4=X[:, 4])
    c = Community(wide; metadata=META)
    @test counts(c) == X && features(c) == FEATURES && metadata(c).group == META.group
    @test counts(Community(wide; features=[:F2, :F1])) == X[:, [2, 1]]
    @test_throws ArgumentError Community(wide; sample=:site)
    mktempdir() do dir
        f = joinpath(dir, "counts.csv")
        CSV.write(f, wide)
        @test counts(Community(f)) == X
    end
end

@testset "BIOM tables" begin
    seqs = [repeat("ACGT", 6), repeat("TTGCA", 5), repeat("GGCAT", 5), repeat("CATG", 6)]
    lineage = [["k__Bacteria", "p__Firmicutes"], ["k__Bacteria", "p__"], String[], ["k__Archaea"]]
    t = BIOMFormat.BIOMTable(permutedims(X), seqs, IDS;
                             observation_metadata=[Dict("taxonomy" => l) for l in lineage],
                             sample_metadata=[Dict("group" => g) for g in META.group])
    c = Community(t)
    @test counts(c) == X
    @test features(c) == ["ASV1", "ASV2", "ASV3", "ASV4"] && sequences(c) == seqs
    @test isequal(taxonomy(c).Phylum, ["Firmicutes", missing, missing, missing])
    @test metadata(c).group == META.group
    top = Community(t; top=2, metadata=META)
    @test sequences(top) == seqs[[1, 3]] && counts(top) == X[:, [1, 3]]
end

@testset "cophenetic distances follow features" begin
    t = read_newick("((F1:1,F2:2):1,(F3:1,F4:1):1);")
    c = Community(X; samples=IDS, features=FEATURES, tree=t)
    @test cophenetic(c)[1, 2] == 3 && cophenetic(c)[1, 3] == 4
    @test cophenetic(c[:, ["F3", "F1"]]) == [0 4; 4 0]
end

@testset "Microbiome.jl round trip" begin
    p = Microbiome.CommunityProfile(C)
    @test Microbiome.abundances(p) == permutedims(X)
    @test Microbiome.samplenames(p) == IDS && Microbiome.featurenames(p) == FEATURES
    back = Community(p)
    @test counts(back) == X && metadata(back).group == META.group
end

@testset "add_taxonomy and add_tree" begin
    tax = DataFrame(Phylum=["A", "A", "B", "B"])
    @test taxonomy(add_taxonomy(C, tax)).Phylum == tax.Phylum
    t = add_tree(C, "((F1:1,F2:2):1,(F3:1,(F4:1,F5:1):1):1);")
    @test ntips(tree(t)) == 4
    @test_throws ArgumentError add_tree(C, "((F1:1,F2:2):1,F3:1);")
end

@testset "no export clashes with DataFrames" begin
    @test (@__MODULE__).metadata === DataFrames.metadata
    @test metadata(C).group == META.group
end

@testset "short names read from CSV (InlineStrings)" begin
    mktempdir() do dir
        write_community(dir, C)
        c = read_community(dir)
        @test samples(c) == IDS && eltype(samples(c)) == String
    end
end

@testset "metadata with its own missing values" begin
    meta = DataFrame(sample=IDS, ph=[7.1, missing, 6.8, 7.0, 6.5, missing])
    c = Community(X; samples=IDS, features=FEATURES, metadata=meta)
    @test isequal(metadata(c).ph, meta.ph)
end

@testset "beta_partition (Baselga multi-site)" begin
    # three sites {a,b,c}, {a,b}, {d}: Σ Sᵢ - S_T = 2, Σ min = 2, Σ max = 6
    P = [1 1 1 0; 1 1 0 0; 0 0 0 1]
    c = Community(P; samples=["s1", "s2", "s3"], features=["a", "b", "c", "d"],
                  metadata=DataFrame(sample=["s1", "s2", "s3"], g=["x", "x", "y"]))
    j = beta_partition(c)
    @test j.total ≈ 0.8 && j.turnover ≈ 2 / 3 && j.nestedness ≈ 0.8 - 2 / 3 && j.relative_turnover ≈ 5 / 6
    s = beta_partition(c; family=:sorensen)
    @test s.total ≈ 2 / 3 && s.turnover ≈ 0.5 && s.nestedness ≈ 1 / 6
    @test_throws ArgumentError beta_partition(c; family=:bray)

    # agrees with the pairwise-loop implementation of DiversityFunctions.jl (Jaccard family)
    function βJAC_loop(x)
        b(x₁, x₂) = sum(x₁ .& .!x₂)
        num = 0; denom = sum(x[1, :]) - sum(sum(x; dims=1) .> 0)
        for i in 2:size(x, 1)
            denom += sum(x[i, :])
            for j in 1:i
                num += b(x[i, :], x[j, :]) + b(x[j, :], x[i, :])
            end
        end
        num / (num + denom)
    end
    Y = rand(MersenneTwister(1), Bool, 9, 12)
    r = Community(Int.(Y); samples=string.("s", 1:9), features=string.("f", 1:12))
    @test beta_partition(r).total ≈ βJAC_loop(Y)

    t = beta_partition(c, :g; dispersion=nothing)
    @test DataFrame(t).g == ["x", "y"] && size(t) == (2, 5)
    @test isnan(DataFrame(t).βJAC[2])                     # one sample: undefined
    @test_throws ArgumentError beta_partition(c[1:1, :])
    @test occursin("| g | βJAC", sprint(show, MIME"text/markdown"(), t))
    v = beta_partition(dune(), :Management)
    @test "dispersion" in names(DataFrame(v)) && all(0 .<= DataFrame(v).βJAC .<= 1)
end
