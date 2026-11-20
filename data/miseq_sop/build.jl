# Builds data/miseq_sop (what `miseq_sop()` loads) from the raw MiSeq SOP reads,
# following the "Reads to ASVs", "Taxonomy" and "Phylogeny" lessons.
#
#   julia -t auto --project=<env> build.jl <silva_nr99_v138.2_toGenus_trainset.fa.gz> [workdir]
#
# <env> needs CommunityEcology, DenoiseDNA, AssignTaxonomy, MSABuilder, PhyloBuilder and
# DataFrames. The SILVA 138.2 training set is at https://zenodo.org/records/14169026.
using CommunityEcology, DenoiseDNA, AssignTaxonomy, DataFrames
using MSABuilder: align
using PhyloBuilder: nj
using PhyloBuilder.PhylogenyIO: writenewick

silva = ARGS[1]
work = length(ARGS) > 1 ? ARGS[2] : mktempdir()

reads = miseq_sop_reads(joinpath(work, "MiSeq_SOP"))
filtF = joinpath.(work, "filtered", reads.sample .* "_F.fastq.gz")
filtR = joinpath.(work, "filtered", reads.sample .* "_R.fastq.gz")
filterAndTrim(reads.forward, filtF, reads.reverse, filtR;
              truncLen=(240, 160), maxEE=(2, 2), rm_phix=true, multithread=true)
errF = learnErrors(filtF; multithread=true)
errR = learnErrors(filtR; multithread=true)
merged = mergePairs(dada(filtF, errF; multithread=true), filtF,
                    dada(filtR, errR; multithread=true), filtR)
seqtab = removeBimeraDenovo(makeSequenceTable(merged; names=reads.sample);
                            method="consensus", multithread=true)

taxa = assignTaxonomy(seqtab, silva; seed=1)
asvs = Community(seqtab; taxonomy=taxa)

mice = prune(filter_samples(asvs, :sample => !=("Mock")))
day = parse.(Int, replace.(samples(mice), "F3D" => ""))
mice = add_metadata(mice, DataFrame(sample=samples(mice), day=day,
                                    time=ifelse.(day .< 100, "Early", "Late")))

write_community(@__DIR__, mice)
writenewick(joinpath(@__DIR__, "tree.nwk"),
            nj(align(features(mice) .=> sequences(mice)); root=:midpoint))
println(miseq_sop())
