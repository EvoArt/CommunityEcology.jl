# Small classic data sets for examples and teaching, as in R's vegan.

const _DATA = joinpath(@__DIR__, "..", "data")

_readdata(file) = CSV.read(joinpath(_DATA, file), DataFrame)

function _community(file, meta; taxonomy=nothing)
    df = _readdata(file)
    return Community(df; metadata=meta, taxonomy=taxonomy)
end

"""
    dune()

Vegetation of 20 dune meadows on Terschelling, the Netherlands (Jongman, ter Braak
& van Tongeren 1987): cover class (0–9) of 30 plant species, with metadata
`A1` (thickness of the A1 soil horizon, cm), `Moisture` (1, 2, 4 or 5),
`Management` (`"BF"` biological farming, `"HF"` hobby farming, `"NM"` nature
conservation, `"SF"` standard farming), `Use` (`"Hayfield"`, `"Haypastu"`,
`"Pasture"`) and `Manure` (0–4). The taxonomy gives each species' genus,
family, order, superorder and subclass.

The same data as `dune`, `dune.env` and `dune.taxon` in R's vegan package.

```julia
c = dune()
permanova(c, @formula(1 ~ Management))
```
"""
dune() = _community("dune.csv", _readdata("dune_env.csv");
                    taxonomy=select(_readdata("dune_taxon.csv"), Not(:species)))

"""
    varespec()

Cover of 44 lichen, moss and vascular plant species in 24 lichen pastures
under Scots pine in eastern Fennoscandia (Väre, Ohtonen & Oksanen 1995),
with soil chemistry as metadata: `N`, `P`, `K`, `Ca`, `Mg`, `S`, `Al`,
`Fe`, `Mn`, `Zn`, `Mo` (mg/kg), `Baresoil` (% cover), `Humdepth` (humus
depth, cm) and `pH`.

The same data as `varespec` and `varechem` in R's vegan package.
"""
varespec() = _community("varespec.csv", _readdata("varechem.csv"))

const _MISEQ_URL = "https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"

"""
    miseq_sop_reads(dir="MiSeq_SOP")

Paired-end 16S rRNA (V4) MiSeq reads from the mothur MiSeq SOP (Kozich et al.
2013): faecal samples from one mouse on days 0–9 and 141–150 after weaning,
plus a mock community of 20 known strains (`"Mock"`). The 35 MB archive is
downloaded and unpacked into `dir` the first time; later calls only list the
files.

Returns a `DataFrame` with one row per sample: `sample`, and the paths of the
`forward` and `reverse` FASTQ files.

```julia
reads = miseq_sop_reads()
filterAndTrim(reads.forward, filtF, reads.reverse, filtR; truncLen=(240, 160))
```
"""
function miseq_sop_reads(dir::AbstractString="MiSeq_SOP")
    if !isdir(dir) || !any(endswith("_R1_001.fastq"), readdir(dir))
        mktempdir() do tmp
            zipfile = Downloads.download(_MISEQ_URL, joinpath(tmp, "miseqsopdata.zip"))
            run(pipeline(`$(p7zip_jll.p7zip()) x -y -o$tmp $zipfile`; stdout=devnull))
            mkpath(dirname(abspath(dir)))
            mv(joinpath(tmp, "MiSeq_SOP"), dir; force=true)
        end
    end
    R1 = sort(filter(endswith("_R1_001.fastq"), readdir(dir; join=true)))
    R2 = sort(filter(endswith("_R2_001.fastq"), readdir(dir; join=true)))
    return DataFrame(sample=first.(split.(basename.(R1), "_")), forward=R1, reverse=R2)
end

"""
    miseq_sop()

The mouse gut samples of the mothur MiSeq SOP ([`miseq_sop_reads`](@ref)),
processed into a community: 19 samples × the amplicon sequence variants
(ASVs) that DenoiseDNA infers from the reads, with their sequences, SILVA
138.2 taxonomy (AssignTaxonomy), a neighbour-joining tree (MSABuilder and
PhyloBuilder), and metadata `day` (days post weaning) and `time`
(`"Early"`, days 0–9; `"Late"`, days 141–150). The mock community is left out.

`data/miseq_sop/build.jl` in this package is the script that made it.
"""
miseq_sop() = read_community(joinpath(_DATA, "miseq_sop"))
