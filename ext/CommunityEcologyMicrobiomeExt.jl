module CommunityEcologyMicrobiomeExt

# Conversions between CommunityEcology's Community (samples × features) and
# Microbiome.jl's CommunityProfile (features × samples), so tables read by
# Microbiome.jl or BiobakeryUtils.jl (e.g. MetaPhlAn profiles) can use the
# analyses here, and the other way round.

using CommunityEcology, DataFrames
import Microbiome
using Microbiome: CommunityProfile, MicrobiomeSample, Taxon

"""
    Community(p::Microbiome.CommunityProfile; kw...)

A community from a Microbiome.jl profile: its abundances transposed to
samples × features, with the profile's sample metadata.
"""
function CommunityEcology.Community(p::CommunityProfile; metadata=nothing, kw...)
    meta = metadata === nothing ? DataFrame(get(p)) : metadata
    return Community(Matrix(permutedims(Microbiome.abundances(p)));
                     samples=Microbiome.samplenames(p), features=Microbiome.featurenames(p),
                     metadata=meta, kw...)
end

"""
    CommunityProfile(c::Community)

A Microbiome.jl `CommunityProfile` (features × samples) with the community's
sample metadata. Features become `Taxon`s without a rank.
"""
function Microbiome.CommunityProfile(c::Community)
    meta = select(metadata(c), Not(:sample))
    smpls = [MicrobiomeSample(s, NamedTuple(row)) for (s, row) in zip(samples(c), eachrow(meta))]
    return CommunityProfile(permutedims(counts(c)), Taxon.(features(c)), smpls)
end

end
