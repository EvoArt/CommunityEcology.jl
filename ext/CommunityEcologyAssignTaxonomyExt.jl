module CommunityEcologyAssignTaxonomyExt

using CommunityEcology
import AssignTaxonomy

function CommunityEcology._classify(c::Community, reference; kw...)
    isempty(sequences(c)) && throw(ArgumentError("the community has no sequences to classify"))
    return add_taxonomy(c, AssignTaxonomy.assignTaxonomy(sequences(c), reference; kw...))
end

end
