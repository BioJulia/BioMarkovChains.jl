module DiscreteMarkovChainsExt

using BioMarkovChains
using DiscreteMarkovChains:
    DiscreteMarkovChains, DiscreteMarkovChain,
    is_absorbing, is_ergodic, is_regular, is_reversible

## Boolean functions ##

import DiscreteMarkovChains: is_absorbing, is_ergodic, is_regular, is_reversible

function DiscreteMarkovChains.is_absorbing(bmc::BioMarkovChain)
    return is_absorbing(DiscreteMarkovChain(bmc.tpm))
end

function DiscreteMarkovChains.is_ergodic(bmc::BioMarkovChain)
    return is_ergodic(DiscreteMarkovChain(bmc.tpm))
end

function DiscreteMarkovChains.is_regular(bmc::BioMarkovChain)
    return is_regular(DiscreteMarkovChain(bmc.tpm))
end

function DiscreteMarkovChains.is_reversible(bmc::BioMarkovChain)
    return is_reversible(DiscreteMarkovChain(bmc.tpm))
end

## 

end