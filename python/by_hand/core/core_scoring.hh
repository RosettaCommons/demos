// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef core_scoring
#define core_scoring

/// Why aren't these #include <...>'s?  Aren't we avoiding #include ""'s?

#include "all.fwd.hh"

#include "utility/pointer/owning_ptr.hh"

#include "core/types.hh"

#include "core/init.hh"

#include "core/scoring/types.hh"
#include "core/scoring/hbonds/types.hh"
#include "core/scoring/CachedData.hh"
#include "core/scoring/constants.hh"
#include "core/scoring/ContextGraph.hh"
#include "core/scoring/ContextGraphFactory.hh"
#include "core/scoring/ContextGraphTypes.hh"
#include "core/scoring/dunbrack/RotamerLibrary.hh"
#include "core/scoring/Energies.hh"
#include "core/scoring/EnergyGraph.hh"
#include "core/scoring/EnergyMap.hh"
#include "core/scoring/etable/atom_pair_energy_inline.hh"
#include "core/scoring/etable/count_pair/CountPair1BC3.hh"
#include "core/scoring/etable/count_pair/CountPair1BC4.hh"
#include "core/scoring/etable/count_pair/CountPairAll.hh"
#include "core/scoring/etable/count_pair/CountPairCrossover3.hh"
#include "core/scoring/etable/count_pair/CountPairCrossover4.hh"
#include "core/scoring/etable/count_pair/CountPairFunction.hh"
#include "core/scoring/etable/Etable.hh"
#include "core/scoring/etable/EtableEnergy.hh"
#include "core/scoring/etable/EtableOptions.hh"
#include "core/scoring/etable/etrie/CountPairData_1_1.hh"
#include "core/scoring/etable/etrie/CountPairData_1_2.hh"
#include "core/scoring/etable/etrie/CountPairData_1_3.hh"
#include "core/scoring/etable/etrie/EtableAtom.hh"
#include "core/scoring/etable/etrie/TrieCountPair1BC3.hh"
#include "core/scoring/etable/etrie/TrieCountPair1BC4.hh"
#include "core/scoring/etable/etrie/TrieCountPairAll.hh"
#include "core/scoring/hbonds/constants.hh"
#include "core/scoring/hbonds/create_poly.hh"
#include "core/scoring/hbonds/hbonds.hh"
#include "core/scoring/hbonds/hbonds_geom.hh"
#include "core/scoring/hbonds/HBondSet.hh"
#include "core/scoring/methods/ChainbreakEnergy.hh"
#include "core/scoring/methods/ContextDependentTwoBodyEnergy.hh"
#include "core/scoring/methods/ContextIndependentOneBodyEnergy.hh"
#include "core/scoring/methods/ContextIndependentTwoBodyEnergy.hh"
#include "core/scoring/methods/DunbrackEnergy.hh"
#include "core/scoring/methods/EnergyMethod.hh"
#include "core/scoring/methods/EnergyMethodOptions.hh"
#include "core/scoring/methods/HBondEnergy.hh"
#include "core/scoring/methods/Methods.hh"
#include "core/scoring/methods/PairEnergy.hh"
#include "core/scoring/methods/RamachandranEnergy.hh"
#include "core/scoring/NeighborList.hh"
#include "core/scoring/PairEPotential.hh"
#include "core/scoring/ProteinTorsion.hh"
#include "core/scoring/Ramachandran.hh"
#include "core/scoring/residue_pair_energy.hh"
#include "core/scoring/sasa.hh"
#include "core/scoring/ScoreFunction.hh"
#include "core/scoring/ScoreFunctionInfo.hh"
#include "core/scoring/ScoreFunctionParameter.hh"
#include "core/scoring/ScoreFunctionVariant.hh"
#include "core/scoring/ScoreType.hh"
#include "core/scoring/ScoringManager.hh"
#include "core/scoring/TenANeighborGraph.hh"
#include "core/scoring/trie/RotamerDescriptor.hh"
#include "core/scoring/trie/RotamerTrie.hh"
#include "core/scoring/trie/RotamerTrieBase.hh"
#include "core/scoring/trie/trie_vs_path.hh"
#include "core/scoring/trie/trie_vs_trie.hh"
#include "core/scoring/trie/TrieCollection.hh"
#include "core/scoring/trie/TrieCountPairBase.hh"

#endif
