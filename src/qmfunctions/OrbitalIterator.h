/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2020 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#pragma once

#include "orbital_utils.h"

namespace mrchem {

class OrbitalIterator final {
public:
    OrbitalIterator(OrbitalVector &Phi, bool sym = false);

    BankAccount orbBank;
    bool next(int max_recv = -1);
    bool bank_next(int max_recv = 1);

    int idx(int i) const { return (this->received_orbital_index)[i]; }
    Orbital &orbital(int i) { return (this->received_orbitals)[i]; }

    int get_size() const { return this->received_orbitals.size(); }
    int get_idx_sent(int i) const { return (this->sent_orbital_index)[i]; }
    int get_rank_sent(int i) const { return (this->sent_orbital_mpirank)[i]; }
    int get_step(int i) const { return (this->rcv_step)[i]; }
    int get_sent_size() const { return this->sent_orbital_index.size(); }

protected:
    const bool symmetric;
    int iter;
    int received_counter; // number of orbitals fetched during this iteration
    int sent_counter;     // number of orbitals sent during this iteration
    OrbitalVector *orbitals;
    OrbitalVector received_orbitals;
    std::vector<int> received_orbital_index;  // index of the orbitals received
    std::vector<int> sent_orbital_index;      // indices of the orbitals sent out
    std::vector<int> sent_orbital_mpirank;    // mpi rank (=rankID) of the orbitals sent out
    std::vector<int> rcv_step;                // for symmetric treatment: if the send was at step=0 or 1
    std::vector<std::vector<int>> orb_ix_map; // indices in the original orbital vector for each MPI
};

} // namespace mrchem
