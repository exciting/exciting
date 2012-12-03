# Define the list of reactions here:
# [('N2', -1), ('N', 2), ('reaction_energy', None), ('reaction_id', 1)] denotes:
# reaction_energy = -1 * N2 + 2 *N
# reaction_id is an integer or string to identify the reaction

reactions = [  [('N2', -1), ('N', 2), ('reaction_id', 1)],
               [('H2O', -1), ('H2', 1), ('O2', 0.5), ('reaction_id', 'H2O decomposition')],
               [('H2O', -1), ('H', 2), ('O', 1), ('reaction_id', 3)],
            ]

# reference reaction energies
reference = {
    # 'reaction_id' : reaction energy
    reactions[0][-1][1]: 9.93722251591,
    reactions[1][-1][1]: -0.496408928568,
    reactions[2][-1][1]: 9.14067431054,
    }
