"""
PyVRNA: A custom Python wrapper for ViennaRNA
Authors: Ayaan Hossain (ain.hoss07@gmail.com), Alex Reis
Trivia: This wrapper was engineered by Ayaan with lots of love and care during his
        third rotation in the BG program at Salis Lab between Feb-Mar 2017.
        ^Oh please we're still working on this wrapper.
"""

## For Cross-compatibility between Python 2 and Python 3


# Required for PyVRNA
import RNA
import math

from   itertools   import chain
from   collections import deque, namedtuple
from   operator    import itemgetter

# Required for ViennaRNA
import os
import re
import random
import string

_ACTIVE_PARAMETER_FILE = None

"""
NOTES:

Allowed characters in constraint strings:
.    (no pressure to base pair)
x    (must not base pair)
|    (must base pair, in either direction)
<, > (must base pair, in < or > direction)
(, ) (base marked ( must pair with its ) base)
i    (must pair with base in same molecule, intramolecular pairing)
e    (must pair with base in another molecule, intermolecular pairing)
+    (marked bases form a g-quadruplex)

Adding aptamers:
aptamers            = list of aptamer sequences in sequence (str)
aptamer_constraints = list of constraints (strs)
dG_ligands          = list of ligand binding free energies (double)

For a given aptamer, the sequence and constraint can be separated by "&" if
they are not a part of a contiguous sequence.
"""


class PyVRNA(object):
    """
    This class abstracts for RNAcentroid, RNAcofold, RNAeval, RNAfold and RNAsubopt programs in ViennaRNA suite, and is designed
    to be used by the various programs developed in Salis Lab.

    NOTE:  None of the (non-)sequence information is stored in PyVRNA object, because the parameters are extrinsic to the energy model
           that is abstracted by the PyVRNA object. All functions in this object operate on supplied sequence(s) and other parameters.
           The philosophy behind this design is to build a model once, at the beginning, and have it operate on all tge sequences and
           other parameters supplied in future, to retain consistency throughout. Also, testing inputs add overhead, and is not always
           necessary, so it must be turned on only when operating on artifically generated sequences or parameters that need testing.

    Usage: object_name = PyVRNA(temperature=float, dangles=int, gquad=bool, parameter_file=string, test_inputs=bool)
           object_name.function_name(parameters)
    """
    def __init__(self,
        temperature=37.0,
        dangles=2,
        noGU=False,
        noLP=False,
        gquad=False,
        parameter_file="rna_andronescu2007.par",
        duplex_adjustment=True,
        test_inputs=False,
        enforce_constraints=False,
        pyindex=False):
        """
        Initializes an PyVRNA object, after validating all parameters.

        Usage: energy_model = PyVRNA(temperature=int, dangles=int, gquad=bool, parameter_file=string, test_inputs=bool)
        """
        # Setup variables for helping in assertions
        self.parameter_files     = ["dna_mathews1999.par",
                                    "rna_turner1999.par",
                                    "dna_mathews2004.par",
                                    "rna_turner2004.par",
                                    "rna_andronescu2007.par"]
        self.parameter_directory =  "/opt/venv/share/ViennaRNA/" #"/usr/local/share/ViennaRNA/"

        # Assert all supplied variables for correctness
        assert 0 < temperature,                                'temperature must be greater than 0.'
        assert 0 <= dangles <= 3 and isinstance(dangles, int), 'dangles must be one of the three integers: 0 (none), 1 (some), 2 (all), or 3 (all + coaxial stacking of multi-branch loops)'
        assert isinstance(noGU, bool),                         'noGU must be a boolean value.'
        assert isinstance(noLP, bool),                         'noLP must be a boolean value.'
        assert isinstance(gquad, bool),                        'gquad must be a boolean value.'
        assert parameter_file in self.parameter_files,         ''.join(['parameter file must be ', " or ".join(self.parameter_files), '. ', 'Input was: {}'.format(parameter_file)])
        if parameter_file + ".par" in self.parameter_files:
            parameter_file += ".par"
        assert isinstance(duplex_adjustment, bool),            'duplex_adjustment must be a boolean value.'
        assert isinstance(test_inputs, bool),                  'test_inputs must be a boolean value.'
        assert isinstance(enforce_constraints, bool),          'enforce_constraints must be a boolean value.'
        assert isinstance(pyindex, bool),                      'pyindex must be a boolean value.'

        # Setup model default settings and variables
        self.settings             = RNA.md()
        self.settings.temperature = temperature   # Temperature in Celsius; default=37.0 (float)
        self.settings.dangles     = dangles       # Dangling end energies (0,1,2); see RNAlib documentation; default=2 (int)
        self.settings.gquad       = gquad         # Incorporate G-Quadruplex formation into structure prediction; default=off/0/False (bool/int)
        self.settings.noGU        = noGU          # Toggles GU wobble for RNA folding processes ... probably need to test this
        self.settings.noLP        = noLP
        RNA.cvar.temperature      = temperature   # Global setting of temperature
        RNA.cvar.dangles          = dangles       # Global setting of dangles
        RNA.cvar.gquad            = gquad         # Global setting of gquad
        RNA.cvar.noGU             = noGU         # Global setting of gquad
        RNA.cvar.noLP             = noLP

        '''
        Other settings variables available for manual overriding
        self.settings.betaScale       = 1.0       # Set the scaling of the Boltzmann factors; default=1.0 (float)
        self.settings.special_hp      = True      # Include tabulated free energies for special hairpin loops (Tri-, Tetra-, or Hexa-loops); default=on/1/True (bool/int)
        self.settings.noLP            = False     # Produce structures without lonely pairs (helices of length 1); default=off/0/False (bool/int)
        self.settings.noGU            = False     # Do not allow GU wobble pairs; default=off/0/False (bool/int)
        self.settings.noGUclosure     = False     # Do not allow GU pairs at the end of helices; default=off/0/False (bool/int)
        self.settings.logML           = False     # Recompute free energies of multi-branch loops using a logarithmic model; default=off/0/False (bool/int)
        self.settings.circ            = False     # Assume a circular (instead of linear) RNA molecule; default=off/0/False (bool/int)
        self.settings.canonicalBPonly = False     # Remove non-canonical base pairs from the structure constraint; default=off/0/False (bool/int)
        self.settings.uniq_ML         = False     # Create additional matrix for unique multi-branch loop prediction; default=off/0/False (bool/int)
        self.settings.energy_set      = 0         # Energy set, rarely used, see --energyModel; default=0 (int)
        self.settings.backtrack       = 1         # Whether to backtrack secondary structures; default=1 (int)
        self.settings.backtrack_type  = 'F'       # Set backtrack type, i.e. which DP matrix is used; default='F' (str)
        self.settings.compute_bpp     = True      # Compute base pair probabilities after partition function computation; default=on/1/True (bool/int)
        self.settings.max_bp_span     = -1        # Maximum base pair span; default=-1 (int)
        self.settings.min_loop_size   = 3         # Minimal loop size; default=3 (int)
        self.settings.window_size     = -1        # Window size for sliding window structure prediction approaches; default=-1 (int)
        self.settings.oldAliEn        = False     # Use old energy model for comparative structure prediction; default=off/0/False  (bool/int)
        self.settings.ribo            = False     # Use Ribosum Scoring in comparative structure prediction; default=off/0/False (bool/int)
        self.settings.cv_fact         = 1.0       # Co-variance scaling factor used in comparative structure prediction; default=1.0 (float)
        self.settings.nc_fact         = 1.0       # Unknown (not well defined in docs); likely for comparative structure prediction; default=1.0 (float)
        self.settings.sfact           = 1.07      # Scaling factor used to avodi under-/overflows in partition function computation; default=1.07 (float)
        '''

        # Setup model constraint variables
        self.enforce_constraints        = enforce_constraints
        self.constraints_options        = 16760832 # Not a mysterious value, see the macro values defined in '.../ViennaRNA-2.4.10/src/ViennaRNA/constraints/hard.h'
        self.pyindex                    = pyindex  # Global boolean for Python based 0-indexing or 1-indexing

        # Setup model parameter file and test_input variables
        self.parameter_file             = parameter_file
        parameter_path                  = self.parameter_directory + parameter_file
        global _ACTIVE_PARAMETER_FILE
        if _ACTIVE_PARAMETER_FILE != parameter_path:
            RNA.read_parameter_file(parameter_path)
            _ACTIVE_PARAMETER_FILE = parameter_path
        self.test_inputs                = test_inputs

        # Calculate dG_init adjustment for associating strands with specified temperature
        self.duplex_adjustment = duplex_adjustment
        self.dG_init_adjustment = self._calc_dG_init_adjustment()

        # Setup result structures as namedtuples
        self.PyVRNA_fold_result         = namedtuple('PyVRNA_fold_result',     'structure energy')
        self.PyVRNA_inverse_result      = namedtuple('PyVRNA_inverse_result',  'sequence distance')
        self.PyVRNA_centroid_result     = namedtuple('PyVRNA_centroid_result', 'structure energy distance')
        self.PyVRNA_ensemble_result     = namedtuple('PyVRNA_ensemble_result', 'structure energy')
        self.PyVRNA_bp_result           = namedtuple('PyVRNA_bp_result',      'length bpx bpy pkx pky gquad')

    def reset_profile(self):
        return None

    def snapshot_profile(self):
        return {}

    def cache_stats(self):
        return {}

    def _profile_add(self, name, start_ns):
        return None

    def _model_key(self):
        settings = self.settings
        return (
            self.parameter_file,
            settings.temperature,
            settings.dangles,
            bool(settings.gquad),
            bool(settings.noGU),
            bool(settings.noLP),
            bool(self.duplex_adjustment),
            bool(self.enforce_constraints),
            bool(self.pyindex),
        )

    def _freeze(self, value):
        if value is None:
            return None
        if isinstance(value, (tuple, list)):
            return tuple(self._freeze(v) for v in value)
        if isinstance(value, dict):
            return tuple(sorted((k, self._freeze(v)) for k, v in value.items()))
        return value

    def _sequence_cache_key(self, kind, *parts):
        return (kind, self._model_key()) + tuple(self._freeze(part) for part in parts)

    def _has_paired_nt(self, structure):
        return ('(' in structure) or (')' in structure) or ('[' in structure) or (']' in structure) or ('+' in structure)

    def _RNAsubopt_impl(self, sequences, constraints, delta_energy, aptamers, aptamer_constraints, dG_ligands, wrap_namedtuple):
        fc_obj = RNA.fold_compound("&".join(map(str, sequences)), self.settings)
        if constraints:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db("".join(constraints))
            else:
                fc_obj.hc_add_from_db("".join(constraints), self.constraints_options)

        if aptamers:
            for aptamer, fld, dG in zip(aptamers, aptamer_constraints, dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer, fld, dG, 0)

        solutions = fc_obj.subopt(delta_energy * 100)
        duplex_adjustment = self.duplex_adjustment and len(sequences) > 1
        dG_init_adjustment = self.dG_init_adjustment

        if wrap_namedtuple:
            make_result = self.PyVRNA_fold_result
            subopt_list = []
            for solution in solutions:
                structure = solution.structure
                energy = solution.energy
                if duplex_adjustment:
                    energy = energy + dG_init_adjustment if self._has_paired_nt(structure) else 0.0
                subopt_list.append(make_result(structure=structure, energy=energy))
        else:
            subopt_list = []
            for solution in solutions:
                structure = solution.structure
                energy = solution.energy
                if duplex_adjustment:
                    energy = energy + dG_init_adjustment if self._has_paired_nt(structure) else 0.0
                subopt_list.append((structure, energy))

        subopt_list.sort(key=itemgetter(1))
        return subopt_list

    def _test_sequences(self, sequence_list, func_name):
        """
        Tests if the sequences supplied to a particular function is valid for that function.

        NOTE:  This function is explictly called from other functions in this class that operates on sequences, when test_inputs=True.
               No external calls are required in order to ensure that sequences are valid.

        Usage: energy_model._test_sequences(sequence_list=list, func_name=string) # Validates sequences for func_name function
        """
        # Check if there are proper number of items in sequence_list
        len_sequence_dict = {'RNAcentroid': [1, 1], 'RNAcofold': [2, 2], 'RNAeval': [1, 2], 'RNAfold': [1, 1], 'RNAsubopt': [1, 2]}
        assert len_sequence_dict[func_name][0] <= len(sequence_list) <= len_sequence_dict[func_name][1], ''.join(['in ', func_name, ': ', ['at most', 'exactly'][not len_sequence_dict[func_name][0] < len_sequence_dict[func_name][1]], ' ', str(max(len_sequence_dict[func_name])), ' sequence(s) are accepted.'])

        # Check if first item is a string
        assert isinstance(sequence_list[0], str), 'in '+func_name+': sequence[_1] must be a string.'

        # Check if second item (if applicable) is also a string
        second_seq_funcs_set = {'RNAcofold', 'RNAeval', 'RNAsubopt'}
        if func_name in second_seq_funcs_set and len(sequence_list) > 1:
            assert isinstance(sequence_list[1], str), 'in {}: {}sequence[_2] must be a string.'.format(func_name, ['(optional) ', ''][func_name == 'RNAcofold'])

        # Check if all sequences are valid
        valid_charset = set('ATGCUatgcu&')
        for sequence in sequence_list:
            if not sequence is None:
                assert set(sequence) <= valid_charset, ''.join(['in ', func_name, ': sequence string(s) must only contain ATGCUatgcu& characters only.'])

    def _test_non_sequences(self, non_sequence_list, sequence_list, func_name, test_for='structure'):
        """
        Tests if the non-sequences supplied to a particular function is valid for that function. Non-sequences must be 'structure' or 'constraint'

        NOTE:  This function is explictly called from other functions in this class that operates on non-sequences, when test_inputs=True.
               No external calls are required in order to ensure that non-sequences are valid.

        Usage: energy_model._test_non_sequences(non_sequence_list=list, sequence_list=list, func_name=string, test_for=string) # Validates non-sequences for func_name function
        """

        if not non_sequence_list:
            return

        # Check if test_for is structures or constraints
        assert test_for == 'structure' or test_for == 'constraint', ''.join(['in ', func_name, ": test_for must be 'constraint' or 'structure'"])

        # Check if there are proper number of items in non_sequence_list
        assert 1 <= len(non_sequence_list) <= 2, ''.join(['in ', func_name, ': only one or two structures(s) are accepted.'])

        # Check if structures types are compatible with sequence types
        assert isinstance(non_sequence_list[0], str) or isinstance(non_sequence_list[0], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_1] must be a string or a None object.'])
        if isinstance(non_sequence_list[0], type(None)):
            assert isinstance(non_sequence_list[1], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_2] must be a None object since ', test_for,'[_1] is also None.'])
        if isinstance(sequence_list[1], str):
            assert isinstance(non_sequence_list[1], str) or isinstance(non_sequence_list[1], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_2] must be a string or a None object.'])
        elif isinstance(sequence_list[1], type(None)):
            assert isinstance(non_sequence_list[1], type(None)), ''.join(['in ', func_name, ': ', test_for, '[_2] must be a None object since sequence[_2] is also None.'])

        # Check if structures are valid, and their lengths are same as their sequences
        valid_chars   = '.([+])&' if test_for == 'structure' else '.x([+])&'
        valid_charset = set(valid_chars)
        for i in range(2):
            string_i = str(i+1)
            if not (non_sequence_list[i] is None or sequence_list[i] is None):
                assert len(non_sequence_list[i]) == len(sequence_list[i]), ''.join(['in ', func_name, ': ', test_for,'_', string_i, ' and sequence_', string_i,' must be of same length.'])
                assert set(non_sequence_list[i]) <= valid_charset,         ''.join(['in ', func_name, ': ', test_for,'_', string_i, ' must be a string of ', valid_chars,' characters only.'])

    def _test_aptamer_inputs(self, aptamers, aptamer_constraints, dG_ligands):
        """
        Tests if the aptamer constraints supplied to a particular function are valid for that function.

        NOTE:  This function is explictly called from other functions in this class that operates on non-sequences, when test_inputs=True.
               No external calls are required in order to ensure that non-sequences are valid.

        Usage: energy_model._test_non_sequences(non_sequence_list=list, sequence_list=list, func_name=string, test_for=string) # Validates non-sequences for func_name function
        """

        if not aptamers and not aptamer_constraints:
            return

        assert len(aptamers)==len(aptamer_constraints), "Make sure to provide equal number of aptamers and aptamer_constraint values"
        assert len(aptamers)==len(dG_ligands), "Make sure to provide equal number of aptamers and dG_ligand values"
        for aptamer in aptamers:
            self._test_sequences(aptamer)
        assert all(isinstance(dG,float) for dG in dG_ligands), "Make sure dG_ligand values are floats in kcal/mol"
        self._test_non_sequences(aptamer_constraints,aptamers,'RNAfold','constraint')

    def _test_bp_tuple(self, length, bpx, bpy, pkx, pky, gquad, func_name):
        """
        Tests if all inputs for a bp_tuple are valid.

        Usage: energy_model._test_bp_tuple(bp_tuple.bpx, bp_tuple.bpy, bp_tuple.pkx, bp_tuple.pky, bp_tuple.gquad, func_name=string) # Validates bpx/bpy, pkx/pky and gquad in func_name
        """
        assert isinstance(length,list),                         ''.join(['in ', func_name, ': length must be an integer.'])
        assert length >= max([len(bpx), len(pkx), len(gquad)]), ''.join(['in ', func_name, ': length must be greater than or equal to the length of bpx/bpy, pkx/pky and quad lists.'])
        assert len(bpx) == len(bpy),                            ''.join(['in ', func_name, ': bpx and bpy must have same length.'])
        assert len(pkx) == len(pky),                            ''.join(['in ', func_name, ': pkx and pky must have same length.'])
        for index, gquad_tuple in enumerate(gquad):
            assert len(gquad_tuple) == 4, ''.join(['in ', func_name, ': gquad[', str(index), '] must contain 4 indices.'])

    def _calc_dG_init_adjustment(self):
        """
        Adjustment according to Dirks et al., Thermodynamic Analysis of Interacting Nucleic Acid Strands
        See footnote 13: "Based on dimensional analysis, we define our concentrations as mole fractions rather than
        molarities. Therefore, the free energy of strand association for a complex of L strands is..."

        Based on a partition function analysis of dilute solutions of interacting strands
        in a fixed volume (the "box")

        NOTE:  This function is explictly called from other functions in this class that returns minimum free energy of a duplex.
               No external calls are required in order to ensure that energies are adjusted.
        """

        kB = 0.00198717 # Boltzmann constant in kcal/mol/K
        T = self.settings.temperature
        a = [-3.983035, 301.797, 522528.9, 69.34881, 999.974950]

        # Calculate the number of moles of water per liter (molarity) at temperature (T in deg C)
        # Density of water calculated using data from
        # Tanaka M., Girard, G., Davis, R., Peuto A., Bignell, N.
        # Recommended table for the density of water..., Metrologia, 2001, 38, 301-309
        pH2O = a[4] * (1 - (T+a[0])**2.0*(T+a[1])/a[2]/(T+a[3])) / 18.0152

        return -kB*(T+273.15)*math.log(pH2O)

    def RNAcentroid(self, sequence, constraint=None, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the centroid structure, its energy and distance for a sequence.

        Usage: centroid_result = energy_model.RNAcentroid(sequence, constraint)  # Executes RNAcentroid
               centroid_result.structure                                         # Retrieves centroid structure
               centroid_result.energy                                            # Retrieves centroid energy
               centroid_result.distance                                          # Retrieves centroid distance
               centroid_structure = energy_model.RNAcentroid(sequence).structure # Executes RNAcentroid and retrieves structure only
        """
        # Test if sequence is valid
        if self.test_inputs:
            self._test_sequences(sequence_list=[sequence], func_name='RNAcentroid')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        # Setup ViennaRNA library objects and call centroid function
        fc_obj = RNA.fold_compound(str(sequence), self.settings)
        if not constraint is None:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db(constraint)
            else:
                fc_obj.hc_add_from_db(constraint, self.constraints_options)
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT
        fc_obj.pf()
        (structure, distance) = fc_obj.centroid()

        # Return the centroid structure, energy and distance
        return self.PyVRNA_centroid_result(structure=structure, energy=fc_obj.eval_structure(structure), distance=distance)

    def RNAcofold(self, sequences, constraints=[], aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the minimum free energy structure and its corresponding energy for a cofold.
        """
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAcofold')
            self._test_non_sequences(non_sequence_list=constraints, sequence_list=sequences, func_name='RNAcofold', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        fc_obj = RNA.fold_compound("&".join(map(str,sequences)), self.settings)
        if constraints:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db("".join(constraints))
            else:
                fc_obj.hc_add_from_db("".join(constraints), self.constraints_options)
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0)
        structure, energy = fc_obj.mfe_dimer()

        if self.duplex_adjustment and set('.') < set(structure):
            energy += self.dG_init_adjustment

        structure = structure[:len(sequences[0])] + "&" + structure[len(sequences[0]):]
        return self.PyVRNA_fold_result(structure=structure, energy=energy)
    def RNAduplex(self, sequences):
        """
        Computes the minimum free energy and structure of an RNA duplex (intermolecular base pairs only!).

        *** FOR NOW RNAduplex DOES NOT SUPPORT ADDING CONSTRAINTS ***

        Usage: duplex_result = energy_model.RNAduplex(sequences=[sequence1,sequence2])  # Executes RNAduplex
               duplex_result.structure                                                  # Retrieves mfe structure
               duplex_result.energy                                                     # Retrieves mfe
               duplex_energy = energy_model.RNAduplex(sequence_1, sequence_2).energy    # Executes RNAduplex and retrieves the mfe energy only
        """
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAcofold')

        # Setup ViennaRNA library objects and call mfe_dimer (cofold) function
        fc_obj = RNA.fold_compound("&".join(map(str,sequences)), self.settings)
        fc_obj.hc_add_from_db('e'*sum(map(len,sequences)),4194304)
        structure, energy = fc_obj.mfe_dimer()

        if self.duplex_adjustment and set('.') < set(structure):
            energy += self.dG_init_adjustment

        structure = structure[:len(sequences[0])] + "&" + structure[len(sequences[0]):]

        # Return the duplex structure and energy
        return self.PyVRNA_fold_result(structure=structure, energy=energy)

    def RNAensemble(self, sequence, constraint=None, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the ensemble structure, the Gibbs free energy based on the partition fxn and distance.

        Usage: ensemble_result = energy_model.RNAensemble(sequence, constraint)  # Executes RNAensemble
               ensemble_result.structure                                         # Retrieves ensemble structure
               ensemble_result.energy                                            # Retrieves ensemble energy
               ensemble_structure = energy_model.RNAcentroid(sequence).structure # Executes RNAensemble and retrieves structure only
        """
        # Test if sequence is valid
        if self.test_inputs:
            self._test_sequences(sequence_list=[sequence], func_name='RNAensemble')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        # Setup ViennaRNA library objects and call centroid function
        fc_obj = RNA.fold_compound(str(sequence), self.settings)
        if not constraint is None:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db(constraint)
            else:
                fc_obj.hc_add_from_db(constraint, self.constraints_options)
        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0) # 0 is VRNA_OPTION_DEFAULT
        (structure, energy) = fc_obj.pf()

        # Return the ensemble structure, energy and distance
        return self.PyVRNA_ensemble_result(structure=structure, energy=energy)

    def RNAeval(self, sequences, structures, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the minimum free energy of a given structure and its corresponding sequence or duplex.
        """
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAeval')
            self._test_non_sequences(non_sequence_list=structures, sequence_list=sequences, func_name='RNAeval', test_for='structure')

        energy = RNA.fold_compound("&".join(map(str,sequences)), self.settings).eval_structure("".join(structures))

        if self.duplex_adjustment and len(sequences) > 1 and set('.') < set("".join(structures)):
            energy += self.dG_init_adjustment

        return energy
    def RNAfold(self, sequence, constraint=None, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the minimum free energy structure and its corresponding energy for a sequence.
        """
        if self.test_inputs:
            self._test_sequences(sequence_list=[sequence], func_name='RNAfold')
            self._test_non_sequences(non_sequence_list=[constraint, None], sequence_list=[sequence, None], func_name='RNAfold', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)

        fc_obj = RNA.fold_compound(str(sequence), self.settings)
        if not constraint is None:
            if not self.enforce_constraints:
                fc_obj.hc_add_from_db(constraint)
            else:
                fc_obj.hc_add_from_db(constraint, self.constraints_options)

        if aptamers:
            for aptamer,fld,dG in zip(aptamers,aptamer_constraints,dG_ligands):
                fc_obj.sc_add_hi_motif(aptamer,fld,dG,0)

        structure, energy = fc_obj.mfe()
        return self.PyVRNA_fold_result(structure=structure, energy=energy)
    def RNAinverse(self, sequence, structure):
        """
        Computes a sequence conforming to given minimum free energy structure starting from given sequence and its distance to starting sequence

        Usage: inverse_result = energy_model.RNAinverse(sequence, structure)           # Executes RNAinverse
               inverse_result.structure                                                # Retrieves inverted sequence
               inverse_result.distance                                                 # Retrieves distance from startng sequence
               inverse_sequence = energy_model.RNAfold(sequence, constraint).sequence  # Executes RNAinverse and retrieves the sequence only
        """

        # Test if sequence and constraint is valid
        # Coming soon

        sequence, distance = RNA.inverse_fold(str(sequence), structure)

        # Return the sequence and the distance
        return self.PyVRNA_inverse_result(sequence=sequence, distance=distance)

    def RNAsubopt(self, sequences, constraints=[], delta_energy=5, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        """
        Computes the suboptimal structures for a sequence or duplex with optional constraints within a delta_energy range from its mfe.
        """
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAsubopt')
            self._test_non_sequences(non_sequence_list=constraints, sequence_list=sequences, func_name='RNAsubopt', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)
        delta_energy = int(round(delta_energy))
        return self._RNAsubopt_impl(sequences, constraints, delta_energy, aptamers, aptamer_constraints, dG_ligands, wrap_namedtuple=True)

    def RNAsubopt_tuples(self, sequences, constraints=[], delta_energy=5, aptamers=[], aptamer_constraints=[], dG_ligands=[]):
        if self.test_inputs:
            self._test_sequences(sequence_list=sequences, func_name='RNAsubopt')
            self._test_non_sequences(non_sequence_list=constraints, sequence_list=sequences, func_name='RNAsubopt', test_for='constraint')
            self._test_aptamer_inputs(aptamers,aptamer_constraints,dG_ligands)
        delta_energy = int(round(delta_energy))
        return self._RNAsubopt_impl(sequences, constraints, delta_energy, aptamers, aptamer_constraints, dG_ligands, wrap_namedtuple=False)
    def create_bp_tuple(self, length=[], bpx=[], bpy=[], pkx=[], pky=[], gquad=[]):
        """
        Returns a customized bp_tuple from input bpx, bpy, pkx, pky and gquad lists.

        NOTE:  bp_tuple is a tuple of (length, base_pair_x_list, base_pair_y_list, pseudo_pair_x_list, pseudo_pair_y_list, gquad_list).

        Usage: custom_bp_tuple = energy_model.create_bp_tuple(bpx=some_list, bpy=some_list) # Executes create_bp_tuple and returns the specified bp_tuple
        """
        # Test if the inputs are valid
        if self.test_inputs:
            self._test_bp_tuple(length=length, bpx=bpx, bpy=bpy, pkx=pkx, pky=pky, gquad=gquad, func_name='create_bp_tuple')

        # Create and return the customized bp_tuple
        return self.PyVRNA_bp_result (length=length, bpx=list(bpx), bpy=list(bpy), pkx=list(pkx), pky=list(pky), gquad=list(gquad))

    def vienna2bp(self, vienna_string, pyindex=None):
        """
        Converts the vienna_string (dot bracket string representing a mfe structure) to bp_tuple.
        """
        if pyindex is None:
            pyindex = self.pyindex

        if self.test_inputs:
            assert set(vienna_string) <= set('&.([+])'), 'in vienna2bp: vienna_string must be a string of .([+]) characters only.'
            assert vienna_string.count('&') <= 1, 'Odd input. Should only have 1 & if a cofold, you provided: {}.'.format(vienna_string)

        if '&' in vienna_string:
            split_string = vienna_string.split('&')
            length = [len(x) for x in split_string]
            compact = "".join(split_string)
        else:
            length = [len(vienna_string)]
            compact = vienna_string

        if ('[' not in compact) and ('+' not in compact):
            stack = []
            pairs = []
            for index, char in enumerate(compact, start=1):
                if char == '(':
                    stack.append(index)
                elif char == ')':
                    pairs.append((stack.pop(), index))
            pairs.sort(key=itemgetter(0))
            if pyindex:
                bpx = [i - 1 for i, _ in pairs]
                bpy = [j - 1 for _, j in pairs]
            else:
                bpx = [i for i, _ in pairs]
                bpy = [j for _, j in pairs]
            return self.create_bp_tuple(length=length, bpx=bpx, bpy=bpy)

        bp_tuple = self.create_bp_tuple(length=length)
        bp_stack, pk_stack = [], []
        has_gquad = False
        gquad_matrix = {0:[], 1:[], 2:[], 3:[]}
        string_index = 0
        gquad_index = 0
        bp_pairs = []
        pk_pairs = []

        for index, char in enumerate(compact, start=1):
            if char == '(':
                bp_stack.append(index)
            elif char == '[':
                pk_stack.append(index)
            elif char == ')':
                bp_pairs.append((bp_stack.pop(), index))
            elif char == ']':
                pk_pairs.append((pk_stack.pop(), index))
            elif char == '+':
                has_gquad = True

        bp_pairs.sort(key=itemgetter(0))
        pk_pairs.sort(key=itemgetter(0))
        bp_tuple.bpx.extend(i for i, _ in bp_pairs)
        bp_tuple.bpy.extend(j for _, j in bp_pairs)
        bp_tuple.pkx.extend(i for i, _ in pk_pairs)
        bp_tuple.pky.extend(j for _, j in pk_pairs)

        while has_gquad and string_index < len(compact):
            if compact[string_index] == '+' and (gquad_index == 0 or gquad_index == 2):
                if gquad_index == 0:
                    gquad_matrix[0].append(deque()); gquad_matrix[1].append(deque()); gquad_matrix[2].append(deque()); gquad_matrix[3].append(deque())
                while string_index < len(compact) and compact[string_index] == '+':
                    gquad_matrix[gquad_index][-1].append(string_index+1)
                    string_index += 1
                gquad_index = (gquad_index + 1) % 4
            elif compact[string_index] == '+' and (gquad_index == 1 or gquad_index == 3):
                while string_index < len(compact) and compact[string_index] == '+':
                    gquad_matrix[gquad_index][-1].appendleft(string_index+1)
                    string_index += 1
                gquad_index = (gquad_index + 1) % 4
            else:
                string_index += 1
        else:
            if has_gquad:
                for i in range(len(gquad_matrix[0])):
                    bp_tuple.gquad.extend(list(zip(gquad_matrix[0][i], gquad_matrix[1][i], gquad_matrix[2][i], gquad_matrix[3][i])))

        if pyindex:
            bpx = [i-1 for i in bp_tuple.bpx]
            bpy = [j-1 for j in bp_tuple.bpy]
            pkx = [i-1 for i in bp_tuple.pkx]
            pky = [j-1 for j in bp_tuple.pky]
            gquad = [tuple(j-1 for j in i) for i in bp_tuple.gquad]
            bp_tuple = self.create_bp_tuple(length,bpx,bpy,pkx,pky,gquad)
        return bp_tuple
    def vienna2bplite(self, vienna_string, pyindex=None):

        #UH OH. RNA.ptable_from_string no longer exists in ViennaRNA v2.5
        return self.vienna2bp(vienna_string, pyindex=None)

    def bp2vienna(self, length, bpx=[], bpy=[], pkx=[], pky=[], gquad=[], pyindex=None):
        """
        Creates a customized bp_tuple from input bpx, bpy, pkx, pky and gquad lists and returns the vienna_string encoded by it.

        NOTE:  bp_tuple is a tuple of (length, base_pair_x_list, base_pair_y_list, pseudo_pair_x_list, pseudo_pair_y_list, gquad_list).

        Usage: custom_vienna_string = energy_model.bp2vienna(bpx=some_list, bpy=some_list) # Executes bp2vienna and returns the vienna_string
        """
        # Setup indexing
        if pyindex is None:
            pyindex = self.pyindex

        return self.bptuple2vienna(self.create_bp_tuple(length, bpx, bpy, pkx, pky, gquad))

    def bptuple2vienna(self, bp_tuple, pyindex=None):
        """
        Converts an bp_tuple to vienna_string (dot bracket string representing a mfe structure).
        """
        if pyindex is None:
            pyindex = self.pyindex

        if self.test_inputs:
            self._test_bp_tuple(length=bp_tuple.length, bpx=bp_tuple.bpx, bpy=bp_tuple.bpy, pkx=bp_tuple.pkx, pky=bp_tuple.pky, gquad=bp_tuple.gquad, func_name='bptuple2vienna')

        if pyindex:
            bpx = [i+1 for i in bp_tuple.bpx]
            bpy = [j+1 for j in bp_tuple.bpy]
            pkx = [i+1 for i in bp_tuple.pkx]
            pky = [j+1 for j in bp_tuple.pky]
            gquad = [tuple(j+1 for j in i) for i in bp_tuple.gquad]
            bp_tuple = self.create_bp_tuple(bp_tuple.length,bpx,bpy,pkx,pky,gquad)

        vienna_string_list = ['.'] * sum(bp_tuple.length)
        for i, j in zip(bp_tuple.bpx, bp_tuple.bpy):
            if i >= 0 and j >= 0:
                vienna_string_list[i-1] = '('
                vienna_string_list[j-1] = ')'
            else:
                break
        for i, j in zip(bp_tuple.pkx, bp_tuple.pky):
            if i >= 0 and j >= 0:
                vienna_string_list[i-1] = '['
                vienna_string_list[j-1] = ']'
            else:
                break
        for quad in bp_tuple.gquad:
            for idx in quad:
                if idx >= 0:
                    vienna_string_list[idx-1] = '+'

        if len(bp_tuple.length) == 2:
            vienna_string_list.insert(bp_tuple.length[0],'&')

        return "".join(vienna_string_list)

# Legacy
class ViennaRNA(dict):
    # Legacy
    def __init__(self, Sequence_List, material, Gquad=False):
        # Legacy: Check if sequences in Sequence_List are valid
        exp = re.compile('[ATGCU]',re.IGNORECASE)
        for seq in Sequence_List:
            if exp.match(seq) is None:
                raise ValueError("Invalid letters found in inputted sequences. Only ATGCU allowed. \n Sequence is \"" + str(seq) + "\".")

        # Legacy: Setting object parameters
        self.ran                = 0
        self["sequences"]       = Sequence_List
        self["material"]        = material
        parameter_files         = ["dna_mathews1999", "rna_turner1999", "dna_mathews2004", "rna_turner2004", "rna_andronescu2007"]
        material_par            = material+".par"
        self["RNA_model_param"] = material_par if material in parameter_files else parameter_files[-1]+".par"
        self["Gquad"]           = Gquad
        self["Gquad_param"]     = ["-g"] if Gquad else []

        # New: Refactoring dangle setup for other functions
        self.dangles_dict = {'all':2, 'some':1, 'none': 0}

    # Legacy
    def centroid(self, strands, constraint=None, Temp=37.0, dangles="all", outputPS=False):
        # Legacy: Checks and setup
        if Temp <= 0:        raise ValueError("The specified temperature must be greater than zero.")
        if len(strands) > 1: raise ValueError("Two RNA strands are inputted. ViennaRNA does NOT return Centroid for RNAcofold.")
        self["Centroid_composition"] = strands

        # New: PyVRNA execution
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        structure, energy, distance = energy_model.RNAcentroid(sequence=self["sequences"][0])

        # Legacy: Parsing and storing output
        bp_tuple = energy_model.vienna2bp(structure, pyindex=False)
        self["program"]                 = "Centroid"
        self["totalnt"]                 = bp_tuple.length
        self["Centroid_energy"]         = [energy]
        self['Centroid_bracket_string'] = structure
        self["Centroid_basepairing_x"]  = [bp_tuple.bpx]
        self["Centroid_basepairing_y"]  = [bp_tuple.bpy]

    # Legacy
    def convert_bracket_to_numbered_pairs(self, bracket_string):
        # New: PyVRNA execution
        bp_tuple = PyVRNA(test_inputs=False).vienna2bp(bracket_string, pyindex=False)
        return [bp_tuple.length, bp_tuple.bpx, bp_tuple.bpy, bp_tuple.pkx, bp_tuple.pky]

    # Legacy
    def convert_numbered_pairs_to_bracket(self, strands, bp_x, bp_y, PK_bp_x=[], PK_bp_y=[], Gquad_bp=[]):
        # New: PyVRNA execution
        energy_model = PyVRNA(test_inputs=False)
        bp_tuple = energy_model.PyVRNA_bp_result(length=strands, bpx=bp_x, bpy=bp_y, pkx=PK_bp_x, pky=PK_bp_y, gquad=Gquad_bp)
        return energy_model.bptuple2vienna(bp_tuple,pyindex=False)

    # Legacy
    def energy(self, strands, base_pairing_x, base_pairing_y, Temp=37.0, dangles="all"):
        # Legacy: Checks and setup
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        self["energy_composition"] = strands
        strands                    = list(map(len, self["sequences"]))
        sequences                  = self["sequences"] if len(strands) == 2 else [self["sequences"][0]]

        # New: PyVRNA execution
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        bp_tuple     = energy_model.PyVRNA_bp_result(length=strands, bpx=base_pairing_x, bpy=base_pairing_y, pkx=[], pky=[], gquad=[])
        energy       = energy_model.RNAeval(sequences, energy_model.bptuple2vienna(bp_tuple, pyindex=False).split('&'))

        # Legacy: Parsing and storing output
        self["program"]              = "energy"
        self["energy_energy"]        = [energy]
        self["energy_basepairing_x"] = [base_pairing_x]
        self["energy_basepairing_y"] = [base_pairing_y]
        return energy

    # Legacy
    def mfe(self, strands, constraints=None, Temp=37.0, dangles="all", outputPS=False, duplex=False):
        # Legacy: Checks and setup
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        self["mfe_composition"] = strands

        # New: PyVRNA execution
        constraint = constraints.split('&') if (not constraints is None and '&' in constraints) else constraints
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        if   len(strands) == 1:
            structure, energy = energy_model.RNAfold(sequence=self["sequences"][0], constraint=constraint)
        elif len(strands) == 2:
            structure, energy = energy_model.RNAcofold(self["sequences"], constraints=constraint)
        else:
            raise ValueError("Three RNA strands are inputted. ViennaRNA does NOT return structure and energy for three sequences in RNA(co)fold.")

        # Legacy: Parsing and storing output
        bp_tuple = energy_model.vienna2bp(structure, pyindex=False)
        self["program"]            = "mfe"
        self["totalnt"]            = bp_tuple.length
        self["mfe_energy"]         = [energy]
        self["mfe_bracket_string"] = structure
        self["mfe_basepairing_x"]  = [bp_tuple.bpx]
        self["mfe_basepairing_y"]  = [bp_tuple.bpy]

    # Legacy
    def subopt(self, strands, energy_gap, Temp=37.0, dangles="all", constraints=None, outputPS=False):
        # Legacy: Checks and setup
        if Temp <= 0: raise ValueError("The specified temperature must be greater than zero.")
        self["subopt_composition"]   = strands
        self["subopt_energy"]        = []
        self["subopt_basepairing_x"] = []
        self["subopt_basepairing_y"] = []

        # New: PyVRNA execution
        constraint = constraints.split('&') if (not constraints is None and '&' in constraints) else constraints
        energy_model = PyVRNA(temperature=Temp, dangles=self.dangles_dict[dangles], gquad=self["Gquad"], parameter_file=self["RNA_model_param"], test_inputs=False)
        if   len(strands) == 1:
            results = energy_model.RNAsubopt(sequences=[self["sequences"][0]], constraints=constraint, delta_energy=energy_gap)
        elif len(strands) == 2:
            results = energy_model.RNAsubopt(sequences=self["sequences"], constraints=constraint, delta_energy=energy_gap)
        else:
            raise ValueError("Three RNA strands are inputted. ViennaRNA does NOT return suboptimal structures and energies for three sequences in RNAsubopt.")

        # Legacy: Parsing and storing output
        for structure, energy in results:
            bp_tuple = energy_model.vienna2bp(structure,pyindex=False)
            self["subopt_energy"].append(energy)
            self["subopt_basepairing_x"].append(bp_tuple.bpx)
            self["subopt_basepairing_y"].append(bp_tuple.bpy)
        self["program"]           = "subopt"
        self["totalnt"]           = strands
        self["subopt_NumStructs"] = len(results)


def tests():
    """
    Tests each function in PyVRNA class
    """
    print("### === Testing PyVRNA class === ###")

    # Tests for parameter files
    print('Testing all parameter files...')
    sequence     = "CGCAGGGAUACCCGCG"
    parameter_files = ["dna_mathews1999.par", "rna_turner1999.par", "dna_mathews2004.par", "rna_turner2004.par", "rna_andronescu2007.par"]
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[0])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -3.200000047683716]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[1])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.((....)).)))', -5.5]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[2])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -3.9000000953674316]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[3])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -5.599999904632568]
    energy_model = None
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file=parameter_files[4])
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -4.619999885559082]
    energy_model = None
    print('Good.')

    # Tests for RNAcentroid
    print('Testing RNAcentroid...')
    sequence     = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    cntrd_result = energy_model.RNAcentroid(sequence)
    assert [cntrd_result.structure, cntrd_result.energy, cntrd_result.distance] == ['(((.(.((.......))..))))....', 1.399999976158142, 3.345900802497659]
    energy_model = None

    sequence     = 'CGCAGGGAUACCCGCGGCGCCCAUAGGGACGC'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    cntrd_result = energy_model.RNAcentroid(sequence)
    assert [cntrd_result.structure, cntrd_result.energy, cntrd_result.distance] == ['(((.(((...))))))((((((...))).)))', -13.5, 2.4947844957445353]
    energy_model = None
    print('Good.')

    print( 'Testing RNAcentroid with aptamers...')
    # example 1
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACG"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[theophylline_aptamer],aptamer_constraints=[theophylline_constraint],dG_ligands=[dG_theophylline])
    assert [fold_result.structure, fold_result.energy] == ['.....((((..(.(((((...((((((((.....)))))...)))...))))).)..))))......',-17.0]
    fold_result  = energy_model.RNAcentroid(sequence)
    assert [fold_result.structure, fold_result.energy] == ['.....((((..(.(((((...((((((((.....)))))...)))...))))).)..))))......',-17.0]
    energy_model = None

    # example 2
    sequence = "CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    streptavidin_aptamer =              'CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGG'
    streptavidin_aptamer_structure =    '((xxxxxxxxxxxxx((((xxxxxxxxx))))))'
    dG_streptavidin = -10.15
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[streptavidin_aptamer],aptamer_constraints=[streptavidin_aptamer_structure],dG_ligands=[dG_streptavidin])
    assert [fold_result.structure, fold_result.energy] == ['((.............((((.........))))))................................',-0.4000000059604645]
    fold_result  = energy_model.RNAcentroid(sequence)
    assert [fold_result.structure, fold_result.energy] == ['................((((((...((((...))))....))))))....................',-9.600000381469727]

    # example 3 - add multiple apatmers
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[0.0,0.0]) # what if the aptamer provides no additional free energy?
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))((((((((((((........))).)))))......((....)).............))))........)))))).................',-29.399999618530273]
    fold_result  = energy_model.RNAcentroid(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[dG_theophylline,dG_streptavidin]) # now with the aptamer binding free energies
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))(((...)))(((........)))((((.............((((.........)))))))).......)))))).................',-24.799999237060547]
    print( 'Good.' )

    # Tests for RNAcofold
    print( 'Testing RNAcofold...' )
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    energy_model = PyVRNA()
    fold_result  = energy_model.RNAcofold(sequences)
    assert [fold_result.structure, fold_result.energy] == ['.((.(((...))))).&((((((...))).)))', -10.25 + energy_model.dG_init_adjustment]
    energy_model = None
    sequences    = ["GCGCACAUAGUGACGC","GCGCCCAUAGGGACGC"]
    constraints  = ["....xx....xx....","....xx....xx...."]
    energy_model = PyVRNA(dangles=1, gquad=False, parameter_file='rna_andronescu2007.par')
    fold_result  = energy_model.RNAcofold(sequences,constraints)
    assert [fold_result.structure, fold_result.energy] == ['((((............&))))............', -5.869999885559082 + energy_model.dG_init_adjustment]
    energy_model = None
    print( 'Good.' )

    # Tests for RNAensemble
    print( 'Testing RNAensemble...' )
    sequence     = 'CGCAGGGAUACCCGCGGCGCCCAUAGGGACGCCGCAGGGAUACCCGCGGCGCCCAUAGGGACGC'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    ensemble_result = energy_model.RNAensemble(sequence)
    assert [ensemble_result.structure, ensemble_result.energy] == ['(((.(((...||||||((((((.,.{||.||||||.|||...))))))}||}|}.,.))).)))', -39.47269058227539]
    energy_model = None
    print( 'Good.' )

    print( 'Testing ensemble: RNAensemble with aptamers...' )
    # example 1
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACG"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[theophylline_aptamer],aptamer_constraints=[theophylline_constraint],dG_ligands=[dG_theophylline])
    assert [fold_result.structure, fold_result.energy] == [',,.,.((((.,{.(((((...((((((((,...,)))))...)))...))))).},.)))),,,...',-28.491722106933594]
    fold_result  = energy_model.RNAensemble(sequence)
    assert [fold_result.structure, fold_result.energy] == [',,.,.((((.,{.(((((...((((((((,...,)))))...)))...))))).},.)))),,,...',-19.407447814941406]
    energy_model = None

    # example 2
    sequence = "CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    streptavidin_aptamer =              'CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGG'
    streptavidin_aptamer_structure =    '((xxxxxxxxxxxxx((((xxxxxxxxx))))))'
    dG_streptavidin = -10.15
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[streptavidin_aptamer],aptamer_constraints=[streptavidin_aptamer_structure],dG_ligands=[dG_streptavidin])
    assert [fold_result.structure, fold_result.energy] == ['{{.............{(((,,,...,,,|}}})),,....},,,,,....................',-11.105342864990234]
    fold_result  = energy_model.RNAensemble(sequence)
    assert [fold_result.structure, fold_result.energy] == ['.........,,.....((((((...((((...))))....))))))..,.................',-10.555535316467285]

    # example 3 - add multiple apatmers
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[0.0,0.0]) # what if the aptamer provides no additional free energy?
    assert [fold_result.structure, fold_result.energy] == ['..,.((((((((((((((((......)))).))))))(((({(((((((........))).))))}.,,,,{((....}|}}..}},...,..))))........)))))).,,..............',-33.18623352050781]
    fold_result  = energy_model.RNAensemble(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[dG_theophylline,dG_streptavidin]) # now with the aptamer binding free energies
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......}))).))))))(((...))){((........)))((((.............((((.........)))))))).......)))))).,...............',-35.90715408325195]
    print( 'Good.' )

    # Tests for RNAeval
    print( 'Testing RNAeval...' )
    sequence     = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    energy  = energy_model.RNAeval([sequence], ['(((.(.((.......))..))))....'])
    assert float("{0:0.2f}".format(energy)) == 2.15
    energy_model = None
    sequences    = ['CGCAGGGAUACCCGCG','GCGCCCAUAGGGACGC']
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_andronescu2007.par')
    energy  = energy_model.RNAeval(sequences,['.((.(((...))))).','((((((...))).)))'])
    assert float("{0:0.2f}".format(energy)) == -12.72
    energy_model = None
    print( 'Good.' )

    # Tests for RNAfold
    print( 'Testing RNAfold...' )
    sequence     = "CGCAGGGAUACCCGCG"
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['(((.(((...))))))', -4.619999885559082]
    energy_model = None
    sequence     = "CGCAGGGAUACCCGCGGCGCCCAUAGGGACGC"
    constraints  = '.....xxx........................'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence, constraints)
    assert [fold_result.structure, fold_result.energy] == ['(((.(......).)))((((((...))).)))', -9.300000190734863]
    energy_model = None
    print( 'Good.' )

    print( 'Testing mfe: RNAfold with aptamers...' )
    # example 1
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACG"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence,aptamers=[theophylline_aptamer],aptamer_constraints=[theophylline_constraint],dG_ligands=[dG_theophylline])
    assert [fold_result.structure, fold_result.energy] == ['.....((((.((.(((((...((((((((.....)))))...)))...))))).)).))))......',-26.920000076293945]
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['.....((((.((.(((((...((((((((.....)))))...)))...))))).)).))))......',-17.700000762939453]
    energy_model = None

    # example 2
    sequence = "CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    streptavidin_aptamer =              'CCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGG'
    streptavidin_aptamer_structure =    '((xxxxxxxxxxxxx((((xxxxxxxxx))))))'
    dG_streptavidin = -10.15
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence,aptamers=[streptavidin_aptamer],aptamer_constraints=[streptavidin_aptamer_structure],dG_ligands=[dG_streptavidin])
    assert [fold_result.structure, fold_result.energy] == ['((.............((((.........))))))................................',-10.550000190734863]
    fold_result  = energy_model.RNAfold(sequence)
    assert [fold_result.structure, fold_result.energy] == ['................((((((...((((...))))....))))))....................',-9.600000381469727]

    # example 3 - add multiple apatmers
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCCAGAAUCAUGCAAGUGCGUAAGAUAGUCGCGGGCUUUAAUUACGCUAUAUUAUAUACCCAAUUCU"
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAfold(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[0.0,0.0]) # what if the aptamer provides no additional free energy?
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))((((((((((((........))).))))).(((((((....))))..)))......))))........)))))).................',-31.399999618530273]
    fold_result  = energy_model.RNAfold(sequence,aptamers=[theophylline_aptamer,streptavidin_aptamer], \
        aptamer_constraints=[theophylline_constraint,streptavidin_aptamer_structure],dG_ligands=[dG_theophylline,dG_streptavidin]) # now with the aptamer binding free energies
    assert [fold_result.structure, fold_result.energy] == ['....((((((((((((((((......)))).))))))(((...)))(((........)))((((.............((((.........)))))))).......)))))).................',-34.95000076293945]
    print( 'Good.' )

    # Tests for RNAsubopt
    print( 'Testing RNAsubopt...' )
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    constraints = ["..xx........xx..","..xx........xx.."]
    delta_energy = 10
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    subopt_list = energy_model.RNAsubopt(sequences, constraints, delta_energy)
    assert len(subopt_list) == 10196
    energy_model = None
    sequence   = "CGCAGGGAUACCCGCG"
    delta_energy = 3
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    assert energy_model.RNAsubopt([sequence], delta_energy=3) == [('(((.(((...))))))', -5.599999904632568),
                                                                ('(((.((....)).)))', -5.5),
                                                                ('.((.(((...))))).', -4.699999809265137),
                                                                ('.((.((....)).)).', -4.599999904632568),
                                                                ('(((.((.....)))))', -3.5),
                                                                ('.((.((.....)))).', -2.5999999046325684),
                                                                ('....(((...)))...', -2.5999999046325684)]
    print( 'Good.' )

    print( 'Testing RNAsubopt with aptamers...' )
    sequence = "AGACAUAGCGAUCAAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAUAGCUAGCUACGAUGUAGCUCGGUAUUAUUU"
    theophylline_aptamer  = "GAUACCAG&CCCUUGGCAGC"
    theophylline_constraint = "(...((((&)...)))...)"
    dG_theophylline = -9.22
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    subopt_result = energy_model.RNAsubopt([sequence], delta_energy=10, aptamers=[theophylline_aptamer], aptamer_constraints=[theophylline_constraint], dG_ligands=[dG_theophylline])
    assert subopt_result[0].structure == '......(((.((.(((((......(((((.....))))).........))))).)).)))(((((.....)))))...........'
    assert subopt_result[0].energy == -30.719999313354492
    # Test without aptamer free energy to confirm different predictions
    # subopt_result = energy_model.RNAsubopt([sequence], delta_energy=10)
    # assert subopt_result[0].structure == '......(((.((.(((((...((((((((.....)))))...)))...))))).)).)))(((((.....)))))...........'
    # assert subopt_result[0].energy == -21.5
    energy_model = None
    print( 'Good.' )


    # Tests for vienna2bp and bp2vienna
    print( 'Testing vienna2bp and bp2vienna...' )
    energy_model  = PyVRNA()
    vienna_string = '..++++....++++....++++....++++....++++....++++....++++....++++..'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '....(((...)))...'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '(((.((....)).)))'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '..[.((....)).]..+++..+++..+++..+++..'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '(((((((....++..++..++..++..)))))))'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '[[[[(((((((....)))))))..((((....))))]]]]'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    energy_model  = PyVRNA()
    vienna_string = '...((..((..))..((..))..))....((..((..))..(.)..))..'
    assert energy_model.bptuple2vienna(energy_model.vienna2bp(vienna_string)) == vienna_string
    energy_model  = None
    print( 'Good.' )

def legacy_tests():
    """
    Tests each function in ViennaRNA class
    """

    print( "### === Testing ViennaRNA legacy class === ###" )

    # Tests for centroid
    print( 'Testing centroid...' )
    sequence      = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    energy_model  = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model  = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    cntrd_result  = energy_model.RNAcentroid(sequence)
    vienna_model.centroid(strands=[1], constraint=None, Temp=37.0, dangles="all", outputPS=False)
    assert [cntrd_result.structure, cntrd_result.energy] == [vienna_model['Centroid_bracket_string'], vienna_model["Centroid_energy"][0]] == ['(((.(.((.......))..))))....', 1.399999976158142]
    energy_model  = None
    vienna_model  = None
    sequence      = 'CGCAGGGAUACCCGCG' + 'GCGCCCAUAGGGACGC'
    energy_model  = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model  = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    cntrd_result  = energy_model.RNAcentroid(sequence)
    vienna_model.centroid(strands=[1], constraint=None, Temp=37.0, dangles="all", outputPS=False)
    assert [cntrd_result.structure, cntrd_result.energy] == [vienna_model['Centroid_bracket_string'], vienna_model["Centroid_energy"][0]] == ['(((.(((...))))))((((((...))).)))', -13.5]
    energy_model  = None
    vienna_model  = None
    print( 'Good.' )

    # Tests for energy
    print( 'Testing energy...' )
    sequence     = 'CGACGUAGAUGCUAGCUGACUCGAUGC'
    structure    = '(((.(.((.......))..))))....'
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_andronescu2007', Gquad=False)
    bp_tuple    = energy_model.vienna2bp(structure)
    vienna_model.energy(strands=[1], base_pairing_x=bp_tuple.bpx, base_pairing_y=bp_tuple.bpy, Temp=37.0, dangles="none")
    assert float("{0:0.2f}".format(energy_model.RNAeval([sequence], [structure]))) == float("{0:0.2f}".format(vienna_model["energy_energy"][0])) == 2.15
    energy_model = None
    vienna_model = None
    sequences    = ['CGCAGGGAUACCCGCG','GCGCCCAUAGGGACGC']
    structures  = ['.((.(((...))))).','((((((...))).)))']
    energy_model = PyVRNA()
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    bp_tuple    = energy_model.vienna2bp("".join(structures))
    vienna_model.energy(strands=[0, 1], base_pairing_x=bp_tuple.bpx, base_pairing_y=bp_tuple.bpy, Temp=37.0, dangles="all")
    assert float("{0:0.2f}".format(energy_model.RNAeval(sequences, structures))) == float("{0:0.2f}".format(vienna_model["energy_energy"][0])) == -12.72
    energy_model = None
    vienna_model = None
    print( 'Good.' )

    # Tests for mfe
    print( 'Testing mfe: RNAfold...' )
    sequence     = "CGCAGGGAUACCCGCG"
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_andronescu2007', Gquad=False)
    vienna_model.mfe(strands=[1], constraints=None, Temp=37.0, dangles="none", outputPS=False, duplex=False)
    assert list(energy_model.RNAfold(sequence)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['(((.(((...))))))', -4.619999885559082]
    energy_model = None
    vienna_model = None
    sequence     = "CGCAGGGAUACCCGCGGCGCCCAUAGGGACGC"
    constraints  = '.....xxx........................'
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    vienna_model.mfe(strands=[1], constraints=constraints, Temp=37.0, dangles="all", outputPS=False, duplex=False)
    assert list(energy_model.RNAfold(sequence, constraints)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['(((.(......).)))((((((...))).)))', -9.300000190734863]
    energy_model = None
    vienna_model = None
    print( 'Good.' )

    print( 'Testing mfe: RNAcofold...' )
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    energy_model = PyVRNA()
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    vienna_model.mfe(strands=[0, 1], constraints=None, Temp=37.0, dangles="all", outputPS=False, duplex=False)
    assert list(energy_model.RNAcofold(sequences)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['.((.(((...))))).&((((((...))).)))', -10.25 + energy_model.dG_init_adjustment]
    energy_model = None
    vienna_model = None
    sequences    = ["GCGCACAUAGUGACGC","GCGCCCAUAGGGACGC"]
    constraints  = ["....xx....xx....","....xx....xx...."]
    energy_model = PyVRNA(dangles=1, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    vienna_model.mfe(strands=[0, 1], constraints=constraints[0]+constraints[1], Temp=37.0, dangles="some", outputPS=False, duplex=False)
    assert list(energy_model.RNAcofold(sequences, constraints)) == [vienna_model["mfe_bracket_string"], vienna_model["mfe_energy"][0]] == ['((((............&))))............', -5.869999885559082 + energy_model.dG_init_adjustment]
    energy_model = None
    vienna_model = None
    print( 'Good.' )

    # Tests for RNAsubopt
    print( 'Testing subopt...' )
    sequences    = ["CGCAGGGAUACCCGCG", "GCGCCCAUAGGGACGC"]
    constraints = ["..xx........xx..", "..xx........xx.."]
    delta_energy = 10
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    vienna_model = ViennaRNA(Sequence_List=sequences, material='rna_andronescu2007', Gquad=False)
    vienna_model.subopt(strands=[0, 1], energy_gap=delta_energy, Temp=37.0, dangles="none", constraints="".join(constraints), outputPS=False)
    assert len(energy_model.RNAsubopt(sequences, constraints, delta_energy=delta_energy)) == len(vienna_model["subopt_energy"]) == 10196
    energy_model = None
    vienna_model = None
    sequence     = "CGCAGGGAUACCCGCG"
    delta_energy = 3
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='rna_turner2004.par')
    vienna_model = ViennaRNA(Sequence_List=[sequence], material='rna_turner2004', Gquad=False)
    vienna_model.subopt(strands=[1], energy_gap=delta_energy, Temp=37.0, dangles="all", constraints=None, outputPS=False)
    assert [x[1] for x in energy_model.RNAsubopt([sequence], delta_energy=delta_energy)] == vienna_model["subopt_energy"] == [-5.599999904632568,
                                                                                                                                  -5.5,
                                                                                                                                  -4.699999809265137,
                                                                                                                                  -4.599999904632568,
                                                                                                                                  -3.5,
                                                                                                                                  -2.5999999046325684,
                                                                                                                                  -2.5999999046325684]
    energy_model = None
    vienna_model = None
    print( 'Good.' )

    # Tests for convert_bracket_to_numbered_pairs and convert_numbered_pairs_to_bracket
    print( 'Testing convert_bracket_to_numbered_pairs ...' )
    print( 'Testing convert_bracket_to_numbered_pairs ...' )
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '................................................................'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '....(((...)))...'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '(((.((....)).)))'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '..[.((....)).]......................'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '(((((((....................)))))))'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '[[[[(((((((....)))))))..((((....))))]]]]'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    vienna_model  = ViennaRNA(Sequence_List=[], material='rna_turner2004', Gquad=False)
    vienna_string = '...((..((..))..((..))..))....((..((..))..(.)..))..'
    assert vienna_model.convert_numbered_pairs_to_bracket(*vienna_model.convert_bracket_to_numbered_pairs(vienna_string)) == vienna_string
    vienna_model  = None
    print( 'Good.' )

# ADDING SEPARATE TEST FUNCTIONS FOR NEW PYVRNA METHODS
# BE ADVISED: ASSERTIONS IN "tests()" WILL FAIL FOR NEWER VERSIONS OF VIENNARNA (e.g.v2.4.10)

def test_vienna2bplite():

    model = PyVRNA(temperature=37.0,dangles=0,gquad=True,pyindex=True,parameter_file="rna_turner1999.par")

    from time import time
    speedups = list()

    folds = [
        '.....................................',
        '((((....))))...(((...))).',
        '..((((((((((.((((((.(((((......))))).))))))..))).....((....))((((......)))))))))))...((((((......))))))...',
        '.((.(((...))))).&((((((...))).)))',
        '(((xxxx((((((((xxx(xx(xxxx((x((xxxxx))x))xxxx)xx)xx))))))))xxx)))',
        '(((....((((((((...(..(....((.((.....)).))....)..)..))))))))...)))',
        '(((((...(((...)))..((((.....))))...((((..(((..)))..))))..)))))...(((...)))...',
        '............&..((((....))))',
        '.((((((((..(((.(((((((((..(((.(((((.....))))).)))((((((((((..(((((((...(((((((.((((((((.((((.(((((((((((.((((.........)))))))))..((.....))......)))))))))).))))...((((((........((((..((((..((((...(((....)))...))))..))))........(((((....(((.....)))...)))))....(((((....(((((((((.....)).)))))))..)))))............(((((((.((.(((((....))))).)).)).))))))))).(((((((....((((((((((((.........))))))))).)))(((.(((((.(((.(((........................))).)))(((....(((..........)))....)))......((((....)))).......(((....)))..(((((.....(((((((((((.(((......))).))))).)))))).....)))))......))))))))....(((((((.(((((((((((((((((...))))))))).....(((.(((((.((((((...)))))).)))))..)))))))))))..)).))))))))))))))))))........))))...((((...........))))..)))))))....).)))))).((((((((((((...))))))))))))...((((.(((((...(((((................))))).)))))))))))))))........))))....))))))))).))).......))))))))(((..((((..........(((..(((((((....)))))))..)))((((....))))((....))..((((.((((((.((((((........(((((((((.(((((..((((....)))))))))....((((((.......))))))..............((((...((((.....))))...))))............)))))..........(((....))).......))))((((((((......)).))))))...((((.((.((((((.((((....))))...))))))))))))(((((.((((((.......))))))...)))))....))))))..)))))))))).....)))).)))..',
        '....(((.+..(((...+..(((...)))..+...)))..+.))).....',
    ]

    # Warning: structures with "[" or "]" get parsed as base pairs, not pseudoknots!
    # '....(((.+.[.(((...+..(((..].)))..+...)))..+.))).....'

    print('Testing vienna2bplite for validity and speed...')
    for vienna_string in folds:

        # check validity of result
        bptuple_new = model.vienna2bplite(vienna_string)
        bptuple = model.vienna2bp(vienna_string)

        assert bptuple_new.bpx == bptuple.bpx
        assert bptuple_new.bpy == bptuple.bpy
        assert bptuple_new.length == bptuple.length

        # timeit
        t0 = time()
        for _ in range(1000):
            bptuple_new = model.vienna2bplite(vienna_string)
        tf = time()
        dtnew = tf - t0

        t0 = time()
        for _ in range(1000):
            bptuple = model.vienna2bp(vienna_string)
        tf = time()
        dt = tf - t0

        speedups.append(dt/dtnew)

    print('Average fold-speedup using vienna2bplite: {:.2f}'.format(sum(speedups)/len(speedups)))
    print('Done.')

# WRITTEN & TESTED WITH ViennaRNA v2.4.10 - ACR 08/22/20
def test_RNAduplex():

    epsilon  = 10.0**-5

    # Tests for RNAduplex
    print( 'Testing RNAcofold...' )
    sequences    = ["CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC"]
    energy_model = PyVRNA(parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAduplex(sequences)
    assert fold_result.structure == '....(((.....(((.&))))))..........'
    assert abs(fold_result.energy - (-7.50 + energy_model.dG_init_adjustment)) < epsilon # where -7.50 is the command line value

    energy_model = None
    sequences    = ["AGCUAGCUUAGGCGCGAGGCGCGCGGGAGC","AAAAAAGGGGGGAAUAUAUCCCCCCAAAAA"]
    energy_model = PyVRNA(parameter_file='rna_turner2004.par')
    fold_result  = energy_model.RNAduplex(sequences)
    assert fold_result.structure == '........................((((..&..................))))........'
    assert (fold_result.energy - (-6.30 + energy_model.dG_init_adjustment)) < epsilon

    energy_model = None
    energy_model = PyVRNA(dangles=0, parameter_file='rna_andronescu2007.par')
    fold_result  = energy_model.RNAduplex(sequences)
    assert fold_result.structure == '........................((((..&..................))))........'
    assert (fold_result.energy - (-3.32 + energy_model.dG_init_adjustment)) < epsilon

    print( 'Good.' )

if __name__ == '__main__':
    # Tests for RNAsubopt
    print( 'Testing RNAsubopt...' )
    sequences    = ["GGGUGCCGGGGC","ACCTCCTTA"]
    constraints = []
    delta_energy = 8.0
    energy_model = PyVRNA(dangles=0, gquad=False, parameter_file='rna_andronescu2007.par')
    subopt_list = energy_model.RNAsubopt(sequences=sequences, constraints=[], delta_energy=delta_energy)
    print(subopt_list)

    tests()
    # legacy_tests()
    # test_vienna2bplite()
    # test_RNAduplex()
