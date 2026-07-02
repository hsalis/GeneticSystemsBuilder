import os, random, sys, time, traceback, itertools, math, copy, json, string, shutil, re, tempfile, threading, requests, pickle, ast, argparse
import multiprocessing
from datetime import datetime as dt
import numpy as np
import pandas as pd
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from numba import njit, prange

try:
    sys.path.append('NRP_Calculator/')
    import NRP_Calculator.nrpcalc as nrpcalc
except ImportError:
    nrpcalc = None
    raise ImportError("Could not import nrpcalc. The Non-Repetitive Parts Calculator is used to design PCR primers.")

try:
    from PyVRNA import PyVRNA
except ImportError:
    PyVRNA = None
    raise ImportError("Could not import PyVRNA. PyVRNA is used to design PCR primers.")

def random_id(length = 6):
    return "".join([random.choice(string.ascii_uppercase) for n in range(length)])

def _revcomp(seq):
    comp = {'A' : 'T', 'G' : 'C', 'T' : 'A', 'C' : 'G'}
    x = "".join([comp[letter] for letter in seq ])
    return x[::-1]  #reverse

restriction_enzyme_sites = { 'AclI': 'AACGTT', 'HindIII': 'AAGCTT', 'SspI': 'AATATT', 'MluCI': 'AATT', 'Tsp509I': 'AATT', 'PciI': 'ACATGT', 'AgeI': 'ACCGGT',
                             'BspMI': 'ACCTGC', 'BfuAI': 'ACCTGC', 'SexAI': 'ACC[AT]GGT', 'MluI': 'ACGCGT', 'MluI-HF': 'ACGCGT', 'BceAI': 'ACGGC', 'HpyCH4IV': 'ACGT', 'HpyCH4III': 'AC[ATCG]{1}GT',
                             'BaeI': 'AC[ATCG]{4}GTA[CT]C', 'BsaXI': 'AC[ATCG]{5}CTCC', 'AflIII': 'AC[AG][CT]GT', 'SpeI': 'ACTAGT', 'BsrI': 'ACTGG', 'BmrI': 'ACTGGG', 'BglII': 'AGATCT',
                             'AfeI': 'AGCGCT', 'AluI': 'AGCT', 'StuI': 'AGGCCT', 'ScaI': 'AGTACT', 'ClaI': 'ATCGAT', 'BspDI': 'ATCGAT', 'PI-SceI': 'ATCTATGTCGGGTGCGGAGAAAGAGGTAAT',
                             'NsiI': 'ATGCAT', 'NsiI-HF': 'ATGCAT', 'AseI': 'ATTAAT', 'SwaI': 'ATTTAAAT', 'CspCI': 'CAA[ATCG]{5}GTGG', 'MfeI': 'CAATTG', 'BssSI': 'CACGAG', 'BssS_alpha_I': 'CACGAG',
                             'BmgBI': 'CACGTC', 'PmlI': 'CACGTG', 'DraIII': 'CAC[ATCG]{3}GTG', 'AleI': 'CAC[ATCG]{4}GTG', 'EcoP15I': 'CAGCAG', 'PvuII': 'CAGCTG', 'AlwNI': 'CAG[ATCG]{3}CTG',
                             'BtsIMutI': 'CAGTG', 'TspRI': 'CA[CG]TG', 'NdeI': 'CATATG', 'NlaIII': 'CATG', 'CviAII': 'CATG', 'FatI': 'CATG', 'MslI': 'CA[CT][ATCG]{4}[AG]TG', 'FspEI': 'CC',
                             'XcmI': 'CCA[ATCG]{9}TGG', 'BstXI': 'CCA[ATCG]{6}TGG', 'PflMI': 'CCA[ATCG]{5}TGG', 'BccI': 'CCATC', 'NcoI': 'CCATGG', 'BseYI': 'CCCAGC', 'FauI': 'CCCGC',
                             'SmaI': 'CCCGGG', 'XmaI': 'CCCGGG', 'TspMI': 'CCCGGG', 'Nt.CviPII': 'CC[AGT]', 'LpnPI': 'CC[AGT]G', 'AciI': 'CCGC', 'SacII': 'CCGCGG', 'BsrBI': 'CCGCTC',
                             'MspI': 'CCGG', 'HpaII': 'CCGG', 'ScrFI': 'CC[ATCG]{1}GG', 'BssKI': 'CC[ATCG]{1}GG', 'StyD4I': 'CC[ATCG]{1}GG', 'BsaJI': 'CC[ATCG]{2}GG', 'BslI': 'CC[ATCG]{7}GG',
                             'BtgI': 'CC[AG][CT]GG', 'NciI': 'CC[CG]GG', 'AvrII': 'CCTAGG', 'MnlI': 'CCTC', 'BbvCI': 'CCTCAGC', 'Nb.BbvCI': 'CCTCAGC', 'Nt.BbvCI': 'CCTCAGC', 'SbfI': 'CCTGCAGG',
                             'Bpu10I': 'CCT[ATCG]{1}AGC', 'Bsu36I': 'CCT[ATCG]{1}AGG', 'EcoNI': 'CCT[ATCG]{5}AGG', 'HpyAV': 'CCTTC', 'BstNI': 'CC[AT]GG', 'PspGI': 'CC[AT]GG', 'StyI': 'CC[AT][AT]GG',
                             'BcgI': 'CGA[ATCG]{6}TGC', 'PvuI': 'CGATCG', 'BstUI': 'CGCG', 'EagI': 'CGGCCG', 'RsrII': 'CGG[AT]CCG', 'BsiEI': 'CG[AG][CT]CG', 'BsiWI': 'CGTACG', 'BsmBI': 'CGTCTC',
                             'Hpy99I': 'CG[AT]CG', 'MspA1I': 'C[AC]GC[GT]G', 'MspJI': 'C[ATCG]{2}[AG]', 'SgrAI': 'C[AG]CCGG[CT]G', 'BfaI': 'CTAG', 'BspCNI': 'CTCAG', 'XhoI': 'CTCGAG', 'PaeR7I': 'CTCGAG',
                             'TliI': 'CTCGAG', 'EarI': 'CTCTTC', 'AcuI': 'CTGAAG', 'PstI': 'CTGCAG', 'BpmI': 'CTGGAG', 'DdeI': 'CT[ATCG]{1}AG', 'SfcI': 'CT[AG][CT]AG', 'AflII': 'CTTAAG', 'BpuEI': 'CTTGAG',
                             'SmlI': 'CT[CT][AG]AG', 'AvaI': 'C[CT]CG[AG]G', 'BsoBI': 'C[CT]CG[AG]G', 'MboII': 'GAAGA', 'BbsI': 'GAAGAC', 'XmnI': 'GAA[ATCG]{4}TTC', 'BsmI': 'GAATGC', 'Nb.BsmI': 'GAATGC',
                             'EcoRI': 'GAATTC', 'HgaI': 'GACGC', 'AatII': 'GACGTC', 'ZraI': 'GACGTC', 'Tth111I': 'GAC[ATCG]{3}GTC', 'PflFI': 'GAC[ATCG]{3}GTC', 'PshAI': 'GAC[ATCG]{4}GTC',
                             'AhdI': 'GAC[ATCG]{5}GTC', 'DrdI': 'GAC[ATCG]{6}GTC', 'Eco53kI': 'GAGCTC', 'SacI': 'GAGCTC', 'BseRI': 'GAGGAG', 'PleI': 'GAGTC', 'Nt.BstNBI': 'GAGTC', 'MlyI': 'GAGTC',
                             'HinfI': 'GA[ATCG]{1}TC', 'EcoRV': 'GATATC', 'MboI': 'GATC', 'DpnII': 'GATC', 'Sau3AI': 'GATC', 'BfuCI': 'GATC', 'DpnI': 'GATC', 'BsaBI': 'GAT[ATCG]{4}ATC',
                             'TfiI': 'GA[AT]TC', 'BsrDI': 'GCAATG', 'Nb.BsrDI': 'GCAATG', 'BbvI': 'GCAGC', 'BtsI': 'GCAGTG', 'Bts_alpha_I': 'GCAGTG', 'Nb.BtsI': 'GCAGTG', 'BstAPI': 'GCA[ATCG]{5}TGC',
                             'SfaNI': 'GCATC', 'SphI': 'GCATGC', 'NmeAIII': 'GCCGAG', 'NaeI': 'GCCGGC', 'NgoMIV': 'GCCGGC', 'BglI': 'GCC[ATCG]{5}GGC', 'AsiSI': 'GCGATCGC', 'BtgZI': 'GCGATG',
                             'HinP1I': 'GCGC', 'HhaI': 'GCGC', 'BssHII': 'GCGCGC', 'NotI': 'GCGGCCGC', 'Fnu4HI': 'GC[ATCG]{1}GC', 'Cac8I': 'GC[ATCG]{2}GC', 'MwoI': 'GC[ATCG]{7}GC', 'NheI': 'GCTAGC',
                             'BmtI': 'GCTAGC', 'SapI': 'GCTCTTC', 'BspQI': 'GCTCTTC', 'Nt.BspQI': 'GCTCTTC', 'BlpI': 'GCT[ATCG]{1}AGC', 'TseI': 'GC[AT]GC', 'ApeKI': 'GC[AT]GC', 'Bsp1286I': 'G[AGT]GC[ACT]C',
                             'AlwI': 'GGATC', 'Nt.AlwI': 'GGATC', 'BamHI': 'GGATCC', 'FokI': 'GGATG', 'BtsCI': 'GGATG', 'HaeIII': 'GGCC', 'PhoI': 'GGCC', 'FseI': 'GGCCGGCC', 'SfiI': 'GGCC[ATCG]{5}GGCC',
                             'NarI': 'GGCGCC', 'KasI': 'GGCGCC', 'SfoI': 'GGCGCC', 'PluTI': 'GGCGCC', 'AscI': 'GGCGCGCC', 'EciI': 'GGCGGA', 'BsmFI': 'GGGAC', 'ApaI': 'GGGCCC', 'PspOMI': 'GGGCCC',
                             'Sau96I': 'GG[ATCG]{1}CC', 'NlaIV': 'GG[ATCG]{2}CC', 'KpnI': 'GGTACC', 'Acc65I': 'GGTACC', 'BsaI': 'GGTCTC', 'HphI': 'GGTGA', 'BstEII': 'GGT[ATCG]{1}ACC', 'AvaII': 'GG[AT]CC',
                             'BanI': 'GG[CT][AG]CC', 'BaeGI': 'G[GT]GC[AC]C', 'BsaHI': 'G[AG]CG[CT]C', 'BanII': 'G[AG]GC[CT]C', 'RsaI': 'GTAC', 'CviQI': 'GTAC', 'BstZ17I': 'GTATAC', 'BciVI': 'GTATCC',
                             'SalI': 'GTCGAC', 'Nt.BsmAI': 'GTCTC', 'BsmAI': 'GTCTC', 'BcoDI': 'GTCTC', 'ApaLI': 'GTGCAC', 'BsgI': 'GTGCAG', 'AccI': 'GT[AC][GT]AC', 'Hpy166II': 'GT[ATCG]{2}AC',
                             'Tsp45I': 'GT[CG]AC', 'HpaI': 'GTTAAC', 'PmeI': 'GTTTAAAC', 'HincII': 'GT[CT][AG]AC', 'BsiHKAI': 'G[AT]GC[AT]C', 'ApoI': '[AG]AATT[CT]', 'NspI': '[AG]CATG[CT]',
                             'BsrFI': '[AG]CCGG[CT]', 'BstYI': '[AG]GATC[CT]', 'HaeII': '[AG]GCGC[CT]', 'CviKI-1': '[AG]GC[CT]', 'EcoO109I': '[AG]GG[ATCG]{1}CC[CT]', 'PpuMI': '[AG]GG[AT]CC[CT]',
                             'I-CeuI': 'TAACTATAACGGTCCTAAGGTAGCGAA', 'SnaBI': 'TACGTA', 'I-SceI': 'TAGGGATAACAGGGTAAT', 'BspHI': 'TCATGA', 'BspEI': 'TCCGGA', 'MmeI': 'TCC[AG]AC', 'Taq_alpha_I': 'TCGA',
                             'NruI': 'TCGCGA', 'NruI-HF': 'TCGCGA', 'Hpy188I': 'TC[ATCG]{1}GA', 'Hpy188III': 'TC[ATCG]{2}GA', 'XbaI': 'TCTAGA', 'BclI': 'TGATCA', 'HpyCH4V': 'TGCA', 'FspI': 'TGCGCA',
                             'PI-PspI': 'TGGCAAACAGCTATTATGGGTATTATGGGT', 'MscI': 'TGGCCA', 'BsrGI': 'TGTACA', 'BsrGI-HF': 'TGTACA', 'MseI': 'TTAA', 'PacI': 'TTAATTAA', 'PsiI': 'TTATAA', 'BstBI': 'TTCGAA',
                             'DraI': 'TTTAAA', 'PspXI': '[ACG]CTCGAG[CTG]', 'BsaWI': '[AT]CCGG[AT]', 'BsaAI': '[CT]ACGT[AG]', 'EaeI': '[CT]GGCC[AG]',
                             'PaqCI': 'CACCTGC', 'AarI': 'CACCTGC', 'SrfI': 'GCCCGGGC', 'Esp3I' : 'CGTCTC',
                            }

SPLIT_SOLUTION_BY_COMBO_SET: dict = {}   # key -> {'positions': List[int], 'meta': {...}}

try:
    current_directory = os.path.dirname(__file__)
except:
    current_directory = os.getcwd()

class ScreeningRequestThread(threading.Thread):
    def __init__(self, base_url, data):
        threading.Thread.__init__(self)
        self.base_url = base_url
        self.data = data
        self.response = None
        self.status_code = None

    def run(self):
        try:
            response = requests.post(self.base_url + "/v1/screen", data=json.dumps(self.data), headers={"Content-Type": "application/json"})
            self.response = response.json()
            self.status_code = response.status_code
        except Exception as e:
            print(traceback.format_exc())
            self.response = str(e)
            self.status_code = None

    def get_response(self):
        return (self.response, self.status_code)

def createLongAssembly(data_inputs, pool):

    ## assembly_levels == 2 only for now [will use recursion later for levels 3, 4 ...]
    (oligoPoolSpecification, sequence, levelOneEnzyme, levelTwoEnzyme, circular, assembly_id, assembly_name, combinatorial_set, no_overhangs_position_range_list, circular_anchor, verbose) = data_inputs
    landingPadOverhangs =  ('TAAA', 'TTAG')

    maximum_oligo_length = oligoPoolSpecification['maximum_oligo_length']

    ## Level Two Assembly (BbsI)
    levelTwoCircular = circular
    (Success, assembly, error_message) = createAssembly( (oligoPoolSpecification, sequence, oligoPoolSpecification['LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH'], oligoPoolSpecification['LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH'], levelTwoEnzyme, levelTwoCircular, assembly_id, assembly_name, combinatorial_set, no_overhangs_position_range_list, circular_anchor, verbose) )
    if not Success:
        return (False, [assembly], [error_message])

    numFragmentsLongAssembly = len(assembly['fragment_list'])

    levelTwoFragmentList = []
    for (n, levelTwoFragment) in enumerate(assembly['fragment_list']):
        fragmentSequence = levelTwoFragment['fragment_sequence']
        fragmentWithCutSites = multiLevelFragmentConstructor(oligoPoolSpecification, fragmentSequence, levelTwoEnzyme, landingPadOverhangs)

        if no_overhangs_position_range_list is None:
            no_overhangs_ranges_in_fragment = None
        else:
            # --- compute forbidden ranges in no_overhangs_position_range_list mapped into fragment-local coordinates ---
            seqlen   = len(assembly['sequence'])
            frag_len = len(fragmentSequence)
            fs       = levelTwoFragment['begin'] % seqlen                    # fragment start on original
            fe_ext   = fs + frag_len - 1                                     # end (possibly past seqlen-1)
            # fragment covers one or two linear intervals on the original sequence:
            frag_intervals = [(fs, min(seqlen - 1, fe_ext))] + ([(0, fe_ext % seqlen)] if fe_ext >= seqlen else [])
            # intersect each excluded [a,b] with the fragment interval(s) and map into [0..frag_len-1]
            no_overhangs_ranges_in_fragment = [
                [ (max(a, u) - fs) % seqlen, (min(b, v) - fs) % seqlen ]
                for (a, b) in no_overhangs_position_range_list
                for (u, v) in frag_intervals
                if max(a, u) <= min(b, v)
            ]

        fragment = {'name' : f"{assembly_name};L2-frag-{n}",
                    'combinatorial_set' : combinatorial_set,
                    'assembly_id' : assembly_id + float(n+1)/100.0,
                    'nucleotide_sequence' : fragmentWithCutSites,
                    'no_overhangs_position_range_list' : no_overhangs_ranges_in_fragment,
                    'circular_anchor' : None,
                   }
        levelTwoFragmentList.append(fragment)

    if verbose: print(f'levelTwoFragmentList: {levelTwoFragmentList}')
    ## Level One Assembly
    primer_sequence_length = oligoPoolSpecification.get('primer_sequence_length', oligoPoolSpecification['DEFAULT_PRIMER_SEQUENCE_LENGTH'])
    overhang_length = oligoPoolSpecification['OVERHANG_LENGTH'][levelOneEnzyme]
    flanking_sequence_length = len(landingPadOverhangs[0]+landingPadOverhangs[1])+len(oligoPoolSpecification['CUT_SITE_SEQUENCES'][levelOneEnzyme][0]) + len(oligoPoolSpecification['CUT_SITE_SEQUENCES'][levelOneEnzyme][1]) + 2 * primer_sequence_length + overhang_length
    maximum_fragment_length = maximum_oligo_length - flanking_sequence_length

    levelOneCircular = False
    inputList = [(oligoPoolSpecification, fragment['nucleotide_sequence'], oligoPoolSpecification['MINIMUM_FRAGMENT_LENGTH'], maximum_fragment_length, levelOneEnzyme, levelOneCircular,
                  fragment['assembly_id'], fragment['name'], fragment['combinatorial_set'], fragment['no_overhangs_position_range_list'], fragment['circular_anchor'], verbose) for fragment in levelTwoFragmentList]

    if pool is not None:
        longAssemblyList = pool.map(createAssembly, inputList)
        pool.close()
        pool.join()
    else:
        longAssemblyList = list(map(createAssembly, inputList))

    outputAssemblyList = []
    if all([success for (success, assembly, error_message) in longAssemblyList]):
        for (n, (success, assembly, error_message)) in enumerate(longAssemblyList):
            assembly['long_assembly_fragment_number'] = n
            outputAssemblyList.append(assembly)
    else:
        errorAssemblyList = []
        error_messages = []
        for (n, (success, assembly, error_message)) in enumerate(longAssemblyList):
            if not success:
                errorAssemblyList.append(assembly)
                error_messages.append(error_message)
        return (False, errorAssemblyList, error_messages)

    return (True, outputAssemblyList, None)

def getOverhangPositions(oligoPoolSpecification, sequence, enzyme, overhang_length, maximum_fragment_length,
                         only_use_overhangs = None,
                         no_overhangs_position_range_list = None,
                         MIN_AVG_FIDELITY = None, position_step = None,
                         ):

    seqlen = len(sequence)

    if position_step is None:
        if seqlen < 2000:
            position_step = 1
        elif seqlen < 4000:
            position_step = 2
        elif seqlen < 6000:
            position_step = 3
        elif seqlen < 8000:
            position_step = 4
        elif seqlen < 10000:
            position_step = 5
        else:
            position_step = 6

    if MIN_AVG_FIDELITY is None:
        minNumFragments = int(seqlen / maximum_fragment_length) + 1
        if minNumFragments < 4:
            MIN_AVG_FIDELITY = 0.85
        elif minNumFragments < 8:
            MIN_AVG_FIDELITY = 0.80
        elif minNumFragments < 12:
            MIN_AVG_FIDELITY = 0.70
        elif minNumFragments < 16:
            MIN_AVG_FIDELITY = 0.65
        else:
            MIN_AVG_FIDELITY = 0.60

    if only_use_overhangs is None:
        availableOverhangs = [overhang for (overhang, avg_fidelity) in oligoPoolSpecification['RANKED_OVERHANGS'][enzyme].items() if avg_fidelity >= MIN_AVG_FIDELITY]
    else:
        availableOverhangs = only_use_overhangs

    allowableOverhangPositions = {}
    shifted_p1nt_overhang_positions = {}
    shifted_n1nt_overhang_positions = {}

    if position_step > 1:
        position_list = list(range(0, seqlen-overhang_length+1, position_step)) + [seqlen-overhang_length]
    else:
        position_list = list(range(0, seqlen-overhang_length+1))

    # Remove positions found in no_overhangs_position_range_list
    if no_overhangs_position_range_list is not None:
        position_list = [p for p in position_list if (p == 0) or (p == seqlen-overhang_length) or all(not (a <= p <= b) for a, b in no_overhangs_position_range_list)]

    #print(f'position_list = {position_list}')

    for position in position_list:
        candidate_overhang = sequence[position:position + overhang_length]
        if position < seqlen - overhang_length:
            overhang_p1nt_shifted = sequence[position + 1:position + overhang_length + 1]
            shifted_p1nt_overhang_positions[position] = overhang_p1nt_shifted
        else:
            shifted_p1nt_overhang_positions[position] = None
        if position > 0:
            overhang_n1nt_shifted = sequence[position - 1:position + overhang_length - 1]
            shifted_n1nt_overhang_positions[position] = overhang_n1nt_shifted
        else:
            shifted_n1nt_overhang_positions[position] = None

        if candidate_overhang in availableOverhangs or (position == 0) or (position == seqlen-overhang_length): allowableOverhangPositions[position] = candidate_overhang

    def calculateWeights(selected_positions, candidate_positions):
        weights = []
        selected_overhang_set = [ (allowableOverhangPositions[position], position) for position in selected_positions]
        for (n, candidate_position) in enumerate(candidate_positions):
            candidate_overhang_set = selected_overhang_set + [(allowableOverhangPositions[candidate_position], candidate_position)]
            weights.append(calculateLigationFidelity(oligoPoolSpecification, candidate_overhang_set, enzyme, shifted_p1nt_overhang_positions, shifted_n1nt_overhang_positions))
        return weights

    lig_matrix, rev_idx, allow_idx, shift_p, shift_n, overhang_to_idx = prepare_numba_ligation_data(
        oligoPoolSpecification=oligoPoolSpecification,
        enzyme=enzyme,
        allowableOverhangPositions=allowableOverhangPositions,
        shifted_p1nt_overhang_positions=shifted_p1nt_overhang_positions,
        shifted_n1nt_overhang_positions=shifted_n1nt_overhang_positions,
        revcomp_func=_revcomp,
    )

    def calculateWeights_fast_numba(
        selected_positions,
        candidate_positions,
        SHIFT_PROBABILITY=0.001,
    ):
        # Build base set arrays in the same order
        base_pos = np.asarray(selected_positions, dtype=np.int32)
        base_overhang_idx = np.empty(base_pos.size, dtype=np.int32)

        for i, pos in enumerate(base_pos):
            idx = allow_idx[pos]
            if idx < 0:
                raise KeyError(f"Position {pos} not found in allowableOverhangPositions / allow_idx mapping.")
            base_overhang_idx[i] = idx

        cand_pos = np.asarray(candidate_positions, dtype=np.int32)

        weights = calculateWeights_numba(
            base_overhang_idx, base_pos,
            cand_pos,
            allow_idx,
            lig_matrix, rev_idx,
            shift_p, shift_n,
            SHIFT_PROBABILITY
        )

        # return python list
        return weights.tolist()

    return (allowableOverhangPositions, shifted_p1nt_overhang_positions, shifted_n1nt_overhang_positions, calculateWeights_fast_numba)

def createAssembly(data_inputs):

    (oligoPoolSpecification, sequence, minimum_fragment_length, maximum_fragment_length, enzyme, circular, assembly_id, assembly_name, combinatorial_set, no_overhangs_position_range_list, circular_anchor, verbose) = data_inputs

    t1 = time.time()

    overhang_length = oligoPoolSpecification['OVERHANG_LENGTH'][enzyme]
    (allowableOverhangPositions, shifted_p1nt_overhang_positions, shifted_n1nt_overhang_positions, calculateWeights) = getOverhangPositions(oligoPoolSpecification, sequence, enzyme, overhang_length, maximum_fragment_length, no_overhangs_position_range_list=no_overhangs_position_range_list)
    if circular_anchor is not None:
        anchors = [circular_anchor]
    else:
        anchors = None

    seqlen = len(sequence)

    t1 = time.time()
    from GeneticSystemsCalculator_Build_Helpers import split_line_adaptive_weights, split_circle_adaptive_weights

    # ---------- Try cache first ----------

    # Include constraints in the key so we don't reuse incompatible splits.
    cache_key = (combinatorial_set, enzyme, circular, overhang_length, minimum_fragment_length, maximum_fragment_length, seqlen)
    S = None
    cached = SPLIT_SOLUTION_BY_COMBO_SET.get(cache_key)

    def evaluate_initial(cached_positions, allowableInitialOverhangPositions):
        """Evaluate a cached split on the current sequence via DP's initial_solution mode."""
        # keep only positions present for this sequence
        internal_positions = [p for p in cached_positions if p in allowableInitialOverhangPositions]
        # Build initial_solution as list[(position, overhang_4mer)]
        initial_solution = [(p, allowableInitialOverhangPositions[p]) for p in internal_positions]
        # Call the DP with initial_solution; if the circular version doesn't support it yet, fall back.
        if circular:
            try:
                return split_circle_adaptive_weights(
                    positions=list(allowableInitialOverhangPositions.keys()),
                    sequences=list(allowableInitialOverhangPositions.values()),
                    L_total=seqlen-overhang_length,
                    L_min=minimum_fragment_length,
                    L_max=maximum_fragment_length,
                    calculateWeights=calculateWeights,
                    include_last_edge_cost=True,
                    beam_size=2,
                    beam_by_node=True,
                    initial_solution=initial_solution,   # input cached initial solution
                )
            except TypeError:
                pass  # helper lacks initial_solution; fall through to solve fresh
        else:
            try:
                return split_line_adaptive_weights(
                    positions=list(allowableInitialOverhangPositions.keys()),
                    sequences=list(allowableInitialOverhangPositions.values()),
                    L_total=seqlen-overhang_length,
                    L_min=minimum_fragment_length,
                    L_max=maximum_fragment_length,
                    calculateWeights=calculateWeights,
                    include_last_edge_cost=False,
                    beam_size=2,
                    beam_by_node=True,
                    state_compressor=None,
                    initial_solution=initial_solution,   # input cached initial solution
                )
            except TypeError:
                pass
        return None

    if cached:
        S = evaluate_initial(cached.get('positions', []), allowableOverhangPositions)
        if S is not None and S.reached:
            print('Reused cached split for combinatorial_set:', combinatorial_set, '->', S)
        else:
            print('Evaluated cached split for combinatorial_set: ', combinatorial_set, ' but it was not feasible ->', S)

    # ---------- If no usable cache, solve as before ----------
    if S is None or not S.reached:
        if circular:
            S = split_circle_adaptive_weights(
                positions=list(allowableOverhangPositions.keys()),
                sequences=list(allowableOverhangPositions.values()),
                L_total=seqlen-overhang_length,
                L_min=minimum_fragment_length,
                L_max=maximum_fragment_length,
                calculateWeights=calculateWeights,
                include_last_edge_cost=True,
                beam_size=2,
                beam_by_node=True,
                anchors = anchors,
            )
            if S.reached and verbose: print('split_circle_adaptive_weights solution: ', S)
        else:
            S = split_line_adaptive_weights(
                positions=list(allowableOverhangPositions.keys()),
                sequences=list(allowableOverhangPositions.values()),
                L_total=seqlen-overhang_length,
                L_min=minimum_fragment_length,
                L_max=maximum_fragment_length,
                calculateWeights=calculateWeights,
                include_last_edge_cost=False,
                beam_size=2,
                beam_by_node=True,
                state_compressor=None,
            )
            if S.reached and verbose: print('split_line_adaptive_weights solution: ', S)

    # ---------- Fallback second try (your existing logic) ----------
    if (not S.reached) or (S.sum_dynamic_cost > 0.01):
        (allowableOverhangPositions, shifted_p1nt_overhang_positions, shifted_n1nt_overhang_positions, calculateWeights) = getOverhangPositions(oligoPoolSpecification, sequence, enzyme, overhang_length, maximum_fragment_length, no_overhangs_position_range_list=no_overhangs_position_range_list, MIN_AVG_FIDELITY=0.60, position_step=1)
        minimum_fragment_length = 100
        cache_key = (combinatorial_set, enzyme, circular, overhang_length, minimum_fragment_length, maximum_fragment_length, seqlen)
        cached = SPLIT_SOLUTION_BY_COMBO_SET.get(cache_key)
        S = None

        if cached:
            S = evaluate_initial(cached.get('positions', []), allowableOverhangPositions)
            if S and S.reached and verbose:
                print('Reused cached split for combinatorial_set:', combinatorial_set, '->', S)

        if not S or not S.reached:
            if circular:
                S = split_circle_adaptive_weights(
                    positions=list(allowableOverhangPositions.keys()),
                    sequences=list(allowableOverhangPositions.values()),
                    L_total=seqlen-overhang_length,
                    L_min=minimum_fragment_length,
                    L_max=maximum_fragment_length,
                    calculateWeights=calculateWeights,
                    include_last_edge_cost=True,
                    beam_size=2,
                    beam_by_node=True,
                    anchors = anchors,
                )
                if verbose: print('2nd try split_circle_adaptive_weights solution: ', S)
            else:
                S = split_line_adaptive_weights(
                    positions=list(allowableOverhangPositions.keys()),
                    sequences=list(allowableOverhangPositions.values()),
                    L_total=seqlen-overhang_length,
                    L_min=minimum_fragment_length,
                    L_max=maximum_fragment_length,
                    calculateWeights=calculateWeights,
                    include_last_edge_cost=False,
                    beam_size=2,
                    beam_by_node=True,
                    state_compressor=None,
                )
                if verbose:
                    print(f'Could not find solution on 2nd try [split_line_adaptive_weights solution {S}]')
                    print(f'sequence: {sequence}')

        if not S.reached:
            assembly = {'id': assembly_id, 'name': assembly_name, 'sequence': sequence, 'is_circular': circular,
                        'assembly_enzyme': enzyme, 'combinatorial_set': combinatorial_set}
            return (False, assembly, f'Assembly {assembly_name} [{assembly_id}]: No build solution found!')

    # ---------- Update cache with this solution (store internal cut positions only) ----------
    # (We cache only if it truly reached.)
    if S.reached:
        # positions list includes 0 and seqlen; cache only internal cuts
        internal_positions = [p for p in S.positions if p not in (0, seqlen)]
        SPLIT_SOLUTION_BY_COMBO_SET[cache_key] = {
            'positions': internal_positions,
            'meta': {
                'ligation_fidelity': 1.0 - S.sum_dynamic_cost,
                'N': S.N,
            }
        }

    # ---------- Build assembly ----------
    ligation_fidelity = 1.0 - S.sum_dynamic_cost
    numFragments = S.N
    fragmentList = []
    for n in range(numFragments):
        overhang_5p = S.segments[n].begin_sequence
        overhang_3p = S.segments[n].end_sequence
        begin_pos = S.segments[n].start
        end_pos = S.segments[n].end
        if circular and begin_pos == end_pos and S.segments[n].length == seqlen - overhang_length:
            fragment_sequence = sequence[begin_pos: seqlen] + sequence[0: begin_pos + overhang_length]
        elif begin_pos > end_pos:  # circular fragment
            fragment_sequence = sequence[begin_pos: seqlen] + sequence[0: end_pos + overhang_length]
        else:
            fragment_sequence = sequence[begin_pos: end_pos + overhang_length]

        fragment_with_cut_sites = fragmentConstructor(oligoPoolSpecification, fragment_sequence, enzyme)
        item = {'fragment_sequence': fragment_sequence, 'fragment_with_cut_sites': fragment_with_cut_sites,
                'begin': begin_pos, 'end': end_pos + overhang_length,
                'overhang_5p': overhang_5p, 'overhang_3p': overhang_3p}
        #print(f"Fragment info: {item}")
        fragmentList.append(item)

    if circular:
        first_overhang = None
        last_overhang = None
    else:
        first_overhang = S.segments[0].begin_sequence
        last_overhang = S.segments[-1].end_sequence

    assembly = {'id': assembly_id, 'name': assembly_name, 'sequence': sequence, 'is_circular': circular,
                'assembly_enzyme': enzyme, 'combinatorial_set': combinatorial_set,
                'first_overhang': first_overhang, 'last_overhang': last_overhang,
                'fragment_list': fragmentList, 'ligation_fidelity': ligation_fidelity}

    t2 = time.time()
    print('Completed Assembly {} [{}]: Enzyme {}. Sequence Length: {} nt. {} oligos. '
          'Predicted ligation fidelity: {}. Runtime: {} seconds.'
          .format(assembly['name'], assembly['id'], assembly['assembly_enzyme'],
                  len(assembly['sequence']), len(assembly['fragment_list']),
                  assembly['ligation_fidelity'], round(t2 - t1, 3)))
    return (True, assembly, None)

def calculateLigationFidelity(oligoPoolSpecification, overhangSet, enzyme, shifted_p1nt_overhang_positions, shifted_n1nt_overhang_positions):

    ligationFrequencies = oligoPoolSpecification['ALL_LIGATION_COUNTS'][enzyme]
    SHIFT_PROBABILITY = 0.001

    ligation_fidelity = 1.0
    ligation_fidelity_revcomp = 1.0
    for (overhang, position) in overhangSet:
        total_correct = ligationFrequencies[overhang][overhang]
        total_correct_revcomp = ligationFrequencies[_revcomp(overhang)][_revcomp(overhang)]

        total_incorrect = ligationFrequencies[overhang][_revcomp(overhang)]
        total_incorrect_revcomp = ligationFrequencies[_revcomp(overhang)][overhang]

        if shifted_p1nt_overhang_positions[position] is not None:
            total_incorrect_shifted_p1nt = ligationFrequencies[overhang][_revcomp(shifted_p1nt_overhang_positions[position])]
            total_incorrect_shifted_p1nt_revcomp = ligationFrequencies[_revcomp(overhang)][shifted_p1nt_overhang_positions[position]]
        else:
            total_incorrect_shifted_p1nt = 0
            total_incorrect_shifted_p1nt_revcomp = 0

        if shifted_n1nt_overhang_positions[position] is not None:
            total_incorrect_shifted_n1nt = ligationFrequencies[overhang][_revcomp(shifted_n1nt_overhang_positions[position])]
            total_incorrect_shifted_n1nt_revcomp = ligationFrequencies[_revcomp(overhang)][shifted_n1nt_overhang_positions[position]]
        else:
            total_incorrect_shifted_n1nt = 0
            total_incorrect_shifted_n1nt_revcomp = 0

        for (other_overhang, other_position) in overhangSet:
            if other_overhang != overhang:
                total_incorrect += ligationFrequencies[overhang][other_overhang]
                total_incorrect += ligationFrequencies[overhang][_revcomp(other_overhang)]

                total_incorrect_revcomp += ligationFrequencies[_revcomp(overhang)][_revcomp(other_overhang)]
                total_incorrect_revcomp += ligationFrequencies[_revcomp(overhang)][other_overhang]

                if shifted_p1nt_overhang_positions[other_position] is not None:
                    total_incorrect_shifted_p1nt += ligationFrequencies[overhang][shifted_p1nt_overhang_positions[other_position]]
                    total_incorrect_shifted_p1nt += ligationFrequencies[overhang][_revcomp(shifted_p1nt_overhang_positions[other_position])]

                    total_incorrect_shifted_p1nt_revcomp += ligationFrequencies[_revcomp(overhang)][shifted_p1nt_overhang_positions[other_position]]
                    total_incorrect_shifted_p1nt_revcomp += ligationFrequencies[_revcomp(overhang)][_revcomp(shifted_p1nt_overhang_positions[other_position])]

                if shifted_n1nt_overhang_positions[other_position] is not None:
                    total_incorrect_shifted_n1nt += ligationFrequencies[overhang][shifted_n1nt_overhang_positions[other_position]]
                    total_incorrect_shifted_n1nt += ligationFrequencies[overhang][_revcomp(shifted_n1nt_overhang_positions[other_position])]

                    total_incorrect_shifted_n1nt_revcomp += ligationFrequencies[_revcomp(overhang)][shifted_n1nt_overhang_positions[other_position]]
                    total_incorrect_shifted_n1nt_revcomp += ligationFrequencies[_revcomp(overhang)][_revcomp(shifted_n1nt_overhang_positions[other_position])]

        probability = float(total_correct) / float(total_correct + total_incorrect + SHIFT_PROBABILITY * total_incorrect_shifted_p1nt + SHIFT_PROBABILITY * total_incorrect_shifted_n1nt )
        probability_revcomp = float(total_correct_revcomp) / float(total_correct_revcomp + total_incorrect_revcomp + SHIFT_PROBABILITY * total_incorrect_shifted_p1nt_revcomp + SHIFT_PROBABILITY * total_incorrect_shifted_n1nt_revcomp)

        ligation_fidelity = ligation_fidelity * probability
        ligation_fidelity_revcomp = ligation_fidelity_revcomp * probability_revcomp

    min_ligation_fidelity = min(ligation_fidelity, ligation_fidelity_revcomp)
    return min_ligation_fidelity

def prepare_numba_ligation_data(
    oligoPoolSpecification,
    enzyme,
    allowableOverhangPositions,
    shifted_p1nt_overhang_positions,
    shifted_n1nt_overhang_positions,
    revcomp_func,   # pass your _revcomp here
):
    """
    Build Numba-friendly arrays:
      - lig_matrix[i,j] = ligationFrequencies[overhang_i][overhang_j] as int64
      - rev_idx[i]      = index of revcomp(overhang_i)
      - allow_idx[pos]  = index of allowableOverhangPositions[pos]
      - shift_p[pos]    = index of shifted_p1nt_overhang_positions[pos] or -1
      - shift_n[pos]    = index of shifted_n1nt_overhang_positions[pos] or -1

    Returns:
      (lig_matrix, rev_idx, allow_idx, shift_p, shift_n, overhang_to_idx)
    """

    ligFreq = oligoPoolSpecification["ALL_LIGATION_COUNTS"][enzyme]

    # Collect all overhang keys from ligFreq (includes revcomps if present in your dataset)
    overhangs = sorted(ligFreq.keys())
    overhang_to_idx = {oh: i for i, oh in enumerate(overhangs)}
    N = len(overhangs)

    # Build dense int64 matrix for fast indexed access
    lig_matrix = np.zeros((N, N), dtype=np.int64)
    for oh, row in ligFreq.items():
        i = overhang_to_idx[oh]
        # row is dict-like: {other_overhang: count}
        for oh2, v in row.items():
            j = overhang_to_idx[oh2]
            lig_matrix[i, j] = int(v)

    # Reverse complement index mapping (must exist in ligFreq)
    rev_idx = np.empty(N, dtype=np.int32)
    for oh, i in overhang_to_idx.items():
        rc = revcomp_func(oh)
        try:
            rev_idx[i] = overhang_to_idx[rc]
        except KeyError as e:
            raise KeyError(
                f"revcomp overhang '{rc}' (from '{oh}') not found in ligation frequencies for enzyme={enzyme}"
            ) from e

    # Positions can be sparse; allocate arrays sized to max position
    max_pos = max(allowableOverhangPositions.keys()) if allowableOverhangPositions else 0

    allow_idx = np.full(max_pos + 1, -1, dtype=np.int32)
    for pos, oh in allowableOverhangPositions.items():
        allow_idx[pos] = overhang_to_idx[oh]

    shift_p = np.full(max_pos + 1, -1, dtype=np.int32)
    shift_n = np.full(max_pos + 1, -1, dtype=np.int32)

    # shifted_* dicts may be list/array-like or dict-like; handle both via indexing if possible
    # If they are lists/arrays indexed by position, max_pos must also be <= len-1.
    def get_shift(container, pos):
        try:
            return container[pos]
        except Exception:
            return container.get(pos, None)

    for pos in range(max_pos + 1):
        ohp = get_shift(shifted_p1nt_overhang_positions, pos)
        if ohp is not None:
            shift_p[pos] = overhang_to_idx[ohp]
        ohn = get_shift(shifted_n1nt_overhang_positions, pos)
        if ohn is not None:
            shift_n[pos] = overhang_to_idx[ohn]

    return lig_matrix, rev_idx, allow_idx, shift_p, shift_n, overhang_to_idx

@njit(cache=True)
def _ligation_fidelity_for_candidate(
    lig_matrix, rev_idx,
    base_overhang_idx, base_pos,
    cand_overhang_idx, cand_pos,
    shift_p, shift_n,
    SHIFT_PROBABILITY
):
    """
    Compute ligation fidelity for overhangSet = base + candidate,
    matching the original Python logic and addition order.
    """
    m = base_overhang_idx.size
    total_len = m + 1

    ligation_fidelity = 1.0
    ligation_fidelity_revcomp = 1.0

    for a in range(total_len):
        if a < m:
            over_i = base_overhang_idx[a]
            pos_i = base_pos[a]
        else:
            over_i = cand_overhang_idx
            pos_i = cand_pos

        rc_over_i = rev_idx[over_i]

        # total_correct, total_correct_revcomp
        total_correct = lig_matrix[over_i, over_i]
        total_correct_rc = lig_matrix[rc_over_i, rc_over_i]

        # total_incorrect starts with opposite (revcomp) as in original
        total_incorrect = lig_matrix[over_i, rc_over_i]
        total_incorrect_rc = lig_matrix[rc_over_i, over_i]

        # shifted initial incorrect counts (for this junction's own shifted)
        total_shift_p = 0
        total_shift_p_rc = 0
        sp_i = shift_p[pos_i]
        if sp_i != -1:
            total_shift_p = lig_matrix[over_i, rev_idx[sp_i]]
            total_shift_p_rc = lig_matrix[rc_over_i, sp_i]

        total_shift_n = 0
        total_shift_n_rc = 0
        sn_i = shift_n[pos_i]
        if sn_i != -1:
            total_shift_n = lig_matrix[over_i, rev_idx[sn_i]]
            total_shift_n_rc = lig_matrix[rc_over_i, sn_i]

        # inner loop over other_overhangs (same ordering as input set)
        for b in range(total_len):
            if b < m:
                other_j = base_overhang_idx[b]
                other_pos = base_pos[b]
            else:
                other_j = cand_overhang_idx
                other_pos = cand_pos

            if other_j != over_i:
                # exact same addition order as original Python code:
                total_incorrect += lig_matrix[over_i, other_j]
                total_incorrect += lig_matrix[over_i, rev_idx[other_j]]

                total_incorrect_rc += lig_matrix[rc_over_i, rev_idx[other_j]]
                total_incorrect_rc += lig_matrix[rc_over_i, other_j]

                sp_o = shift_p[other_pos]
                if sp_o != -1:
                    total_shift_p += lig_matrix[over_i, sp_o]
                    total_shift_p += lig_matrix[over_i, rev_idx[sp_o]]

                    total_shift_p_rc += lig_matrix[rc_over_i, sp_o]
                    total_shift_p_rc += lig_matrix[rc_over_i, rev_idx[sp_o]]

                sn_o = shift_n[other_pos]
                if sn_o != -1:
                    total_shift_n += lig_matrix[over_i, sn_o]
                    total_shift_n += lig_matrix[over_i, rev_idx[sn_o]]

                    total_shift_n_rc += lig_matrix[rc_over_i, sn_o]
                    total_shift_n_rc += lig_matrix[rc_over_i, rev_idx[sn_o]]

        # probability (exact same formula)
        denom = float(total_correct + total_incorrect) + SHIFT_PROBABILITY * float(total_shift_p + total_shift_n)
        denom_rc = float(total_correct_rc + total_incorrect_rc) + SHIFT_PROBABILITY * float(total_shift_p_rc + total_shift_n_rc)

        probability = float(total_correct) / denom if denom > 0.0 else 0.0
        probability_rc = float(total_correct_rc) / denom_rc if denom_rc > 0.0 else 0.0

        ligation_fidelity *= probability
        ligation_fidelity_revcomp *= probability_rc

    return ligation_fidelity if ligation_fidelity < ligation_fidelity_revcomp else ligation_fidelity_revcomp

@njit(cache=True, parallel=False)
def calculateWeights_numba(
    base_overhang_idx, base_pos,
    candidate_positions,
    allow_idx,
    lig_matrix, rev_idx,
    shift_p, shift_n,
    SHIFT_PROBABILITY=0.001
):
    """
    Compute weights for each candidate position, matching your original calculateWeights:
      candidate_overhang_set = selected_set + [(allowableOverhangPositions[candidate_position], candidate_position)]
      weight = calculateLigationFidelity(...)
    """
    out = np.empty(candidate_positions.size, dtype=np.float64)

    for t in prange(candidate_positions.size):
        cand_pos = candidate_positions[t]
        cand_over_i = allow_idx[cand_pos]  # overhang index at that position

        # If position not mapped, treat as 0 weight (or raise in Python wrapper)
        if cand_over_i < 0:
            out[t] = 0.0
            continue

        out[t] = _ligation_fidelity_for_candidate(
            lig_matrix, rev_idx,
            base_overhang_idx, base_pos,
            cand_over_i, cand_pos,
            shift_p, shift_n,
            SHIFT_PROBABILITY
        )

    return out

def multiLevelFragmentConstructor(oligoPoolSpecification, fragmentSequence, enzyme, landingPadOverhangs):
    fragmentWithCutSites = landingPadOverhangs[0] + oligoPoolSpecification['CUT_SITE_SEQUENCES'][enzyme][0] + fragmentSequence + oligoPoolSpecification['CUT_SITE_SEQUENCES'][enzyme][1] + landingPadOverhangs[1]
    return fragmentWithCutSites

def fragmentConstructor(oligoPoolSpecification, fragmentSequence, enzyme):
    fragmentWithCutSites = oligoPoolSpecification['CUT_SITE_SEQUENCES'][enzyme][0] + fragmentSequence + oligoPoolSpecification['CUT_SITE_SEQUENCES'][enzyme][1]
    return fragmentWithCutSites

class GoldenGatePool(object):

    OVERHANG_FILENAMES = ['BbsI_HF.csv', 'BsaI_HFv2.csv', 'BsmBI_v2.csv', 'Esp3I.csv', 'SapI.csv']
    ENZYME_LIGATION_FREQUENCIES = ['BsaI_HFv2_LigationFrequency.csv','BsmBI_v2_LigationFrequency.csv','Esp3I_LigationFrequency.csv','BbsI_HF_LigationFrequency.csv','SapI_LigationFrequency.csv']

    STD_ENZYME_NAMES = {'bsai' : 'BsaI',
                       'bbsi' : 'BbsI',
                       'bsmbi' : 'BsmBI',
                       'esp3i' : 'Esp3I',
                       'sapi' : 'SapI'
                       }
    NEB_ENZYME_NAMES = {'bsai' : 'BsaI_HFv2',
                       'bbsi' : 'BbsI_HF',
                       'bsmbi' : 'BsmBI_v2',
                       'esp3i' : 'Esp3I',
                       'sapi' : 'SapI'
                       }

    RECOGNITION_SEQUENCES = {'BsaI_HFv2' : ['GGTCTC', 'GAGACC'],
                            'BsaI'       : ['GGTCTC', 'GAGACC'],
                            'BbsI_HF'    :  ['GAAGAC', 'GTCTTC'],
                            'BbsI'       :  ['GAAGAC', 'GTCTTC'],
                            'BsmBI_v2'   :  ['CGTCTC', 'GAGACG'],
                            'BsmBI'      :  ['CGTCTC', 'GAGACG'],
                            'Esp3I'      :  ['CGTCTC', 'GAGACG'],
                            'SapI'       :  ['GCTCTTC','GAAGAGC']
                            }

    CUT_SITE_SEQUENCES =   {'BsaI_HFv2' :  ['GGTCTCA', 'AGAGACC'],
                            'BbsI_HF'   :  ['GAAGACAA', 'AAGTCTTC'],
                            'BsmBI_v2'  :  ['CGTCTCA', 'AGAGACG'],
                            'Esp3I'     :  ['CGTCTCA', 'AGAGACG'],
                            'SapI'      :  ['GCTCTTCA','AGAAGAGC']
                            }
    SPACER_LENGTHS =  { 'BsaI_HFv2' : 1, 'BbsI_HF' : 2, 'BsmBI_v2' : 1, 'Esp3I' : 1, 'SapI' : 1 }
    OVERHANG_LENGTH = { 'BsaI_HFv2' : 4, 'BbsI_HF' : 4, 'BsmBI_v2' : 4, 'Esp3I' : 4, 'SapI' : 3 }

    MINIMUM_FRAGMENT_LENGTH = 150
    MAXIMUM_OLIGO_LENGTH = 300

    LONG_ASSEMBLY_THRESHOLD = 3000 #bp
    LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH = 1500 #bp
    LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH = 500 #bp
    NUM_FRAGMENTS_PER_LONG_ASSEMBLY = 9 #11

    DEFAULT_PRIMER_SEQUENCE_LENGTH = 18 #nt
    DEFAULT_PRIMER_SEQUENCE_CONSTRAINT = 'N' * (DEFAULT_PRIMER_SEQUENCE_LENGTH-2) + 'WW'
    DEFAULT_TM_TARGET = 60.0 #degC
    TARGET_TM_TOLERANCE = 0.50 #degC
    DEFAULT_TOLERANCE_DELTA_TM = 0.20 #degC
    DEFAULT_TOLERANCE_DG_FOLDING = -1.0 #kcal/mol
    DEFAULT_TOLERANCE_DG_HOMODIMER = -5.0 #kcal/mol
    DEFAULT_TOLERANCE_DG_PRIMER_DIMER = -5.0 #kcal/mol
    DEFAULT_TM_SETTINGS = {'dnac1' : 500.0, 'Na' : 150.0, 'K' : 0.0, 'Mg' : 2.0, 'dNTPs' : 0.2, 'Tris' : 65.0, 'saltcorr' : 1} #Q5 DNA polymerase buffer
    DEFAULT_DNA_MODEL = PyVRNA(dangles=2, gquad=True, parameter_file='dna_mathews2004.par') if PyVRNA is not None else None
    DEFAULT_EXCLUDED_RESTRICTION_ENZYMES = ['BsaI','BbsI','SapI','BsmBI','Esp3I']
    DEFAULT_MAXIMUM_REPEAT_LENGTH = 10 #nt
    DEFAULT_PADDING_END = '5p' # choices are ['5p','3p']

    BsaIF = 'GGTCTCG'   #GGTCTCNXXXX
    BsaIR = 'GGAGACC'   #XXXXNGAGACC

    BbsIF = 'GAAGACAA'      #GAAGACNNXXXX
    BbsIR = 'TTGTCTTC'      #XXXXNNGTCTTC

    landingPadOverhangs =  ('TAAA', 'TTAG')

    def __init__(self, verbose = False, parallel = False, doScreening = True):

        self.verbose = verbose
        self.parallel = parallel
        self.doScreening = doScreening

        self.screening_threads = []

        # Initialize data structures to store mappings between combinatorial set, well_number, and assembly names.
        self.used_well_numbers = []
        self.used_combinatorial_sets = {}
        self.combo_set_to_wells = {}
        self.well_to_combo_set = {}
        self.well_to_assembly_names = {}
        self.getNextWell = self.getNextWellGenerator()

        self.loadOverhangs()

        if self.parallel:
            self.number_processes = int(os.environ.get("NUM_WORKERS", 8))

            # Set start method only if not already set
            if multiprocessing.get_start_method(allow_none=True) is None:
                multiprocessing.set_start_method('spawn')  # Safe to call only once. 'fork' caused Vienna RNA segfaults at high CPU load.

            self.pool = multiprocessing.Pool(self.number_processes)
            print("Running in Parallel mode with {} CPUs.".format(self.number_processes))
        else:
            print("Running in Serial mode.")
            self.pool = None

    def loadOverhangs(self):
        self.ALL_LIGATION_COUNTS = {}
        self.RANKED_OVERHANGS = {}
        self.ALL_OVERHANGS = {}

        for filename in self.ENZYME_LIGATION_FREQUENCIES:
            with open(os.path.join(current_directory,'data/' + filename), 'r') as fc:
                header = fc.readline().replace('\n','').replace('\r','')
                overhangs = [_revcomp(x) for x in header.split(',')[1:]]

                ligationFrequencies = {}
                for (n, line) in enumerate(fc.readlines()):
                    line = line.replace('\n','').replace('\r','')
                    words = line.split(',')

                    overhang_5p = words[0]
                    overhang_3p = overhangs[n]
                    ligation_counts = [int(word) for word in words[1:]]

                    assert overhang_5p == overhang_3p, "Mismatched overhangs {}, {} in filename {}".format(overhang_3p, overhang_5p, filename)

                    ligationFrequencies[overhang_5p] = {}
                    for (overhang_3p, counts) in zip(overhangs, ligation_counts):
                        ligationFrequencies[overhang_5p][overhang_3p] = counts

                enzyme = filename.replace('_LigationFrequency.csv','')
                self.ALL_LIGATION_COUNTS[enzyme] = ligationFrequencies

                overhangAverageFidelities = {}
                overhangs = list(self.ALL_LIGATION_COUNTS[enzyme].keys())

                # Remove palidromic overhangs. Palidromic overhangs allow two of the same overhang to ligate when DNA orientation is flipped.
                if   self.OVERHANG_LENGTH[enzyme] == 2:
                    palidromic_overhangs = ["".join(x) + _revcomp("".join(x)) for x in itertools.product(['A','T','G','C'],repeat = 1)]
                elif self.OVERHANG_LENGTH[enzyme] == 4:
                    palidromic_overhangs = ["".join(x) + _revcomp("".join(x)) for x in itertools.product(['A','T','G','C'],repeat = 2)] #['AATT', 'ACGT', 'AGCT', 'ATAT', 'CATG', 'CCGG', 'CGCG', 'CTAG', 'GATC', 'GCGC', 'GGCC', 'GTAC', 'TATA', 'TCGA', 'TGCA', 'TTAA']
                elif self.OVERHANG_LENGTH[enzyme] == 6:
                    palidromic_overhangs = ["".join(x) + _revcomp("".join(x)) for x in itertools.product(['A','T','G','C'],repeat = 3)]
                else:
                    palidromic_overhangs = []

                for pal_overhang in palidromic_overhangs:
                    overhangs.remove(pal_overhang)

                for overhang_5p in overhangs:
                    total_correct = self.ALL_LIGATION_COUNTS[enzyme][overhang_5p][overhang_5p]
                    total_incorrect = 0
                    for overhang_3p in overhangs:
                        if overhang_3p != overhang_5p: total_incorrect +=  self.ALL_LIGATION_COUNTS[enzyme][overhang_5p][overhang_3p]

                    overhangAverageFidelities[overhang_5p] = float(total_correct) / float(total_correct + total_incorrect)

                self.ALL_OVERHANGS[enzyme] = overhangs
                self.RANKED_OVERHANGS[enzyme] = overhangAverageFidelities

    def getNextWellGenerator(self):
        counter = 0
        while True:
            self.used_well_numbers.append(counter)
            yield counter
            counter += 1

    def updateCombinatorialSet(self, input_combinatorial_set):

        if input_combinatorial_set is None:
            next_combinatorial_set = max([-1] + list(self.used_combinatorial_sets.values())) + 1
            self.used_combinatorial_sets[None] = next_combinatorial_set
            return next_combinatorial_set
        elif input_combinatorial_set in self.used_combinatorial_sets:
            return self.used_combinatorial_sets[input_combinatorial_set]
        else:
            next_combinatorial_set = max([-1] + list(self.used_combinatorial_sets.values())) + 1
            self.used_combinatorial_sets[input_combinatorial_set] = next_combinatorial_set
            return next_combinatorial_set

    def findSequenceDifferences(self, inputSpec, sequenceSpecList):

        if inputSpec['no_overhangs_position_range_list'] is not None:
            no_overhangs_position_range_list = inputSpec['no_overhangs_position_range_list'].copy()
        else:
            no_overhangs_position_range_list = []

        combinatorial_set = inputSpec['combinatorial_set']
        sequenceList = [str(spec['nucleotide_sequence']) for spec in sequenceSpecList if spec['combinatorial_set'] == combinatorial_set]

        # if the assembly is a one-pot library, but the # of sequences in the library is less than the cutoff
        ALIGNMENT_CUTOFF = 1000
        if len(sequenceList) > 1 and len(sequenceList) < ALIGNMENT_CUTOFF:
            from GeneticSystemsCalculator_Build_Helpers import find_nonmatching_ranges_from_unaligned, aligned_ranges_to_sequence_ranges

            # Carry out a multiple sequence alignment and identify (begin, end) ranges where mismatches or gaps occur
            aligned_ranges, aln_unique, inverse = find_nonmatching_ranges_from_unaligned(sequenceList)

            # Convert aligned-coordinate ranges into ungapped coordinates for the inputted sequence
            target_sequence = str(inputSpec["nucleotide_sequence"]).upper()
            target_index = next(
                i for i, seq in enumerate(sequenceList)
                if str(seq).upper() == target_sequence
            )
            seq_nonmatching_ranges = aligned_ranges_to_sequence_ranges(
                aligned_ranges,
                aln_unique,
                inverse,
                target_index=target_index,
            )
            return no_overhangs_position_range_list + seq_nonmatching_ranges
        else:
            return no_overhangs_position_range_list

    def assignWells(self, assembly):

        assembly_name = assembly['name']
        assembly_id = assembly['id']
        combinatorial_set = assembly['combinatorial_set']

        if 'long_assembly_fragment_number' in assembly:
            long_assembly_fragment_number = assembly['long_assembly_fragment_number']
        else:
            long_assembly_fragment_number = None

        if long_assembly_fragment_number is None: # this is a level one assembly
            if combinatorial_set in self.combo_set_to_wells:
                well_number = self.combo_set_to_wells[combinatorial_set][0]
                self.well_to_assembly_names[well_number].append(assembly_name)
            else:
                well_number = next(self.getNextWell)
                self.combo_set_to_wells[combinatorial_set] = [well_number]
                self.well_to_combo_set[well_number] = combinatorial_set
                self.well_to_assembly_names[well_number] = [assembly_name]
        else:
            # this assumes long assembly fragments are sorted in numerical order [long_assembly_fragment_number of 0 is seen first]
            # this is true by construction
            if combinatorial_set in self.combo_set_to_wells:
                starting_well_number = self.combo_set_to_wells[combinatorial_set][0]
                well_number = starting_well_number + long_assembly_fragment_number
                if not (well_number in self.used_well_numbers):
                    well_number = next(self.getNextWell)

                self.combo_set_to_wells[combinatorial_set].append(well_number)
                if well_number in self.well_to_combo_set:
                    self.well_to_combo_set[well_number] = combinatorial_set
                    self.well_to_assembly_names[well_number].append(assembly_name)
                else:
                    self.well_to_combo_set[well_number] = combinatorial_set
                    self.well_to_assembly_names[well_number] = [assembly_name]

            else:
                well_number = next(self.getNextWell)
                self.combo_set_to_wells[combinatorial_set] = [well_number]
                self.well_to_combo_set[well_number] = combinatorial_set
                self.well_to_assembly_names[well_number] = [assembly_name]

        assembly['well_number'] = well_number
        return assembly

    def addMultipleAssemblies(self, oligoPoolSpecification, sequenceSpecList):

        if 'library_specifications' in oligoPoolSpecification and len(oligoPoolSpecification['library_specifications']) > 0:
            lastID = oligoPoolSpecification['library_specifications'][-1]['id']
        else:
            lastID = 0

        oligoPoolSpecification['CUT_SITE_SEQUENCES'] = self.CUT_SITE_SEQUENCES
        oligoPoolSpecification['OVERHANG_LENGTH'] = self.OVERHANG_LENGTH
        oligoPoolSpecification['ALL_LIGATION_COUNTS'] = self.ALL_LIGATION_COUNTS
        oligoPoolSpecification['RANKED_OVERHANGS'] = self.RANKED_OVERHANGS

        oligoPoolSpecification['MINIMUM_FRAGMENT_LENGTH'] = oligoPoolSpecification.get('MINIMUM_FRAGMENT_LENGTH', self.MINIMUM_FRAGMENT_LENGTH)
        oligoPoolSpecification['DEFAULT_PRIMER_SEQUENCE_LENGTH'] = oligoPoolSpecification.get('DEFAULT_PRIMER_SEQUENCE_LENGTH', self.DEFAULT_PRIMER_SEQUENCE_LENGTH)
        oligoPoolSpecification['LONG_ASSEMBLY_THRESHOLD'] = oligoPoolSpecification.get('LONG_ASSEMBLY_THRESHOLD', self.LONG_ASSEMBLY_THRESHOLD)
        oligoPoolSpecification['LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH'] = oligoPoolSpecification.get('LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH', self.LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH)
        oligoPoolSpecification['LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH'] = oligoPoolSpecification.get('LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH', self.LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH)
        oligoPoolSpecification['NUM_FRAGMENTS_PER_LONG_ASSEMBLY'] = oligoPoolSpecification.get('NUM_FRAGMENTS_PER_LONG_ASSEMBLY', self.NUM_FRAGMENTS_PER_LONG_ASSEMBLY)

        maximum_oligo_length = oligoPoolSpecification.get('maximum_oligo_length', self.MAXIMUM_OLIGO_LENGTH)
        primer_sequence_length = oligoPoolSpecification.get('primer_sequence_length', self.DEFAULT_PRIMER_SEQUENCE_LENGTH)

        for spec in sequenceSpecList:
            spec['sequence'] = spec['nucleotide_sequence'].upper()
            enzyme = spec.get('enzyme', spec.get('assembly_enzyme', 'BsaI_HFv2'))
            spec['no_overhangs_position_range_list'] = spec.get('no_overhangs_position_range_list', None)

            overhang_length = oligoPoolSpecification['OVERHANG_LENGTH'][enzyme]
            flanking_sequence_length = len(oligoPoolSpecification['CUT_SITE_SEQUENCES'][enzyme][0]) + len(oligoPoolSpecification['CUT_SITE_SEQUENCES'][enzyme][1]) + 2 * primer_sequence_length + overhang_length
            maximum_fragment_length = maximum_oligo_length - flanking_sequence_length
            spec['maximum_fragment_length'] = maximum_fragment_length

            if spec['is_circular']:
                circular_overhang_anchor = spec['sequence'][0:overhang_length]
                if circular_overhang_anchor in self.RANKED_OVERHANGS[enzyme] and self.RANKED_OVERHANGS[enzyme][circular_overhang_anchor] > 0.85:
                    spec['circular_anchor'] = 0
                else:
                    spec['circular_anchor'] = None
            else:
                spec['circular_anchor'] = None

            # combinatorial_set defines which assemblies take place together inside the same well
            # if two assemblies have the same value of combinatorial_set, their assemblies take place in the same well (combinatorial assembly)
            if 'combinatorial_set' in spec:
                spec['combinatorial_set'] = self.updateCombinatorialSet(spec['combinatorial_set'])
            else:
                spec['combinatorial_set'] = self.updateCombinatorialSet(None)

            if oligoPoolSpecification['level_two_assembly_enzyme'] is not None and len(spec['sequence']) >= oligoPoolSpecification['LONG_ASSEMBLY_THRESHOLD']:
                spec['long_assembly'] = True
                spec['level_two_assembly_enzyme'] = spec.get('level_two_assembly_enzyme', oligoPoolSpecification['level_two_assembly_enzyme'])
                # assembly level = ceiling( log(sequenceLength/fragmentLength) / log(numFragmentsPerAssembly)) + 1
                spec['assembly_levels'] = math.ceil(math.log(len(spec['sequence']) / oligoPoolSpecification['LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH'])
                                               / math.log(oligoPoolSpecification['NUM_FRAGMENTS_PER_LONG_ASSEMBLY'])) + 1
            else:
                spec['long_assembly'] = False
                spec['assembly_levels'] = 1

            # When multiple assemblies take place in the same well, they share the same combinatorial_set
            # We can calculate no_overhangs_position_range_list based on nucleotide sequence differences across assemblies with the same combinatorial_set
            # This will prevent the selection of overhangs in one-pot assemblies where the overhang positions, sequences might change across assemblies
            spec['no_overhangs_position_range_list'] = self.findSequenceDifferences(spec, sequenceSpecList)

        if self.doScreening: self.initialize_screen(sequenceSpecList)

        if self.verbose: print(f'sequenceSpecList: {sequenceSpecList}')
        # Regular Assemblies First (assembly_levels = 1)
        inputList = [(oligoPoolSpecification, spec['sequence'], oligoPoolSpecification['MINIMUM_FRAGMENT_LENGTH'], spec['maximum_fragment_length'],
                      spec['assembly_enzyme'], spec['is_circular'], lastID + n, spec['name'], spec['combinatorial_set'], spec['no_overhangs_position_range_list'],
                      spec['circular_anchor'],
                      self.verbose) for (n, spec) in enumerate(sequenceSpecList) if spec['long_assembly'] == False]

        if self.parallel:
            assemblyList = self.pool.map(createAssembly, inputList)
            self.pool.close()
            self.pool.join()
        else:
            assemblyList = list(map(createAssembly, inputList))

        successfulLevelOneAssemblies = [assembly for (Success, assembly, error_message) in assemblyList if Success]
        errorMessages = [(assembly, error_message) for (n, (Success, assembly, error_message)) in enumerate(assemblyList) if not Success]

        # Assign well_number for successful level one assemblies
        finalLevelOneAssemblyList = [self.assignWells(assembly) for assembly in successfulLevelOneAssemblies]

        # Long Assemblies Second (assembly_levels >= 2)
        inputList = [(oligoPoolSpecification, spec['sequence'], spec['assembly_enzyme'], spec['level_two_assembly_enzyme'], spec['is_circular'], lastID + n, spec['name'], spec['combinatorial_set'], spec['no_overhangs_position_range_list'], spec['circular_anchor'],
                      self.verbose) for (n, spec) in enumerate(sequenceSpecList) if spec['long_assembly'] == True]

        finalLevelTwoAssemblyList = []
        for input_specifications in inputList:
            (Success, assemblyList, error_message_list) = createLongAssembly(input_specifications, self.pool)
            if Success:
                for assembly in assemblyList:
                    #Assign well to assembly
                    assembly = self.assignWells(assembly)
                    finalLevelTwoAssemblyList.append(assembly)
            else:
                for (errorAssembly, error_message) in zip(assemblyList, error_message_list):
                    errorMessages.append((errorAssembly, error_message))

        finalSuccessfulAssemblies = finalLevelOneAssemblyList + finalLevelTwoAssemblyList

        if 'library_specifications' in oligoPoolSpecification:
            oligoPoolSpecification['library_specifications'] += finalSuccessfulAssemblies
        else:
            oligoPoolSpecification['library_specifications'] = finalSuccessfulAssemblies

        if 'failed_assemblies' in oligoPoolSpecification:
            oligoPoolSpecification['failed_assemblies'] += errorMessages
        else:
            oligoPoolSpecification['failed_assemblies'] = errorMessages

        return oligoPoolSpecification

    def selectOligoLength(self, oligoPoolSpecification):

        primer_sequence_length = oligoPoolSpecification.get('primer_sequence_length', self.DEFAULT_PRIMER_SEQUENCE_LENGTH)
        maxFragmentLength = 0
        minFragmentLength = 10000000
        assemblyList = oligoPoolSpecification['library_specifications']
        for assembly in assemblyList:
            fragmentList = assembly['fragment_list']
            for fragment in fragmentList:
                maxFragmentLength = max(maxFragmentLength, len(fragment['fragment_with_cut_sites']))
                minFragmentLength = min(minFragmentLength, len(fragment['fragment_with_cut_sites']))

        maximumOligoLength = maxFragmentLength + 2 * primer_sequence_length
        minimumOligoLength = minFragmentLength + 2 * primer_sequence_length
        oligoPoolSpecification['maximumOligoLength'] = maximumOligoLength
        oligoPoolSpecification['minimumOligoLength'] = minimumOligoLength
        return oligoPoolSpecification

    def designPadding(self, oligoPoolSpecification):

        primer_sequence_length = oligoPoolSpecification.get('primer_sequence_length', self.DEFAULT_PRIMER_SEQUENCE_LENGTH)
        padding_end = oligoPoolSpecification.get('padding_end', self.DEFAULT_PADDING_END)
        maximumOligoLength = oligoPoolSpecification.get('maximum_oligo_length', self.MAXIMUM_OLIGO_LENGTH)
        excluded_restriction_enzymes = ['BsaI','BbsI','Esp3I','SapI']
        excludeSequenceList = [restriction_enzyme_sites[enzyme] for enzyme in excluded_restriction_enzymes] + ['AAAA','TTTT','GGGG','CCCC']
        excludeSequenceLength = max([len(x) for x in excludeSequenceList])

        # reuse same padding sequence for the same-sequence oligos (enables removal of same-sequence oligos from oligopool)
        padding_sequences_dict = {}

        def check_for_excluded_sequences(seq, edge_sequence = '', padding_end = None):
            # Check padding sequences for excluded sites
            # Forward Strand
            for regexp in excludeSequenceList:
                if padding_end == '5p':
                    extended_sequence = seq + edge_sequence[0:(len(regexp)-1)]
                elif padding_end == '3p':
                    extended_sequence = edge_sequence[-(len(regexp)-1):] + seq
                else:
                    extended_sequence = seq

                if re.search(regexp, extended_sequence):
                    return False

            # Reverse Strand
            for regexp in excludeSequenceList:
                if padding_end == '5p':
                    extended_sequence = seq + edge_sequence[0:(len(regexp)-1)]
                elif padding_end == '3p':
                    extended_sequence = edge_sequence[-(len(regexp)-1):] + seq
                else:
                    extended_sequence = seq

                if re.search(regexp, _revcomp(extended_sequence)):
                    return False

            return True

        librarySpecs = oligoPoolSpecification['library_specifications']
        numAssemblies = len(librarySpecs)
        for assembly in librarySpecs:
            libraryID = assembly['id']
            fragmentList = assembly['fragment_list']
            for fragment in fragmentList:
                paddingLength = maximumOligoLength - (len(fragment['fragment_with_cut_sites']) + 2 * primer_sequence_length)
                assert paddingLength >= 0, "A fragment sequence from assembly [id: {}] is too long [{} nt]! No room for primer binding sites [2 x {} nt]".format(libraryID, len(fragment['fragment_with_cut_sites']), primer_sequence_length)

                if fragment['fragment_with_cut_sites'] in padding_sequences_dict:
                    paddingSequence = padding_sequences_dict[fragment['fragment_with_cut_sites']]
                    fragment['padding_sequence'] = paddingSequence

                    if padding_end == '3p':
                        fragment['fragment_with_cut_sites_and_padding'] = fragment['fragment_with_cut_sites'] + paddingSequence

                    elif padding_end == '5p':
                        fragment['fragment_with_cut_sites_and_padding'] = paddingSequence + fragment['fragment_with_cut_sites']
                    else: #default is 5p
                        fragment['fragment_with_cut_sites_and_padding'] = paddingSequence + fragment['fragment_with_cut_sites']

                else:
                    CHUNK_SIZE = max(10, excludeSequenceLength+1)
                    MAX_EDGE_FIX_TRIES = 100  # safety to avoid infinite loops
                    if padding_end == '3p':
                        edge_chunk_index = 0
                    elif padding_end == '5p':
                        edge_chunk_index = -1
                    else: #5p
                        edge_chunk_index = -1

                    while True:
                        chunks = []
                        paddingSequence = ""

                        # 1. Build padding sequence in chunks
                        while len(paddingSequence) < paddingLength:
                            remaining = paddingLength - len(paddingSequence)
                            this_chunk_size = min(CHUNK_SIZE, remaining)

                            # test chunks until one passes the excluded-sequence test
                            while True:
                                chunk_seq = "".join(random.choice(['A', 'T', 'G', 'C']) for _ in range(this_chunk_size))

                                # test just this chunk
                                if check_for_excluded_sequences(paddingSequence + chunk_seq):
                                    chunks.append(chunk_seq)
                                    paddingSequence = "".join(chunks)
                                    #print(f'Designed correct chunk: {chunk_seq}')
                                    break
                                else:
                                    # chunk failed; try a different random chunk
                                    #print(f'Designed incorrect chunk, repeating: {chunk_seq}')
                                    continue

                        # 2. Check extended sequence; if it passes, we're done
                        if check_for_excluded_sequences(paddingSequence, fragment['fragment_with_cut_sites'], padding_end):
                            Passed = True
                            break

                        # 3. Extended sequence failed: try to fix the edge by modifying the edge chunk
                        # Does fragment_with_cut_sites contain an excluded sequence? If so, remove it.
                        while True:
                            new_chunk = "".join(random.choice(['A', 'T', 'G', 'C']) for _ in range(len(chunks[edge_chunk_index])))
                            chunks[edge_chunk_index] = new_chunk
                            paddingSequence = "".join(chunks)

                            if check_for_excluded_sequences(paddingSequence, edge_sequence = fragment['fragment_with_cut_sites'], padding_end = padding_end):
                                break  # accept this new chunk

                    fragment['padding_sequence'] = paddingSequence

                    if padding_end == '3p':
                        fragment['fragment_with_cut_sites_and_padding'] = fragment['fragment_with_cut_sites'] + paddingSequence
                    elif padding_end == '5p':
                        fragment['fragment_with_cut_sites_and_padding'] = paddingSequence + fragment['fragment_with_cut_sites']
                    else: #default is 5p
                        fragment['fragment_with_cut_sites_and_padding'] = paddingSequence + fragment['fragment_with_cut_sites']

                    # Save padding sequence to dictionary using fragment sequence as key
                    padding_sequences_dict[fragment['fragment_with_cut_sites']] = paddingSequence

        return oligoPoolSpecification

    def designPrimers(self, oligoPoolSpecification, verbose = True):
        excluded_restriction_enzymes = ['BsaI','BbsI','Esp3I','SapI'] #oligoPoolSpecification.get('excluded_restriction_enzymes', self.DEFAULT_EXCLUDED_RESTRICTION_ENZYMES)
        Tm_target = oligoPoolSpecification.get('PCR_Tm_target', self.DEFAULT_TM_TARGET)
        Tm_target_range = [ Tm_target - self.TARGET_TM_TOLERANCE, Tm_target + self.TARGET_TM_TOLERANCE ]
        tolerance_delta_Tm = oligoPoolSpecification.get('tolerance_delta_Tm', self.DEFAULT_TOLERANCE_DELTA_TM)
        tolerance_dG_folding = oligoPoolSpecification.get('tolerance_dG_folding', self.DEFAULT_TOLERANCE_DG_FOLDING)
        tolerance_dG_homodimer = oligoPoolSpecification.get('tolerance_dG_homodimer', self.DEFAULT_TOLERANCE_DG_HOMODIMER)
        tolerance_primer_dimer_folding_energy = oligoPoolSpecification.get('tolerance_primer_dimer_folding_energy', self.DEFAULT_TOLERANCE_DG_PRIMER_DIMER)
        Tm_settings = oligoPoolSpecification.get('PCR_Tm_settings', self.DEFAULT_TM_SETTINGS)
        DNA_folding_model = oligoPoolSpecification.get('DNA_folding_model', self.DEFAULT_DNA_MODEL)
        if DNA_folding_model is None:
            raise ImportError("Could not import PyVRNA. Ensure ViennaRNA's Python module is available before running primer design.")
        primer_sequence_length = oligoPoolSpecification.get('primer_sequence_length', self.DEFAULT_PRIMER_SEQUENCE_LENGTH)
        primer_sequence_constraint = oligoPoolSpecification.get('primer_sequence_constraint', self.DEFAULT_PRIMER_SEQUENCE_CONSTRAINT)
        structureConstraint = "x" * len(primer_sequence_constraint)
        Lmax = oligoPoolSpecification.get('primer_maximum_repeat_length', self.DEFAULT_MAXIMUM_REPEAT_LENGTH)
        padding_end = oligoPoolSpecification.get('padding_end', self.DEFAULT_PADDING_END)

        librarySpecs = oligoPoolSpecification['library_specifications']

        # Multiple assemblies could take place inside one well
        # Count the number of wells needed across all assemblies
        # The number of wells will be the number of unique primer binding sites and primers
        uniqueWells = list(set([spec['well_number'] for spec in librarySpecs]))
        numAssemblies = len(uniqueWells)
        numCandidatesPerRound = max(1000, numAssemblies * 10)

        # a library is a set of assemblies
        # each assembly uses one pair of PCR primers and one TypeIIS enzyme

        # Generate list of DNA assembly enzymes used to create a list of excluded sequences
        excludeSequenceList = [restriction_enzyme_sites[enzyme] for enzyme in excluded_restriction_enzymes] + ['AAAA','TTTT','GGGG','CCCC']
        excludeCutSiteList = [restriction_enzyme_sites[enzyme] for enzyme in excluded_restriction_enzymes]
        excludeSequenceLength = max([len(x) for x in excludeSequenceList])

        ## Need to create background sequence covering all fragment_with_cut_sites across all assemblies
        backgroundSequenceList = oligoPoolSpecification.get('background_sequence_list', [])
        for assembly in librarySpecs:
            fragmentList = assembly['fragment_list']
            for fragment in fragmentList:
                backgroundSequenceList.append(str(fragment['fragment_with_cut_sites_and_padding']).upper())

        background_directory = tempfile.mkdtemp()
        background = nrpcalc.background(path=background_directory, Lmax = Lmax)
        background.multiadd(backgroundSequenceList)

        from GeneticSystemsCalculator_ModelFunctions import createGlobalModelFunction

        globalRules = {'rule_Tm_bounds' : True}
        ruleInputs = { 'MeltingTemperature' : {'TM_BOUNDS' : Tm_target_range, 'TM_SETTINGS' : Tm_settings},
                     }
        verbosity = 1
        globalModelFunctionInputs = {  'sequence_length' : primer_sequence_length,
                                       'upstream_sequence' : '',
                                       'downstream_sequence' : '',
                                       'activeRules' : globalRules,
                                       'ruleInputs' : ruleInputs,
                                       'verbosity' : verbosity,
                                       }
        globalModelFunction = createGlobalModelFunction(globalModelFunctionInputs)

        selectedPrimers = []
        numFails = 0
        while True:
            candidatePrimers = []
            results = nrpcalc.maker(seq_constr = primer_sequence_constraint,
                                struct_constr = structureConstraint,
                                part_type ='DNA',
                                Lmax = Lmax,
                                target_size = numCandidatesPerRound,
                                internal_repeats = False,
                                background = background,
                                struct_type='mfe',
                                seed=None,
                                synth_opt=True,
                                local_model_fn=None,
                                global_model_fn=globalModelFunction,
                                jump_count=10, fail_count=1000, output_file=None, verbose=False)

            if len(results) < 2:
                print('** PRIMER DESIGN ** Could not generate non-repetitive primers. Stopping.')
                return None

            candidateSequences = [seq for (num, seq) in list(results.items())]

            # Create Function Check for any Excluded Sequences
            def hasExcludedSequencesFcn(seq, assembly = None):
                # Check sequence for excluded sites
                for regexp in excludeSequenceList:
                    if re.search(regexp, seq): return True
                    if re.search(regexp, _revcomp(seq)): return True

                if assembly is not None:
                    fragmentSequences = [fragment['fragment_with_cut_sites_and_padding'] for fragment in assembly['fragment_list']]

                    for regexp in excludeCutSiteList:
                        for context in fragmentSequences:
                            context_length = len(regexp)-1
                            extended_sequence = seq + context[0:context_length]
                            if re.search(regexp, extended_sequence):
                                #print(f'{regexp} is inside {seq} + {context[0:context_length]} ?')
                                return True

                            extended_sequence = context[-context_length:] + seq
                            if re.search(regexp, extended_sequence):
                                #print(f'{regexp} is inside {context[-context_length:]} + {seq}?')
                                return True
                return False

            for seq in candidateSequences:
                if not hasExcludedSequencesFcn(seq):
                    fold = DNA_folding_model.RNAfold(seq)
                    if fold.energy > tolerance_dG_folding:
                        homodimer_fold = DNA_folding_model.RNAcofold(sequences = [ seq , seq])
                        if homodimer_fold.energy > tolerance_dG_homodimer:
                            Tm = mt.Tm_NN(seq, **Tm_settings)
                            incompatibleAssemblies = []
                            for (n, assembly) in enumerate(librarySpecs):
                                if hasExcludedSequencesFcn(seq, assembly):
                                    incompatibleAssemblies.append(n)
                                    #print('seq: ', seq)
                            item = {'sequence' : seq, 'Tm' : Tm, 'dG_folding' : fold.energy, 'structure' : fold.structure,
                                    'dG_homodimer' : homodimer_fold.energy, 'homodimer_structure' : homodimer_fold.structure,
                                    'incompatibleAssemblies' : copy.copy(incompatibleAssemblies)
                                    }
                            candidatePrimers.append( item )

            #Identify Sets of Primer Binding Sites with Similar Melting Temperatures
            acceptablePrimers = {}
            compatiblePrimers = {}

            numSteps = int( (Tm_target_range[1] - Tm_target_range[0]) / tolerance_delta_Tm) + 1
            for Tm_target in np.linspace(Tm_target_range[0], Tm_target_range[1], numSteps):
                acceptablePrimers[Tm_target] = []
                compatiblePrimers[Tm_target] = []

                for candidate in candidatePrimers:
                    if abs(candidate['Tm'] - Tm_target) < tolerance_delta_Tm:
                        acceptablePrimers[Tm_target].append( candidate )

                if len(acceptablePrimers[Tm_target]) >= 2:

                    primerPairs = itertools.combinations(acceptablePrimers[Tm_target],2)
                    for primerPair in primerPairs:
                        # Calculate primer dimer interactions for all pairs of forward & reverse primers
                        fold = DNA_folding_model.RNAcofold(sequences = [ primerPair[0]['sequence'] , primerPair[1]['sequence']])
                        structure_5p, structure_3p = fold.structure.split('&')[0], fold.structure.split('&')[1]
                        no_terminal_binding_5p = all([x == '.' for x in structure_5p[-3:]])
                        no_terminal_binding_3p = all([x == '.' for x in structure_3p[0:3]])
                        if fold.energy > self.DEFAULT_TOLERANCE_DG_PRIMER_DIMER and no_terminal_binding_5p and no_terminal_binding_3p:
                            primerPair[0]['dG_primer_dimer'] = fold.energy
                            primerPair[1]['dG_primer_dimer'] = fold.energy
                            primerPair[0]['structure_primer_dimer'] = fold.structure
                            primerPair[1]['structure_primer_dimer'] = fold.structure

                            compatiblePrimers[Tm_target].append( primerPair )


            #Select Compatible Primers with the Highest Melting Temperatures
            usedPrimerSequences = []

            sortedByTm = sorted([Tm for (Tm, primerPairs) in list(compatiblePrimers.items()) if len(primerPairs) > 0], reverse=True)
            for Tm_index in range(len(sortedByTm)):
                selected_Tm = sorted([Tm for (Tm, primerPairs) in list(compatiblePrimers.items()) if len(primerPairs) > 0], reverse=True)[Tm_index]
                primerPairs = compatiblePrimers[selected_Tm]
                for (primer1, primer2) in primerPairs:
                    if primer1['sequence'] in usedPrimerSequences or primer2['sequence'] in usedPrimerSequences:
                        pass # primer already used in a primer pair
                    else:
                        selectedPrimers.append( (primer1, primer2) )
                        usedPrimerSequences.append( primer1['sequence'] )
                        usedPrimerSequences.append( primer2['sequence'] )
                        background.add(primer1['sequence'])
                        background.add(primer2['sequence'])

            if verbose:
                maxAcceptablePrimersAtSameTM = max([len(x) for x in list(acceptablePrimers.values())])
                maxCompatiblePrimersAtSameTM = max([len(x) for x in list(compatiblePrimers.values())])
                print('# Non-repetitive, target Tm Primers: {}. # Candidate Primers: {}. # Max Acceptable Primers at same Tm: {}. # Max Compatible Primers at same Tm: {}.'.format(len(candidateSequences), len(candidatePrimers), maxAcceptablePrimersAtSameTM, maxCompatiblePrimersAtSameTM))

            print('** PRIMER DESIGN ** Designed {} / {} primers **'.format(len(selectedPrimers), numAssemblies))

            # Pair primers with assemblies, but exclude any primers with incompatibleAssemblies
            if len(selectedPrimers) > numAssemblies:
                needMorePrimers = False
                usedPrimers = []

                for well in uniqueWells:

                    primerID = 0
                    while primerID < len(selectedPrimers):

                        # This will make primerID the next in the unused primer series
                        if primerID in usedPrimers:
                            primerID += 1
                            continue

                        primerPair = selectedPrimers[primerID]

                        primerOK = True
                        for (n, assembly) in enumerate(librarySpecs):
                            if assembly['well_number'] == well:
                                if (n in primerPair[0]['incompatibleAssemblies']) or (n in primerPair[1]['incompatibleAssemblies']):
                                    print('Primer pair #{} is incompatible for assembly #{}. Seq 5p: {}. Seq 3p: {}'.format(primerID, n, primerPair[0]['sequence'],primerPair[1]['sequence']))
                                    primerOK = False
                                    break

                        if primerOK:
                            for (n, assembly) in enumerate(librarySpecs):
                                if assembly['well_number'] == well:
                                    assembly['primer_5p'] = primerPair[0]
                                    assembly['primer_3p'] = primerPair[1]
                                    fragmentList = assembly['fragment_list']
                                    for fragment in fragmentList:
                                        fragment['oligo_sequence'] = primerPair[0]['sequence'] + fragment['fragment_with_cut_sites_and_padding'] + _revcomp(primerPair[1]['sequence'])

                            usedPrimers.append(primerID)

                            break # breaks while loop

                        else:
                            primerID += 1
                            continue

                    if primerID >= len(selectedPrimers):
                        needMorePrimers = True
                        print(f'Used up all primers: {usedPrimers}\n\n\n\n\n\n\n')
                        break

            else:
                needMorePrimers = True

            if not needMorePrimers: break

        #Remove background
        background.drop()

        try:
            shutil.rmtree(background_directory)
        except:
            pass

        return oligoPoolSpecification

    def qualityControl(self, oligoPoolSpecification):

        if not 'failed_assemblies' in oligoPoolSpecification: oligoPoolSpecification['failed_assemblies'] = []

        librarySpecs = oligoPoolSpecification['library_specifications']
        uniqueWells = list(set([spec['well_number'] for spec in oligoPoolSpecification['library_specifications']]))

        for assembly in librarySpecs:
            id = assembly['id']
            name = assembly['name']
            enzyme = assembly['assembly_enzyme']
            assemblyEnzymeSites = self.RECOGNITION_SEQUENCES[enzyme]

            for fragment in assembly['fragment_list']:
                oligo_sequence = fragment['oligo_sequence']
                for site in assemblyEnzymeSites:
                    siteCounts = oligo_sequence.count(site)
                    if siteCounts == 1: pass
                    else:
                        error_message = f'OLIGO in ASSEMBLY {name} [{id}] contains {siteCounts} {enzyme} sites {site}, when only 1 is allowed.'
                        print(error_message)
                        oligoPoolSpecification['failed_assemblies'].append( (assembly, error_message) )

        for (counter, well) in enumerate(uniqueWells):
            assemblyList = [assembly for assembly in librarySpecs if assembly['well_number'] == well]
            numFragmentsPerOnePotAssembly = [len(assembly['fragment_list']) for assembly in assemblyList ]

            pairedOverhangs_5p = {}
            pairedOverhangs_3p = {}

            for assembly in assemblyList:
                id = assembly['id']
                name = assembly['name']
                fragmentList = assembly['fragment_list']

                for (oligoNum, fragment) in enumerate(fragmentList):
                    overhang_5p = fragment['overhang_5p']
                    overhang_3p = fragment['overhang_3p']

                    if overhang_5p in pairedOverhangs_5p:
                        overhang_3p = pairedOverhangs_5p[overhang_5p]
                        expected_overhang_5p = pairedOverhangs_3p[overhang_3p]
                        if overhang_5p != expected_overhang_5p:
                            error_message= f"*** ERROR *** OLIGO SEQUENCE in ASSEMBLY {name} [{id}] has mismatched 5p overhangs. Expected {expected_overhang_5p} but actual is {overhang_5p}"
                            print(error_message)
                            oligoPoolSpecification['failed_assemblies'].append( (assembly, error_message) )
                    else:
                        pairedOverhangs_5p[overhang_5p] = overhang_3p

                    if overhang_3p in pairedOverhangs_3p:
                        overhang_5p = pairedOverhangs_3p[overhang_3p]
                        expected_overhang_3p = pairedOverhangs_5p[overhang_5p]
                        if overhang_3p != expected_overhang_3p:
                            error_message = f"*** ERROR *** OLIGO SEQUENCE in ASSEMBLY {name} [{id}] has mismatched 3p overhangs. Expected {expected_overhang_3p} but actual is {overhang_3p}"
                            print(error_message)
                            oligoPoolSpecification['failed_assemblies'].append( (assembly, error_message) )
                    else:
                        pairedOverhangs_3p[overhang_3p] = overhang_5p

        return oligoPoolSpecification

    def calculateOligopoolCost(self, maximumOligoLength, numOligos, provider):

        # Calculate the oligopool synthesis cost (no discount applied)
        if provider == 'twist':
            oligoThresholds = [100, 500, 1000, 2000, 6000, 12000, 18000, 24000, 30000, 36000, 42000, 48000, 54000, 60000, 72000, 84000, 96000, 120000, 150000, 180000, 210000, 240000, 300000, 360000, 420000, 480000, 600000, 696000]
            if maximumOligoLength > 150 and maximumOligoLength <= 200:
                tierCosts = [520.0, 1040.0, 1560.0, 2080.0,3120.0, 4056.0, 5273.0, 6855.0, 8912.0, 9149.0, 10064.0, 11070.0, 12177.0, 13395.0, 14065.0, 14557.0, 14953.0, 15119.0, 17270.0, 18997.0, 20123.0, 22135.0, 27000.0, 32400.0, 37800.0, 43200.0, 52920.0, 62640.0]
                assert len(oligoThresholds) == len(tierCosts), "Uh Oh!"
                for (tier, threshold) in enumerate(oligoThresholds):
                    if numOligos <= threshold:
                        oligoCost = tierCosts[tier]
                        break

            elif maximumOligoLength > 200 and maximumOligoLength <= 250:
                tierCosts = [689.0, 1379.0, 2068.0, 2757.0, 4136.0, 5148.0, 6694.0, 8702.0, 11315.0, 11615.0, 12775.0, 14051.0, 15456.0, 17003.0, 17854.0, 18477.0, 18981.0, 19190.0, 21919.0, 24111.0, 25541.0, 28095.0, 33000.0, 39600.0, 46200.0, 52800.0, 64680.0, 76560.0]
                assert len(oligoThresholds) == len(tierCosts), "Uh Oh!"
                for (tier, threshold) in enumerate(oligoThresholds):
                    if numOligos <= threshold:
                        oligoCost = tierCosts[tier]
                        break

            elif maximumOligoLength > 250 and maximumOligoLength <= 300:
                tierCosts = [1030.0, 2060.0, 3090.0, 4121.0, 6181.0, 7694.0, 10004.0, 13006.0, 16910.0, 17359.0, 19093.0, 21000.0, 23100.0, 25411.0, 26684.0, 27614.0, 28367.0, 28680.0, 32758.0, 36035.0, 38172.0, 41989.0, 49320.0, 59184.0, 69048.0, 78912.0, 96667.0, 114422.0]
                assert len(oligoThresholds) == len(tierCosts), "Uh Oh!"
                for (tier, threshold) in enumerate(oligoThresholds):
                    if numOligos <= threshold:
                        oligoCost = tierCosts[tier]
                        break
            else:
                oligoCost = 0.0

        elif provider == 'genscript':
            if maximumOligoLength < 150:
                if numOligos < 12472:
                    oligoCost = 2200.0
                elif numOligos < 91766:
                    oligoCost = 5500.0
                else:
                    numChips = int(numOligos / 91766) + 1
                    oligoCost = 5500.0 * numChips
            elif maximumOligoLength > 150 and maximumOligoLength <= 170:
                if numOligos < 12472:
                    oligoCost = 2400.0
                elif numOligos < 91766:
                    oligoCost = 6000.0
                else:
                    numChips = int(numOligos / 91766) + 1
                    oligoCost = 6000.0 * numChips
            else:
                print("Genscript does not currently sell oligos over 170 nt long")
                oligoCost = None
        else:
            raise Exception("Currently, the oligopool provider must be either twist or genscript")

        return oligoCost

    def calculateCosts(self, oligoPoolSpecification, provider = 'twist'):
        primer_sequence_length = oligoPoolSpecification.get('primer_sequence_length', self.DEFAULT_PRIMER_SEQUENCE_LENGTH)

        costsPerPCR = 3.06  # water, Q5 DNA poly, buffer, plastic plate, PCR clean up well
        costsPerPrimerNt = 0.11 # IDT 500 pmole primer cost per nt
        costsPerAssembly = 5.89 # water, T4 DNA ligase, BsaI enzyme, PCR clean up well
        costsPerPrimerPair = costsPerPrimerNt * primer_sequence_length * 2

        assemblyList = oligoPoolSpecification['library_specifications']
        uniqueWells = list(set([spec['well_number'] for spec in assemblyList]))

        numAssemblies = len(uniqueWells)
        numPCRs = len(uniqueWells)

        costs = {'cost_per_PCR' : costsPerPCR, 'total_PCR_cost' : numPCRs * costsPerPCR, 'numAssemblies' : numAssemblies,
                 'cost_per_assembly' :  costsPerAssembly, 'total_assembly_cost' : numAssemblies * costsPerAssembly,
                }

        # Tally the number of oligos in the pool
        numOligos = 0
        maximumOligoLength = 0
        oligoSequences = []
        for assembly in assemblyList:
            primer_5p = assembly['primer_5p']
            primer_3p = assembly['primer_3p']
            fragmentList = assembly['fragment_list']
            maximumOligoLength = max(maximumOligoLength, max([len(fragment['oligo_sequence']) for fragment in fragmentList]))
            numOligos += len(fragmentList)
            oligoSequences += [fragment['oligo_sequence'] for fragment in fragmentList]
        numUniqueOligos = len(list(set(oligoSequences)))

        print('Unique # Oligos: {} vs. # Oligos: {}'.format(numUniqueOligos, numOligos))

        costs['num_oligos'] = numOligos
        costs['num_unique_oligos'] = numUniqueOligos
        costs['maximum_oligo_length'] = maximumOligoLength
        costs['provider'] = provider

        oligoCost = self.calculateOligopoolCost(maximumOligoLength, numOligos, provider)
        uniqueOligoCost = self.calculateOligopoolCost(maximumOligoLength, numUniqueOligos, provider)

        costs['total_cost_oligopool_synthesis'] = oligoCost
        costs['total_cost_unique_oligopool_synthesis'] = uniqueOligoCost
        costs['total_cost'] = costs['total_PCR_cost'] + costs['total_assembly_cost'] + costs['total_cost_oligopool_synthesis']
        costs['total_cost_unique_oligopool'] = costs['total_PCR_cost'] + costs['total_assembly_cost'] + costs['total_cost_unique_oligopool_synthesis']
        if numAssemblies > 1:
            costs['total_cost_per_assembly'] = costs['total_cost'] / numAssemblies
            costs['total_cost_unique_oligopool_per_assembly'] = costs['total_cost_unique_oligopool'] / numAssemblies
        else:
            costs['total_cost_per_assembly'] = 0.0
            costs['total_cost_unique_oligopool_per_assembly'] = 0.0

        oligoPoolSpecification['material_costs'] = costs

        return oligoPoolSpecification

    def initialize_screen(self, sequenceSpecList):
        base_url = "https://www.denovodna.com/screening"

        caller_id = random_id()

        for (n, item) in enumerate(sequenceSpecList):
            header = '> {} {}'.format(caller_id, n)
            seq = item['nucleotide_sequence'].upper()
            fasta = '{}\n{}\n'.format(header, seq)
            data = {
                "fasta": fasta,
                "region": "all",
                "provider_reference": "denovodna",
                }

            thread = ScreeningRequestThread(base_url, data)
            thread.start()
            self.screening_threads.append(thread)
            time.sleep(0.1)

    def cull_designs_by_screening_result(self, oligoPoolSpecification):
        culledResults = []
        assemblyList = oligoPoolSpecification['library_specifications']

        for thread in self.screening_threads:

            thread.join()
            (response, status_code) = thread.get_response()

            if status_code == 200:
                if response['synthesis_permission'] == 'denied':
                    for hit in response['hits_by_record']:
                        fasta_header = str(hit['fasta_header']).strip()
                        assemblyID = int(fasta_header.split(' ')[1])

                        culled_assembly = assemblyList[assemblyID]
                        culledResults.append({'name' : culled_assembly['name'], 'id' : culled_assembly['id']})
                        print('Hazardous Sequence Identified: {} [#{}]. Culling from assembly list.'.format(culled_assembly['name'], culled_assembly['id']))
            else:
                raise Exception('*** ERROR *** Could not complete DNA screening process. Server status_code was {} with response {}'.format(status_code, response))
        culledIDs = [culled['id'] for culled in culledResults]
        assemblyList = [assembly for assembly in assemblyList if not assembly['id'] in culledIDs]

        print('Hazardous Sequence Screening Complete. {} sequences were flagged for removal.'.format(len(culledResults)))
        oligoPoolSpecification['library_specifications'] = assemblyList
        oligoPoolSpecification['screening_results'] = culledResults

        self.screening_threads = []

        return oligoPoolSpecification

    def finalize(self, oligoPoolSpecification):

        if self.doScreening:
            print('** Screening Results **')
            oligoPoolSpecification = self.cull_designs_by_screening_result(oligoPoolSpecification)
        else:
            oligoPoolSpecification['screening_results'] = []

        print('** Finding Maximum Oligo Lengths **')
        oligoPoolSpecification = self.selectOligoLength(oligoPoolSpecification)

        print('** Designing Padding Sequences **')
        oligoPoolSpecification = self.designPadding(oligoPoolSpecification)

        print('** Designing Primers and Primer Binding Sites **')
        oligoPoolSpecification = self.designPrimers(oligoPoolSpecification)

        print('** Running Quality Control **')
        oligoPoolSpecification = self.qualityControl(oligoPoolSpecification)

        print('** Calculating Material Costs **')
        oligoPoolSpecification = self.calculateCosts(oligoPoolSpecification)

        return oligoPoolSpecification

    def print(self, oligoPoolSpecification):

        librarySpecs = oligoPoolSpecification['library_specifications']
        for assembly in librarySpecs:
            primer_5p = assembly['primer_5p']
            primer_3p = assembly['primer_3p']
            fragmentList = assembly['fragment_list']

            print('Assembly {}: Enzyme {}. Sequence Length: {} nt. {} oligos. Predicted ligation fidelity: {}'.format(assembly['id'], assembly['assembly_enzyme'], len(assembly['sequence']), len(assembly['fragment_list']), assembly['ligation_fidelity']))
            print('Primers: TM 5p {} / {} 3p. dG_fold 5p {} / {} 3p. dG_homodimer 5p {} / {} 3p. dG_primer_dimer 5p {} / {} 3p.'.format(round(primer_5p['Tm'],1),
                  round(primer_3p['Tm'],1), round(primer_5p['dG_folding'],2), round(primer_3p['dG_folding'],2), round(primer_5p['dG_homodimer'],2), round(primer_3p['dG_homodimer'],2),
                  round(primer_5p['dG_primer_dimer'],2), round(primer_3p['dG_primer_dimer'],2)))

            for (n, fragment) in enumerate(fragmentList):
                print('#oligo {}: [{}:{}] {}/{} fragment [{} nt]: {}'.format(n, fragment['begin'], fragment['end'], fragment['overhang_5p'], fragment['overhang_3p'], len(fragment['fragment_sequence']), fragment['fragment_sequence']))

    def convertCounterToWell(self, counter, plate = '96'):
        if plate == '96':
            plateNumber = int(counter / 96) + 1
            counter = counter % 96
            columns = ['A','B','C','D','E','F','G','H']
            numRows = 12
            numColumns = len(columns)
        elif plate == '384':
            plateNumber = int(counter / 384) + 1
            counter = counter % 384
            columns = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
            numRows = 24
            numColumns = len(columns)

        columnPosition = int(counter / numRows)
        row = counter % numRows
        well_name = '{}{}'.format(columns[columnPosition], row+1)
        return (well_name, plateNumber)

    def export(self, oligoPoolSpecification, excelFilename, oligoExcelFilename, primerPlateMapExcelFilename):

        print('*** EXPORT BEGINNING ***')

        inputData = []
        inputDataHeaders = ['name', 'well_number', 'assembly_id', 'enzyme', 'is_circular', 'nucleotide_sequence', 'sequence_length']
        overviewData = []
        overviewDataHeaders = ['name','assembly_id','plate_number','well_position','primer_Tm_5p','primer_Tm_3p',
                               'enzyme','first_overhang','last_overhang','is_circular',
                               'number_oligos','predicted_ligation_fidelity',
                               ]
        oligoData = []
        uniqueOligoData = []
        oligoDataHeaders = ['name','sequence']

        assemblyData = []
        assemblyDataHeaders = ['assembly_id','name','enzyme','oligo_number','overhang_5p','overhang_3p','begin','end','region','fragment_with_cut_sites','padding_sequence','oligo_sequence']

        ## Prepare Primer Data and Plate Maps for Primers
        primerData = []
        primerDataHeaders = ['assembly_id', 'plate_number','well_position', 'primer_5p_sequence', 'primer_3p_sequence',
                             'TM_5p', 'TM_3p', 'dG_fold_5p','dG_fold_3p','dG_homodimer_5p','dG_homodimer_3p', 'dG_primer_dimer',
                             'structure_5p','structure_3p','homodimer_structure_5p','homodimer_structure_3p', 'structure_primer_dimer'
                             ]

        primerPlateMapData_5p = []
        primerPlateMapData_3p = []
        primerPlateMapHeaders = ['well_position', 'name', 'sequence']

        assemblyList = oligoPoolSpecification['library_specifications']
        uniqueWells = list(set([spec['well_number'] for spec in oligoPoolSpecification['library_specifications']]))

        for (counter, well) in enumerate(uniqueWells):
            assemblyList = [assembly for assembly in oligoPoolSpecification['library_specifications'] if assembly['well_number'] == well]

            for assembly in assemblyList:
                primer_5p = assembly['primer_5p']
                primer_3p = assembly['primer_3p']
                fragmentList = assembly['fragment_list']
                (well_name, plateNumber) = self.convertCounterToWell(counter)

                inputData.append({'name' : assembly['name'], 'well_number' : assembly['well_number'], 'assembly_id' : assembly['id'], 'enzyme' : assembly['assembly_enzyme'], 'is_circular' : assembly['is_circular'], 'nucleotide_sequence' : assembly['sequence'], 'sequence_length' : len(assembly['sequence']) })
                overviewData.append({'name' : assembly['name'], 'assembly_id' : assembly['id'], 'well_position' : well_name, 'plate_number' : plateNumber,
                                     'primer_Tm_5p' : round(primer_5p['Tm'],2), 'primer_Tm_3p' : round(primer_3p['Tm'],2),
                                     'enzyme' : assembly['assembly_enzyme'], 'first_overhang' : assembly['first_overhang'], 'last_overhang' : assembly['last_overhang'],
                                     'is_circular' : assembly['is_circular'], 'number_oligos' : len(fragmentList), 'predicted_ligation_fidelity' : assembly['ligation_fidelity'],
                                     })

                for (oligoNum, fragment) in enumerate(fragmentList):
                    uniqueName = 'A{}-O{}'.format(assembly['id'], oligoNum)

                    oligoData.append({'name' : uniqueName, 'sequence' : fragment['oligo_sequence']})
                    if not (fragment['oligo_sequence'] in [x['sequence'] for x in uniqueOligoData]): uniqueOligoData.append({'name' : uniqueName, 'sequence' : fragment['oligo_sequence']})

                    assemblyData.append({'assembly_id' : assembly['id'], 'name' : assembly['name'], 'enzyme' : assembly['assembly_enzyme'], 'oligo_number' : oligoNum,
                                         'overhang_5p' : fragment['overhang_5p'], 'overhang_3p' : fragment['overhang_3p'],
                                         'begin' : fragment['begin'], 'end' : fragment['end'], 'region' : fragment['fragment_sequence'],
                                         'fragment_with_cut_sites' : fragment['fragment_with_cut_sites'], 'padding_sequence' : fragment['padding_sequence'],
                                         'oligo_sequence' : fragment['oligo_sequence'],
                                    })

            # Find first assembly carried out in each unique well (if each well is one assembly, then no difference)
            assembly = assemblyList[0]

            primer_5p = assembly['primer_5p']
            primer_3p = assembly['primer_3p']

            primerData.append({'assembly_id' : assembly['id'], 'primer_5p_sequence' : primer_5p['sequence'], 'primer_3p_sequence' : primer_3p['sequence'],
                               'well_position' : well_name, 'plate_number' : plateNumber,
                               'TM_5p' : round(primer_5p['Tm'],2), 'TM_3p' : round(primer_3p['Tm'],2),
                               'dG_fold_5p' : round(primer_5p['dG_folding'],2), 'dG_fold_3p' : round(primer_3p['dG_folding'],2),
                               'dG_homodimer_5p' : round(primer_5p['dG_homodimer'],2), 'dG_homodimer_3p' : round(primer_3p['dG_homodimer'],2),
                               'dG_primer_dimer' : round(primer_5p['dG_primer_dimer'],2),
                               'structure_5p' : primer_5p['structure'], 'structure_3p' : primer_3p['structure'],
                               'homodimer_structure_5p' : primer_5p['homodimer_structure'], 'homodimer_structure_3p' : primer_3p['homodimer_structure'],
                               'structure_primer_dimer' : primer_5p['structure_primer_dimer']
                               })

            primerPlateMapData_5p.append({'plate_number' : plateNumber, 'well_position' : well_name, 'name' : 'p5_' + assembly['name'], 'sequence' : primer_5p['sequence']})
            primerPlateMapData_3p.append({'plate_number' : plateNumber, 'well_position' : well_name, 'name' : 'p3_' + assembly['name'], 'sequence' : primer_3p['sequence']})

        materialCosts = [oligoPoolSpecification['material_costs']]
        materialCostHeaders = ['cost_per_PCR', 'cost_per_assembly', 'numAssemblies', 'total_cost_per_PCR', 'total_cost_per_assembly',
                               'num_oligos', 'maximum_oligo_length', 'provider', 'total_cost_oligopool_synthesis', 'total_cost', 'total_cost_per_assembly',
                               'num_unique_oligos', 'total_cost_unique_oligopool_synthesis', 'total_cost_unique_oligopool', 'total_cost_unique_oligopool_per_assembly',
                               ]

        screening_results = [oligoPoolSpecification['screening_results']]
        screeningHeaders = ['id','name']


        df_inputData = pd.DataFrame(inputData, columns = inputDataHeaders)
        df_overviewData = pd.DataFrame(overviewData, columns=overviewDataHeaders)
        df_oligoData = pd.DataFrame(oligoData, columns=oligoDataHeaders)
        df_uniqueOligoData = pd.DataFrame(uniqueOligoData, columns=oligoDataHeaders)
        df_primerData = pd.DataFrame(primerData, columns=primerDataHeaders)
        df_assemblyData = pd.DataFrame(assemblyData, columns=assemblyDataHeaders)
        df_materialCosts = pd.DataFrame(materialCosts, columns=materialCostHeaders)
        if len(oligoPoolSpecification['screening_results']) > 0: df_screening = pd.DataFrame(screening_results, columns=screeningHeaders)

        writer = pd.ExcelWriter(excelFilename, engine='openpyxl')
        df_inputData.to_excel(writer, sheet_name='inputs', header=inputDataHeaders)
        df_overviewData.to_excel(writer, sheet_name='overview', header=overviewDataHeaders)
        df_oligoData.to_excel(writer, sheet_name='oligos', header=oligoDataHeaders)
        df_uniqueOligoData.to_excel(writer, sheet_name='unique_oligos', header=oligoDataHeaders)
        df_primerData.to_excel(writer, sheet_name='primers', header=primerDataHeaders)
        df_assemblyData.to_excel(writer, sheet_name='assemblies', header=assemblyDataHeaders)
        df_materialCosts.to_excel(writer, sheet_name='material costs', header=materialCostHeaders)
        if len(oligoPoolSpecification['screening_results']) > 0: df_screening.to_excel(writer, sheet_name='screening results', header=screeningHeaders)

        writer.close()
        print('*** MASTER EXCEL FILE EXPORTED TO {}'.format(excelFilename))

        unique_oligo_filename = ".".join(oligoExcelFilename.split('.')[:-1]) + '_unique.xlsx'

        df_oligoData.to_excel(oligoExcelFilename, index=False, header=oligoDataHeaders)
        df_uniqueOligoData.to_excel(unique_oligo_filename, index=False, header=oligoDataHeaders)
        print('*** OLIGO EXCEL FILE EXPORTED TO {}'.format(oligoExcelFilename))
        print('*** UNIQUE OLIGO EXCEL FILE EXPORTED TO {}'.format(unique_oligo_filename))

        df_primerPlateMapData_5p = pd.DataFrame(primerPlateMapData_5p)
        df_primerPlateMapData_3p = pd.DataFrame(primerPlateMapData_3p)

        unique_plate_numbers = df_primerPlateMapData_5p['plate_number'].unique()
        for plate_number in unique_plate_numbers:
            sliced_df_primerPlateMapData_5p = df_primerPlateMapData_5p[df_primerPlateMapData_5p['plate_number'] == plate_number]
            sliced_df_primerPlateMapData_3p = df_primerPlateMapData_3p[df_primerPlateMapData_3p['plate_number'] == plate_number]

            prefix = primerPlateMapExcelFilename.split('.')[0]
            output_file_5p = "{}_5p_{}.xls".format(prefix, plate_number)
            output_file_3p = "{}_3p_{}.xls".format(prefix, plate_number)
            sliced_df_primerPlateMapData_5p.to_excel(output_file_5p, sheet_name='Plate{}_5p'.format(plate_number), index=False, columns=primerPlateMapHeaders, header=primerPlateMapHeaders, engine='openpyxl')
            sliced_df_primerPlateMapData_3p.to_excel(output_file_3p, sheet_name='Plate{}_3p'.format(plate_number), index=False, columns=primerPlateMapHeaders, header=primerPlateMapHeaders, engine='openpyxl')
            print('*** PRIMER PLATE MAP EXCEL FILE EXPORTED TO {}'.format(output_file_5p))
            print('*** PRIMER PLATE MAP EXCEL FILE EXPORTED TO {}'.format(output_file_3p))
        print('*** EXPORT COMPLETED ***')

        output = {'input_data' : inputData, 'overview_data' : overviewData, 'oligo_data' : oligoData, 'primer_data' : primerData,
                  'assembly_data' : assemblyData, 'material_costs' : materialCosts,
                  'primer_plate_maps_5p' : primerPlateMapData_5p, 'primer_plate_maps_3p' : primerPlateMapData_3p,
                  'combination_set_to_wells' : self.combo_set_to_wells, 'well_to_combination_set' : self.well_to_combo_set,
                  'well_to_assembly_names' : self.well_to_assembly_names,
                }
        return output

    def output(self, oligoPoolSpecification):

        inputData = []
        overviewData = []
        oligoData = []
        uniqueOligoData = []
        primerData = []
        assemblyData = []

        assemblyList = oligoPoolSpecification['library_specifications']
        uniqueWells = list(set([spec['well_number'] for spec in oligoPoolSpecification['library_specifications']]))

        for (counter, well) in enumerate(uniqueWells):
            assemblyList = [assembly for assembly in oligoPoolSpecification['library_specifications'] if assembly['well_number'] == well]
            for assembly in assemblyList:

                primer_5p = assembly['primer_5p']
                primer_3p = assembly['primer_3p']
                fragmentList = assembly['fragment_list']

                inputData.append({'name' : assembly['name'], 'assembly_id' : assembly['id'], 'enzyme' : assembly['assembly_enzyme'], 'is_circular' : assembly['is_circular'], 'nucleotide_sequence' : assembly['sequence'], 'sequence_length' : len(assembly['sequence']) })
                overviewData.append({'name' : assembly['name'], 'assembly_id' : assembly['id'],
                                     'primer_Tm_5p' : primer_5p['Tm'], 'primer_Tm_3p' : primer_3p['Tm'],
                                     'enzyme' : assembly['assembly_enzyme'], 'first_overhang' : assembly['first_overhang'], 'last_overhang' : assembly['last_overhang'],
                                     'is_circular' : assembly['is_circular'], 'number_oligos' : len(fragmentList), 'predicted_ligation_fidelity' : assembly['ligation_fidelity'],
                                     })

                for (oligoNum, fragment) in enumerate(fragmentList):
                    uniqueName = 'A{}-O{}'.format(assembly['id'], oligoNum)

                    oligoData.append({'name' : uniqueName, 'sequence' : fragment['oligo_sequence']})
                    if not (fragment['oligo_sequence'] in [x['sequence'] for x in uniqueOligoData]): uniqueOligoData.append({'name' : uniqueName, 'sequence' : fragment['oligo_sequence']})

                    assemblyData.append({'assembly_id' : assembly['id'], 'name' : assembly['name'], 'enzyme' : assembly['assembly_enzyme'], 'oligo_number' : oligoNum,
                                         'overhang_5p' : fragment['overhang_5p'], 'overhang_3p' : fragment['overhang_3p'],
                                         'begin' : fragment['begin'], 'end' : fragment['end'], 'region' : fragment['fragment_sequence'],
                                         'fragment_with_cut_sites' : fragment['fragment_with_cut_sites'], 'padding_sequence' : fragment['padding_sequence'],
                                         'oligo_sequence' : fragment['oligo_sequence'],
                                        })

            # Find first assembly carried out in each unique well (if each well is one assembly, then no difference)
            assembly = assemblyList[0]

            if len(assemblyList) > 1:
                lastAssembly = assemblyList[-1]
                assembly_id = '{}:{}'.format(assembly['id'], lastAssembly['id'])
            else:
                assembly_id = str(assembly['id'])

            primer_5p = assembly['primer_5p']
            primer_3p = assembly['primer_3p']

            primerData.append({'assembly_id' : assembly_id, 'primer_5p_sequence' : primer_5p['sequence'], 'primer_3p_sequence' : primer_3p['sequence'],
                               'TM_5p' : round(primer_5p['Tm'],2), 'TM_3p' : round(primer_3p['Tm'],2),
                               'dG_fold_5p' : round(primer_5p['dG_folding'],2), 'dG_fold_3p' : round(primer_3p['dG_folding'],2),
                               'dG_homodimer_5p' : round(primer_5p['dG_homodimer'],2), 'dG_homodimer_3p' : round(primer_3p['dG_homodimer'],2),
                               'dG_primer_dimer' : round(primer_5p['dG_primer_dimer'],2),
                               'structure_5p' : primer_5p['structure'], 'structure_3p' : primer_3p['structure'],
                               'homodimer_structure_5p' : primer_5p['homodimer_structure'], 'homodimer_structure_3p' : primer_3p['homodimer_structure'],
                               'structure_primer_dimer' : primer_5p['structure_primer_dimer']
                               })

        materialCosts = [oligoPoolSpecification['material_costs']]
        screening_results = oligoPoolSpecification['screening_results']

        output = {'input_data' : inputData, 'overview_data' : overviewData, 'oligo_data' : oligoData, 'primer_data' : primerData,
                  'assembly_data' : assemblyData, 'material_costs' : materialCosts, 'screening_results' : screening_results,
                  'failed_assemblies' : oligoPoolSpecification['failed_assemblies'],
                }
        return output

def GibsonAssembly(DNA_sequence, MAX_FRAGMENT_LENGTH = 1500, MIN_FRAGMENT_LENGTH = 400, FRAGMENT_OVERLAP = 30, PRIMER_LENGTH = 18, TM_TARGET = [50.0, 65.0], excluded_motif_list = ['CCCCC','GGGGG','AAAAA','TTTTT']):

    #Tm settings optimized to match NEB TM Calculator for Q5 Polymerase to within 0.7 degC (on average)
    Tm_settings = {'dnac1' : 500.0, 'Na' : 150.0, 'K' : 0.0, 'Mg' : 2.0, 'dNTPs' : 0.2, 'Tris' : 65.0, 'saltcorr' : 1} #Q5 DNA polymerase buffer
    energy_model = PyVRNA(dangles=2, gquad=False, parameter_file='dna_mathews2004.par')

    toleranceDeltaT = 0.20 #degC
    tolerance_dG_folding = -1.5 #kcal/mol
    tolerance_dG_homodimer = -8.0 #kcal/mol
    tolerance_primer_dimer_folding_energy = -5.0 #kcal/mol

    seq = DNA_sequence
    seqlen = len(seq)

    ## Step 1:  Evaluate candidate primer binding sites and identify valid sites.
    acceptablePrimerList = []
    candidatePrimerList = [{'pos' : n, 'seq' : seq[n : n + PRIMER_LENGTH], 'revcomp_seq' : _revcomp(seq[n : n + PRIMER_LENGTH]) } for n in range(seqlen-PRIMER_LENGTH)]
    for candidate in candidatePrimerList:
        if all([candidate['seq'].find(motif) == -1 for motif in excluded_motif_list]):
            Tm = mt.Tm_NN(candidate['seq'], Tm_settings)
            if Tm >= TM_TARGET[0] and Tm <= TM_TARGET[1]:
                candidate['Tm'] = Tm
                primer_fold_fwd = energy_model.RNAfold(candidate['seq'])
                primer_fold_rev = energy_model.RNAfold(candidate['revcomp_seq'])
                if primer_fold_fwd.energy >= tolerance_dG_folding and primer_fold_rev.energy >= tolerance_dG_folding:
                    candidate['dG_folding_fwd'] = primer_fold_fwd.energy
                    candidate['dG_folding_rev'] = primer_fold_rev.energy
                    dimer_fold_fwd = energy_model.RNAcofold(sequences = [ candidate['seq'] , candidate['seq']])
                    dimer_fold_rev = energy_model.RNAcofold(sequences = [ candidate['revcomp_seq'] , candidate['revcomp_seq']])
                    if dimer_fold_fwd.energy >= tolerance_dG_homodimer and dimer_fold_rev.energy >= tolerance_dG_homodimer:
                        candidate['dG_homodimer_fwd'] = dimer_fold_fwd.energy
                        candidate['dG_homodimer_rev'] = dimer_fold_rev.energy
                        acceptablePrimerList.append( candidate )

    ## Step 2:  Identify pairs of compatible (similar Tm) primer binding sites that are (MAX_FRAGMENT_LENGTH-PRIMER_LENGTH) apart
    ##          fragments must have a length of at least MIN_FRAGMENT_LENGTH
    cutoff = [MIN_FRAGMENT_LENGTH - PRIMER_LENGTH, MAX_FRAGMENT_LENGTH - PRIMER_LENGTH]
    numFragments = 1000000
    foundPrimers = False

    for target_Tm in np.arange(TM_TARGET[1], TM_TARGET[0], -0.1):
        selectedPrimers = []
        equalMeltPrimers = [primer for primer in acceptablePrimerList if (abs(primer['Tm'] - target_Tm) <= toleranceDeltaT)]
        if len(equalMeltPrimers) == 0: continue

        firstPrimer = equalMeltPrimers[0]
        leftPrimer = firstPrimer

        selectedPrimers.append(firstPrimer)
        while True:
            rightPrimers = [primer for primer in equalMeltPrimers if ((primer['pos'] - leftPrimer['pos']) <= cutoff[1]) and ((primer['pos'] - leftPrimer['pos']) > cutoff[0])]
            sortedRightPrimers = sorted(rightPrimers, key = lambda x: x['pos'] - leftPrimer['pos'], reverse = True)

            if len(sortedRightPrimers) > 0:
                selectedRightPrimer = sortedRightPrimers[0]
                selectedPrimers.append(selectedRightPrimer)
                leftPrimer = selectedRightPrimer

            else:
                foundSolution = False
                break

            if (seqlen - selectedRightPrimer['pos'] + firstPrimer['pos']) < cutoff[1]:
                foundSolution = True
                break

        if foundSolution and len(selectedPrimers) < numFragments:
            numFragments = len(selectedPrimers)
            finalSelectedPrimers = selectedPrimers
            foundPrimers = True

    ## Step 3:  Generate DNA FRAGMENTS and PRIMERs based on selected primer binding site positions and overlap lengths
    if not foundPrimers: raise Exception("Error:  The TM_TARGET range was too narrow or too high to generate valid primer binding sites. Lower the TM_TARGET temperatures to resolve this error.")
    fragmentList = []
    primerList = []

    for (n, PBS_fwd) in enumerate(finalSelectedPrimers):

        if (n + 1) < len(finalSelectedPrimers):
            PBS_rev = finalSelectedPrimers[n+1]
            fragment = {'name' : 'F_{}'.format(n+1), 'position' : [PBS_fwd['pos'], PBS_rev['pos'] + PRIMER_LENGTH], 'sequence' : seq[PBS_fwd['pos'] : PBS_rev['pos'] + PRIMER_LENGTH]}

        else: #wrap around circular plasmid
            PBS_rev = finalSelectedPrimers[0]
            fragment = {'name' : 'F_{}'.format(n+1), 'position' : [PBS_fwd['pos'], PBS_rev['pos'] + PRIMER_LENGTH], 'sequence' : seq[PBS_fwd['pos'] : seqlen] + seq[0 : PBS_rev['pos'] + PRIMER_LENGTH]}

        primer_begin = PBS_fwd['pos'] - (FRAGMENT_OVERLAP - PRIMER_LENGTH)
        primer_end = PBS_fwd['pos'] + PRIMER_LENGTH

        if primer_begin > 0:
            primer_fwd = seq[primer_begin: primer_end]
        else:
            primer_fwd = seq[primer_begin : ] + seq[0 : primer_end]


        rightBorderPos = PBS_rev['pos'] + FRAGMENT_OVERLAP
        if (rightBorderPos < seqlen):
            primer_rev = _revcomp(seq[ PBS_rev['pos'] : rightBorderPos])
        else:
            primer_rev = _revcomp(seq[ PBS_rev['pos'] : seqlen] + seq[0 : (seqlen - rightBorderPos)])

        selectedFwdPrimer = {'name' : 'P_{}_fwd'.format(n+1), 'sequence' : primer_fwd, 'position' : [primer_begin, primer_end], 'Tm' : PBS_fwd['Tm'], 'dG_folding' : PBS_fwd['dG_folding_fwd'], 'dG_homodimer' : PBS_fwd['dG_homodimer_fwd'], 'strand' : +1 }
        selectedRevPrimer = {'name' : 'P_{}_rev'.format(n+1), 'sequence' : primer_rev, 'position' : [PBS_rev['pos'], rightBorderPos], 'Tm' : PBS_rev['Tm'], 'dG_folding' : PBS_rev['dG_folding_rev'], 'dG_homodimer' : PBS_rev['dG_homodimer_rev'], 'strand' : -1 }

        fragmentList.append( fragment )
        primerList.append(selectedFwdPrimer)
        primerList.append(selectedRevPrimer)

    output = {'fragmentList' : fragmentList, 'primerList' : primerList }
    return output

def parse_comma_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [item.strip() for item in str(value).split(',') if item.strip()]

def parse_optional_none(value):
    if value is None:
        return None
    if isinstance(value, str) and value.strip().lower() in ('', 'none', 'null'):
        return None
    return value

def parse_bool(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)):
        return bool(value)
    if value is None or pd.isna(value):
        return False
    normalized = str(value).strip().lower()
    if normalized in ('1', 'true', 't', 'yes', 'y'):
        return True
    if normalized in ('0', 'false', 'f', 'no', 'n', ''):
        return False
    raise ValueError(f"Invalid boolean value: {value}")

def normalize_build_enzyme(enzyme):
    enzyme = parse_optional_none(enzyme)
    if enzyme is None:
        return None
    enzyme = str(enzyme).strip()
    if enzyme in GoldenGatePool.CUT_SITE_SEQUENCES:
        return enzyme
    lower_enzyme = enzyme.lower().replace('-', '').replace('_', '')
    if lower_enzyme in GoldenGatePool.NEB_ENZYME_NAMES:
        return GoldenGatePool.NEB_ENZYME_NAMES[lower_enzyme]
    raise ValueError(
        "Unsupported assembly enzyme '{}'. Supported enzymes: {}".format(
            enzyme, ", ".join(sorted(GoldenGatePool.CUT_SITE_SEQUENCES.keys()))
        )
    )

def parse_no_overhang_ranges(value):
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == '' or stripped.lower() in ('none', 'null'):
            return None
        try:
            value = ast.literal_eval(stripped)
        except (ValueError, SyntaxError) as exc:
            raise ValueError(f"Invalid no_overhangs_position_range_list: {value}") from exc
    if value is None:
        return None
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"Invalid no_overhangs_position_range_list: {value}")

    ranges = []
    for item in value:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            raise ValueError(f"Invalid no_overhangs_position_range_list item: {item}")
        begin, end = int(item[0]), int(item[1])
        if begin < 0 or end < begin:
            raise ValueError(f"Invalid no-overhang range: {item}")
        ranges.append([begin, end])
    return ranges

def validate_nucleotide_sequence(sequence, name):
    sequence = str(sequence).strip().upper().replace(' ', '').replace('\n', '').replace('\r', '')
    if not sequence:
        raise ValueError(f"Sequence '{name}' is empty.")
    invalid_chars = sorted(set(sequence) - set('ACGT'))
    if invalid_chars:
        raise ValueError(
            "Sequence '{}' contains non-ACGT characters: {}".format(
                name, ''.join(invalid_chars)
            )
        )
    return sequence

def normalize_build_sequence_spec(spec, index, default_assembly_enzyme='BsaI_HFv2', default_level_two_assembly_enzyme='BbsI_HF', default_is_circular=0):
    if not isinstance(spec, dict):
        raise ValueError("JSON sequence specification item #{} must be an object.".format(index))

    if "name" not in spec:
        raise ValueError("Missing required key 'name' in sequence specification item #{}.".format(index))
    if "nucleotide_sequence" in spec:
        sequence_value = spec["nucleotide_sequence"]
    elif "sequence" in spec:
        sequence_value = spec["sequence"]
    else:
        raise ValueError("Missing required key 'nucleotide_sequence' or 'sequence' in sequence specification item #{}.".format(index))

    name = str(spec["name"]).strip()
    nucleotide_sequence = validate_nucleotide_sequence(sequence_value, name)

    assembly_enzyme = normalize_build_enzyme(spec.get("assembly_enzyme", spec.get("enzyme", default_assembly_enzyme)))
    if assembly_enzyme is None:
        raise ValueError("assembly_enzyme cannot be none for sequence '{}'.".format(name))

    level_two_value = spec.get("level_two_assembly_enzyme", default_level_two_assembly_enzyme)
    level_two_assembly_enzyme = normalize_build_enzyme(level_two_value) if level_two_value is not None else None

    return {"name": name,
            "assembly_enzyme": assembly_enzyme,
            "level_two_assembly_enzyme": level_two_assembly_enzyme,
            "combinatorial_set": int(spec["combinatorial_set"]) if "combinatorial_set" in spec and spec["combinatorial_set"] is not None else index,
            "is_circular": int(parse_bool(spec["is_circular"])) if "is_circular" in spec and spec["is_circular"] is not None else int(parse_bool(default_is_circular)),
            "no_overhangs_position_range_list" : parse_no_overhang_ranges(spec.get("no_overhangs_position_range_list", None)),
            "nucleotide_sequence": nucleotide_sequence}

def loadSequencesFromJsonFile(filename, default_assembly_enzyme='BsaI_HFv2', default_level_two_assembly_enzyme='BbsI_HF', default_is_circular=0):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    with open(filename, 'r') as handle:
        payload = json.load(handle)
    if not isinstance(payload, list):
        raise ValueError("JSON build input must contain a list of sequence specification objects.")
    return [normalize_build_sequence_spec(spec,
                                          index,
                                          default_assembly_enzyme=default_assembly_enzyme,
                                          default_level_two_assembly_enzyme=default_level_two_assembly_enzyme,
                                          default_is_circular=default_is_circular)
            for (index, spec) in enumerate(payload)]

def loadSequencesFromFile(filename, file_format='auto', sheet_name=0, default_assembly_enzyme='BsaI_HFv2', default_level_two_assembly_enzyme='BbsI_HF', default_is_circular=0):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    # --------------------------------------------------
    # Load file
    # --------------------------------------------------
    ext = os.path.splitext(filename)[1].lower()
    normalized_format = (file_format or 'auto').lower()
    if normalized_format == 'auto':
        if ext in [".csv", ".txt"]:
            normalized_format = 'csv'
        elif ext in [".xlsx", ".xls"]:
            normalized_format = 'excel'
        elif ext in [".fa", ".fas", ".fasta"]:
            normalized_format = 'fasta'
        elif ext in [".gb", ".gbk", ".genbank"]:
            normalized_format = 'genbank'
        elif ext in [".json"]:
            normalized_format = 'json'
        else:
            raise ValueError("Unsupported file type. Must be CSV, Excel, FASTA, GENBANK, or JSON format.")

    if normalized_format == 'csv':
        df = pd.read_csv(filename)
        df = df.dropna(how='all')
        df = df.convert_dtypes()
    elif normalized_format == 'excel':
        df = pd.read_excel(filename, sheet_name=sheet_name)
    elif normalized_format == 'fasta':
        with open(filename, 'r') as fc:
            records = SeqIO.parse(fc, 'fasta')
            seqList = [{'name' : str(record.id), 'sequence' : str(record.seq)} for record in records]
        df = pd.DataFrame(seqList)
    elif normalized_format == 'genbank':
        with open(filename, 'r') as fc:
            records = SeqIO.parse(fc, 'genbank')
            seqList = [{'name' : str(record.id), 'sequence' : str(record.seq)} for record in records]
        df = pd.DataFrame(seqList)
    elif normalized_format == 'json':
        return loadSequencesFromJsonFile(filename,
                                         default_assembly_enzyme=default_assembly_enzyme,
                                         default_level_two_assembly_enzyme=default_level_two_assembly_enzyme,
                                         default_is_circular=default_is_circular)
    else:
        raise ValueError("Unsupported file format. Must be auto, csv, excel, fasta, genbank, or json.")

    if df.empty:
        return []

    # --------------------------------------------------
    # Normalize column names
    # --------------------------------------------------
    df.columns = [str(c).strip().lower() for c in df.columns]

    # --------------------------------------------------
    # Validate required columns
    # --------------------------------------------------
    if "name" not in df.columns:
        raise ValueError("Missing required column: 'name'")

    if "sequence" in df.columns:
        seq_col = "sequence"
    elif "nucleotide_sequence" in df.columns:
        seq_col = "nucleotide_sequence"
    else:
        raise ValueError("Missing required column: 'sequence' or 'nucleotide_sequence'")

    # --------------------------------------------------
    # Optional columns with defaults
    # --------------------------------------------------
    has_comb_set = "combinatorial_set" in df.columns
    has_is_circ = "is_circular" in df.columns
    has_assembly_enzyme = "assembly_enzyme" in df.columns
    has_enzyme = "enzyme" in df.columns
    has_no_overhangs_position_range_list = "no_overhangs_position_range_list" in df.columns

    sequences = []

    for n, row in df.iterrows():
        name = str(row["name"]).strip()
        nucleotide_sequence = validate_nucleotide_sequence(row[seq_col], name)
        if has_assembly_enzyme and pd.notna(row["assembly_enzyme"]):
            assembly_enzyme = normalize_build_enzyme(row["assembly_enzyme"])
        elif has_enzyme and pd.notna(row["enzyme"]):
            assembly_enzyme = normalize_build_enzyme(row["enzyme"])
        else:
            assembly_enzyme = normalize_build_enzyme(default_assembly_enzyme)

        entry = {
            "name": name,
            "assembly_enzyme": assembly_enzyme,
            "level_two_assembly_enzyme": normalize_build_enzyme(default_level_two_assembly_enzyme) if default_level_two_assembly_enzyme is not None else None,
            "combinatorial_set": int(row["combinatorial_set"]) if has_comb_set and pd.notna(row["combinatorial_set"]) else n,
            "is_circular": int(parse_bool(row["is_circular"])) if has_is_circ and pd.notna(row["is_circular"]) else int(parse_bool(default_is_circular)),
            "no_overhangs_position_range_list" : parse_no_overhang_ranges(row['no_overhangs_position_range_list']) if has_no_overhangs_position_range_list else None,
            "nucleotide_sequence": nucleotide_sequence
        }

        sequences.append(entry)

    return sequences

def loadBuildJsonForExport(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    with open(filename, 'r') as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("GSC build JSON must decode to a JSON object.")
    if 'oligo_pool_specification' not in payload:
        raise ValueError("GSC build JSON is missing required key 'oligo_pool_specification'.")
    oligoPoolSpecification = payload['oligo_pool_specification']
    if not isinstance(oligoPoolSpecification, dict):
        raise ValueError("'oligo_pool_specification' must be a JSON object.")
    if 'library_specifications' not in oligoPoolSpecification:
        raise ValueError("'oligo_pool_specification' is missing required key 'library_specifications'.")
    if not isinstance(oligoPoolSpecification['library_specifications'], list):
        raise ValueError("'library_specifications' must be a list.")
    if len(oligoPoolSpecification['library_specifications']) == 0:
        raise ValueError("'library_specifications' is empty; no build specifications can be exported.")
    oligoPoolSpecification.setdefault('material_costs', {})
    oligoPoolSpecification.setdefault('screening_results', [])
    return (payload, oligoPoolSpecification)

def exportBuildJsonExcelFiles(build_json_filename, output_prefix, verbose = False):
    (payload, oligoPoolSpecification) = loadBuildJsonForExport(build_json_filename)
    build = GoldenGatePool(verbose=verbose, parallel=False, doScreening=False)
    exportData = build.export(
        oligoPoolSpecification,
        output_prefix + '.xlsx',
        output_prefix + '_oligos.xlsx',
        output_prefix + '_primer_plates.xlsx'
    )
    return {'run_id' : payload.get('run_id'),
            'num_sequence_specs' : len(payload.get('sequence_spec_list', [])) if isinstance(payload.get('sequence_spec_list'), list) else None,
            'num_library_specifications' : len(oligoPoolSpecification.get('library_specifications', [])),
            'master_excel_filename' : output_prefix + '.xlsx',
            'oligo_excel_filename' : output_prefix + '_oligos.xlsx',
            'unique_oligo_excel_filename' : output_prefix + '_oligos_unique.xlsx',
            'primer_plate_map_prefix' : output_prefix + '_primer_plates',
            'export_data' : exportData}

def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Build Golden Gate oligo pools from CSV, Excel, FASTA, GenBank, or JSON sequence specifications."
    )
    parser.add_argument('input', nargs='?', help='Input sequence specification file.')
    parser.add_argument('--json-input', default=None, help='JSON input file containing a list of build specification dictionaries.')
    parser.add_argument('--export-build-json', '--build-json', dest='export_build_json', default=None,
                        help='Parse a GSC-exported build.json and write Excel build output files without re-running assembly design.')
    parser.add_argument('-o', '--output-prefix', required=True, help='Output filename prefix.')
    parser.add_argument('--file-format', choices=['auto', 'csv', 'excel', 'fasta', 'genbank', 'json'], default='auto')
    parser.add_argument('--sheet-name', default=0, help='Excel sheet name or zero-based index. Default: 0.')
    parser.add_argument('--default-is-circular', default='0', help='Default circular flag for rows/records without is_circular.')

    parser.add_argument('--pcr-tm-target', type=float, default=60.0)
    parser.add_argument('--excluded-restriction-enzymes', default='BsaI,BbsI,Esp3I')
    parser.add_argument('--maximum-oligo-length', type=int, default=300)
    parser.add_argument('--primer-sequence-length', type=int, default=18)
    parser.add_argument('--padding-end', choices=['5p', '3p'], default='5p')
    parser.add_argument('--assembly-enzyme', default='BsaI_HFv2')
    parser.add_argument('--level-two-assembly-enzyme', default='BbsI_HF')

    parser.add_argument('--minimum-fragment-length', type=int, default=GoldenGatePool.MINIMUM_FRAGMENT_LENGTH)
    parser.add_argument('--long-assembly-threshold', type=int, default=GoldenGatePool.LONG_ASSEMBLY_THRESHOLD)
    parser.add_argument('--long-assembly-maximum-fragment-length', type=int, default=GoldenGatePool.LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH)
    parser.add_argument('--long-assembly-minimum-fragment-length', type=int, default=GoldenGatePool.LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH)
    parser.add_argument('--num-fragments-per-long-assembly', type=int, default=GoldenGatePool.NUM_FRAGMENTS_PER_LONG_ASSEMBLY)

    parser.add_argument('--tolerance-delta-tm', type=float, default=GoldenGatePool.DEFAULT_TOLERANCE_DELTA_TM)
    parser.add_argument('--tolerance-dg-folding', type=float, default=GoldenGatePool.DEFAULT_TOLERANCE_DG_FOLDING)
    parser.add_argument('--tolerance-dg-homodimer', type=float, default=GoldenGatePool.DEFAULT_TOLERANCE_DG_HOMODIMER)
    parser.add_argument('--tolerance-primer-dimer-folding-energy', type=float, default=GoldenGatePool.DEFAULT_TOLERANCE_DG_PRIMER_DIMER)
    parser.add_argument('--primer-sequence-constraint', default=None)
    parser.add_argument('--primer-maximum-repeat-length', type=int, default=GoldenGatePool.DEFAULT_MAXIMUM_REPEAT_LENGTH)

    parser.add_argument('--parallel', action='store_true')
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--screening', action='store_true', help='Enable external DNA screening requests.')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--config-json', default=None, help='JSON object merged into oligoPoolSpecification before CLI values.')
    return parser

def parse_excel_sheet_name(value):
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.isdigit():
            return int(stripped)
        return stripped
    return value

def build_oligo_pool_specification(args):
    if args.maximum_oligo_length <= 0:
        raise ValueError("--maximum-oligo-length must be positive.")
    if args.primer_sequence_length <= 0:
        raise ValueError("--primer-sequence-length must be positive.")
    if args.minimum_fragment_length <= 0:
        raise ValueError("--minimum-fragment-length must be positive.")
    if args.long_assembly_threshold <= 0:
        raise ValueError("--long-assembly-threshold must be positive.")
    if args.long_assembly_maximum_fragment_length <= 0:
        raise ValueError("--long-assembly-maximum-fragment-length must be positive.")
    if args.long_assembly_minimum_fragment_length <= 0:
        raise ValueError("--long-assembly-minimum-fragment-length must be positive.")
    if args.num_fragments_per_long_assembly <= 1:
        raise ValueError("--num-fragments-per-long-assembly must be greater than 1.")

    config = {}
    if args.config_json:
        try:
            config = json.loads(args.config_json)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid --config-json: {exc}") from exc
        if not isinstance(config, dict):
            raise ValueError("--config-json must decode to a JSON object.")

    assembly_enzyme = normalize_build_enzyme(args.assembly_enzyme)
    if assembly_enzyme is None:
        raise ValueError("--assembly-enzyme cannot be none.")
    level_two_assembly_enzyme = normalize_build_enzyme(args.level_two_assembly_enzyme)

    primer_sequence_constraint = args.primer_sequence_constraint
    if primer_sequence_constraint is None:
        primer_sequence_constraint = 'N' * max(args.primer_sequence_length - 2, 0) + 'WW'
    if len(primer_sequence_constraint) != args.primer_sequence_length:
        raise ValueError("--primer-sequence-constraint length must equal --primer-sequence-length.")

    explicit_config = {
        'PCR_Tm_target' : args.pcr_tm_target,
        'excluded_restriction_enzymes' : parse_comma_list(args.excluded_restriction_enzymes),
        'maximum_oligo_length' : args.maximum_oligo_length,
        'primer_sequence_length' : args.primer_sequence_length,
        'padding_end' : args.padding_end,
        'assembly_enzyme' : assembly_enzyme,
        'level_two_assembly_enzyme' : level_two_assembly_enzyme,
        'MINIMUM_FRAGMENT_LENGTH' : args.minimum_fragment_length,
        'LONG_ASSEMBLY_THRESHOLD' : args.long_assembly_threshold,
        'LONG_ASSEMBLY_MAXIMUM_FRAGMENT_LENGTH' : args.long_assembly_maximum_fragment_length,
        'LONG_ASSEMBLY_MINIMUM_FRAGMENT_LENGTH' : args.long_assembly_minimum_fragment_length,
        'NUM_FRAGMENTS_PER_LONG_ASSEMBLY' : args.num_fragments_per_long_assembly,
        'tolerance_delta_Tm' : args.tolerance_delta_tm,
        'tolerance_dG_folding' : args.tolerance_dg_folding,
        'tolerance_dG_homodimer' : args.tolerance_dg_homodimer,
        'tolerance_primer_dimer_folding_energy' : args.tolerance_primer_dimer_folding_energy,
        'primer_sequence_constraint' : primer_sequence_constraint,
        'primer_maximum_repeat_length' : args.primer_maximum_repeat_length,
    }
    config.update(explicit_config)
    return config

def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    try:
        if args.export_build_json is not None:
            summary = exportBuildJsonExcelFiles(args.export_build_json, args.output_prefix, verbose=args.verbose)
            print(json.dumps({key : value for (key, value) in summary.items() if key != 'export_data'}, indent=2, sort_keys=True))
            return 0

        if args.workers is not None:
            if args.workers <= 0:
                raise ValueError("--workers must be positive.")
            os.environ["NUM_WORKERS"] = str(args.workers)

        oligoPoolSpecification = build_oligo_pool_specification(args)
        input_filename = args.json_input if args.json_input is not None else args.input
        if input_filename is None:
            raise ValueError("An input file is required. Provide positional input or --json-input.")
        file_format = 'json' if args.json_input is not None and args.file_format == 'auto' else args.file_format
        sequenceSpecList = loadSequencesFromFile(
            input_filename,
            file_format=file_format,
            sheet_name=parse_excel_sheet_name(args.sheet_name),
            default_assembly_enzyme=oligoPoolSpecification['assembly_enzyme'],
            default_level_two_assembly_enzyme=oligoPoolSpecification['level_two_assembly_enzyme'],
            default_is_circular=args.default_is_circular,
        )
        if len(sequenceSpecList) == 0:
            raise ValueError("No sequences were loaded from the input file.")

        print("# Sequences: {}".format(len(sequenceSpecList)))
        t1 = time.time()

        build = GoldenGatePool(verbose=args.verbose, parallel=args.parallel, doScreening=args.screening)
        oligoPoolSpecification = build.addMultipleAssemblies(oligoPoolSpecification, sequenceSpecList)
        if len(oligoPoolSpecification.get('library_specifications', [])) == 0:
            for (assembly, msg) in oligoPoolSpecification.get('failed_assemblies', []):
                print(msg)
            raise ValueError("No successful assemblies were generated; export skipped.")

        oligoPoolSpecification = build.finalize(oligoPoolSpecification)

        if 'failed_assemblies' in oligoPoolSpecification:
            for (assembly, msg) in oligoPoolSpecification['failed_assemblies']:
                print(msg)

        build.output(oligoPoolSpecification)
        build.export(
            oligoPoolSpecification,
            args.output_prefix + '.xlsx',
            args.output_prefix + '_oligos.xlsx',
            args.output_prefix + '_primer_plates.xlsx'
        )

        t2 = time.time()
        print("GoldenGatePool Elapsed Time: ", t2 - t1, " seconds.")
        return 0
    except Exception as exc:
        parser.error(str(exc))

if __name__ == "__main__":
    sys.exit(main())
