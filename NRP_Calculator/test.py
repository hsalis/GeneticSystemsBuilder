import random, copy, re, traceback

def object2Dictionary(outputData):
    #Everything is an object! This will convert all user-defined classes into dictionaries.
    #Using duck typing.
    import numpy, collections
    from datetime import datetime

    def isinstance_namedtuple(x):
      return (isinstance(x, tuple) and
          isinstance(getattr(x, '__dict__', None), collections.Mapping) and
          getattr(x, '_fields', None) is not None)

    try:

        if isinstance(outputData,collections.OrderedDict):
            return outputData

        for (key, value) in list(outputData.items()):
            try:
                if isinstance(value,collections.OrderedDict):
                    outputData[key] = object2Dictionary(value)
                else:
                    outputData[key] = object2Dictionary(value.__dict__)
            except:
                outputData[key] = object2Dictionary(value)
        return outputData

    except:

        try:
            if isinstance(outputData, list):
                for (i,item) in enumerate(outputData):
                    try:
                        if isinstance(item,collections.OrderedDict):
                            outputData[i] = object2Dictionary(item)
                        else:
                            outputData[i] = object2Dictionary(item.__dict__)
                    except:
                        outputData[i] = object2Dictionary(item)
                return outputData

            if isinstance(outputData, tuple):
                temp = []
                for (i, item) in enumerate(outputData):
                    try:
                        if isinstance(item,collections.OrderedDict):
                            temp.append( object2Dictionary(item) )
                        else:
                            temp.append( object2Dictionary(item.__dict__) )
                    except:
                        temp.append( object2Dictionary(item) )
                outputData = tuple(temp)
                return outputData

            if isinstance_namedtuple(outputData):
                return outputData

            if outputData is None or isinstance(outputData, (int, float, str, bool, complex) ):
                return outputData

            try:
                return outputData.isoformat()
            except:
                pass

            try:
                #If it's a numpy array
                return outputData.tolist()
            except:
                pass

            return object2Dictionary(outputData.__dict__)

        except:
            print("Could not convert to dictionary!", outputData)
            print("Error Traceback:", traceback.format_exc())
            return outputData

def NonRepetitivePartsCalculator_FinderMode(job):

    import nrpcalc

    inputData = copy.copy(job['inputData'])
    Lmax = inputData['max_repeat_size']
    outputData = {}
    outputData['non_repetitive_indices'] = {}

    parts_sequence_list = [str(x).upper() for x in inputData['part_sequence_list']][:]

    while True:
        if len(inputData['background_sequence']) > 0:
            background_sequence = str(inputData['background_sequence']).upper()
            bkg = nrpcalc.background(path='./nrpcalc_bkg/', Lmax = Lmax)
            bkg.add(background_sequence)
        else:
            bkg = None

        try:
            output = nrpcalc.finder(parts_sequence_list,
                            Lmax,
                            internal_repeats=False,
                            background = bkg,
                            verbose=True)
        except:
            bkg.drop()
            raise

        outputData['non_repetitive_indices'][Lmax] = sorted(output.keys())
        bkg.drop()

        if Lmax <= 6 or len(outputData['non_repetitive_indices'][Lmax]) <= 1:
            break
        else:
            Lmax = Lmax - 1

    job['outputData'] = object2Dictionary(outputData)
    return job

def NonRepetitivePartsCalculator_MakerMode(job):

    import nrpcalc
    from pre_process_maker import is_seq_constr_sufficient

    inputData = copy.copy(job['inputData'])
    outputData = {}

    sequenceConstraint = str(inputData['sequence_constraint'].replace('U','T'))
    Lmax = inputData['max_repeat_size']
    allowInternalRepeats = inputData['allow_internal_repeats']

    if inputData['RNA_structural_constraint'] is None or len(inputData['RNA_structural_constraint']) == 0:
        RNA_structural_constraint = "." * len(sequenceConstraint)
        part_type = 'DNA'
    else:
        RNA_structural_constraint = str(inputData['RNA_structural_constraint'])
        part_type = 'RNA'

    if inputData['checkbox_goal_maximize']:
        numParts = 1000000
    else:
        numParts = min(inputData['target_number'], 1000000)

    if len(inputData['background_sequence']) > 0:
        background_sequence = str(inputData['background_sequence']).upper()
        bkg = nrpcalc.background(path='./nrpcalc_bkg/', Lmax = Lmax)
        bkg.add(background_sequence)
    else:
        bkg = None

    #Check if sequence constraint is degenerate enough for the desired homology and toolbox size
    constraintCheckResult = is_seq_constr_sufficient(sequenceConstraint, Lmax-1, toolbox_size=numParts)
    if constraintCheckResult[0]: #everything is fine
        outputData['constraintCheckResult'] = (True, None)
    else: #constraint not sufficiently degenerate
        badConstraintSeqList = []
        for badConstraintSeq in constraintCheckResult[1]:
            badConstraintSeqList.append( sequenceConstraint[badConstraintSeq[0]:badConstraintSeq[1]] )
        outputData['constraintCheckResult'] = (False, badConstraintSeqList)

    from model_functions import modelFunction_GC_Content, modelFunction_Exclude_Sequences

    def multi_model_local_function(seq, modelFunctionInputList):
        for (modelFunction, inputs) in modelFunctionInputList:
            (Passed, index) = modelFunction(seq, inputs)
            if not Passed: return (False, index)
        return (True, None)

    def multi_model_global_function(seq, modelFunctionInputList):
        for (modelFunction, inputs) in modelFunctionInputList:
            Passed = modelFunction(seq, inputs)
            if not Passed: return False
        return True

    GC_percentage_constraint = inputData['GC_percentage']

    if len(inputData['restriction_enzymes_to_avoid']) > 0:
        compiledRegExp = re.compile("|".join([x['sequence'] for x in inputData['restriction_enzymes_to_avoid']]))
        local_model_function = lambda seq: modelFunction_Exclude_Sequences(seq, compiledRegExp)
    else:
        local_model_function = None

    if (GC_percentage_constraint[0] > 0 or GC_percentage_constraint[1] < 100):
        global_model_function = lambda seq: modelFunction_GC_Content(seq, GC_percentage_constraint[0], GC_percentage_constraint[1])
    else:
        global_model_function = None

    try:
        parts_dict = nrpcalc.maker(seq_constr = sequenceConstraint,
                            struct_constr = RNA_structural_constraint,
                            part_type = part_type,
                            Lmax = Lmax,
                            target_size = numParts,
                            internal_repeats=allowInternalRepeats,
                            background=bkg,
                            struct_type='mfe',
                            seed=None,
                            synth_opt=inputData['reduce_synthesis_complexity'],
                            local_model_fn=local_model_function,
                            global_model_fn=global_model_function,
                            jump_count=1000,
                            fail_count=10000,
                            output_file=None,
                            verbose=True)
    except:
        bkg.drop()
        raise

    NonRepetitiveToolbox = list(parts_dict.values())
    outputData['non_repetitive_sequences'] = NonRepetitiveToolbox
    outputData['non_repetitive_indices'] = {}
    bkg.drop()

    if len(NonRepetitiveToolbox) > 0:
        parts_sequence_list = [str(x) for x in NonRepetitiveToolbox][:]
        while True:
            if len(inputData['background_sequence']) > 0:
                background_sequence = str(inputData['background_sequence']).upper()
                bkg = nrpcalc.background(path='./nrpcalc_bkg/', Lmax = Lmax)
                bkg.add(background_sequence)
            else:
                bkg = None

            output = nrpcalc.finder(parts_sequence_list,
                                Lmax,
                                internal_repeats=False,
                                background = bkg,
                                verbose=True)
            outputData['non_repetitive_indices'][Lmax] = sorted(output.keys())
            bkg.drop()

            if Lmax <= 6 or len(outputData['non_repetitive_indices'][Lmax]) <= 1:
                break
            else:
                Lmax = Lmax - 1

    job['outputData'] = object2Dictionary(outputData)
    return job

def NonRepetitivePartsCalculator_CDSMakerMode(job):

    import nrpcalc
    from pre_process_maker import is_seq_constr_sufficient

    inputData = copy.copy(job['inputData'])
    outputData = {}

    AA_seq = inputData['aa_sequence_constraint'].upper()
    sequenceConstraint = "N" * (3 * len(AA_seq) )
    Lmax = inputData['max_repeat_size']
    allowInternalRepeats = inputData['allow_internal_repeats']
    RNA_structural_constraint = "." * len(sequenceConstraint)
    part_type = 'DNA'

    if inputData['checkbox_goal_maximize']:
        numParts = 1000000
    else:
        numParts = min(inputData['target_number'], 1000000)

    if len(inputData['background_sequence']) > 0:
        background_sequence = str(inputData['background_sequence']).upper()
        bkg = nrpcalc.background(path='./nrpcalc_bkg/', Lmax = Lmax)
        bkg.add(background_sequence)
    else:
        bkg = None

    from model_functions import modelFunction_GC_Content, modelFunction_Exclude_Sequences, modelFunction_AASeq

    def multi_model_local_function(seq, modelFunctionInputList):
        for (modelFunction, inputs) in modelFunctionInputList:
            (Passed, index) = modelFunction(seq, inputs)
            if not Passed: return (False, index)
        return (True, None)

    def multi_model_global_function(seq, modelFunctionInputList):
        for (modelFunction, inputs) in modelFunctionInputList:
            Passed = modelFunction(seq, inputs)
            if not Passed: return False
        return True

    GC_percentage_constraint = inputData['GC_percentage']
    if (GC_percentage_constraint[0] > 0 or GC_percentage_constraint[1] < 100):
        global_model_function = lambda seq: modelFunction_GC_Content(seq, GC_percentage_constraint[0], GC_percentage_constraint[1])
    else:
        global_model_function = None

    if len(inputData['restriction_enzymes_to_avoid']) > 0:
        compiledRegExp = re.compile("|".join([x['sequence'] for x in inputData['restriction_enzymes_to_avoid']]))
        modelFunctionInputList = [(modelFunction_Exclude_Sequences, compiledRegExp), (modelFunction_AASeq, AA_seq)]
        local_model_function = lambda seq: multi_model_local_function(seq, modelFunctionInputList)
    else:
        local_model_function = lambda seq: modelFunction_AASeq(seq, AA_seq)

    try:
        parts_dict = nrpcalc.maker(seq_constr = sequenceConstraint,
                            struct_constr = RNA_structural_constraint,
                            part_type = part_type,
                            Lmax = Lmax,
                            target_size = numParts,
                            internal_repeats=allowInternalRepeats,
                            background=bkg,
                            struct_type='mfe',
                            seed=None,
                            synth_opt=inputData['reduce_synthesis_complexity'],
                            local_model_fn=local_model_function,
                            global_model_fn=global_model_function,
                            jump_count=1000,
                            fail_count=10000,
                            output_file=None,
                            verbose=True)
    except:
        bkg.drop()
        raise

    NonRepetitiveToolbox = list(parts_dict.values())
    outputData['non_repetitive_sequences'] = NonRepetitiveToolbox
    outputData['non_repetitive_indices'] = {}
    bkg.drop()

    if len(NonRepetitiveToolbox) > 0:
        parts_sequence_list = [str(x) for x in NonRepetitiveToolbox][:]
        while True:
            if len(inputData['background_sequence']) > 0:
                background_sequence = str(inputData['background_sequence']).upper()
                bkg = nrpcalc.background(path='./nrpcalc_bkg/', Lmax = Lmax)
                bkg.add(background_sequence)
            else:
                bkg = None

            output = nrpcalc.finder(parts_sequence_list,
                                Lmax,
                                internal_repeats=False,
                                background = bkg,
                                verbose=True)
            outputData['non_repetitive_indices'][Lmax] = sorted(output.keys())
            bkg.drop()

            if Lmax <= 6 or len(outputData['non_repetitive_indices'][Lmax]) <= 1:
                break
            else:
                Lmax = Lmax - 1

    job['outputData'] = object2Dictionary(outputData)
    return job

def testFinder():

    job = {}
    job['inputData'] = {}
    job['inputData']['part_sequence_list'] = [ "".join([random.choice(['A','C','G','T']) for n in range(200)]) for m in range(1000)]
    job['inputData']['max_repeat_size'] = 15
    job['inputData']['background_sequence'] = "".join([random.choice(['A','C','G','T']) for n in range(10000)])
    job = NonRepetitivePartsCalculator_FinderMode(job)

    print(job['outputData'])

def testMaker():

    job = {}
    job['inputData'] = {}
    job['inputData']['sequence_constraint'] = 'S'*5 + 'N'*15 + 'TTGACA' + 'N'*6 + 'CCN' + 'N'*8 + 'TATAAT' + 'N'*6
    job['inputData']['RNA_structural_constraint'] = None
    job['inputData']['checkbox_goal_maximize'] = False
    job['inputData']['target_number'] = 10000
    job['inputData']['max_repeat_size'] = 15
    job['inputData']['allow_internal_repeats'] = False
    job['inputData']['reduce_synthesis_complexity'] = True
    job['inputData']['GC_percentage'] = [0.0, 100.0]
    job['inputData']['background_sequence'] = "".join([random.choice(['A','C','G','T']) for n in range(10000)])
    job['inputData']['restriction_enzymes_to_avoid'] = [ {'sequence' : 'TCTAGA', 'enzyme' : 'XbaI'} ]

    job = NonRepetitivePartsCalculator_MakerMode(job)
    print(job['outputData'])

def testCDSMaker():

    from model_functions import AA2Codon

    job = {}
    job['inputData'] = {}
    job['inputData']['aa_sequence_constraint'] = "M" + "".join([random.choice(list(AA2Codon.keys())) for n in range(299)])
    job['inputData']['checkbox_goal_maximize'] = False
    job['inputData']['target_number'] = 10
    job['inputData']['max_repeat_size'] = 20
    job['inputData']['allow_internal_repeats'] = False
    job['inputData']['reduce_synthesis_complexity'] = False
    job['inputData']['GC_percentage'] = [-1, 1000.0]
    job['inputData']['background_sequence'] = "".join([random.choice(['A','C','G','T']) for n in range(10000)])
    job['inputData']['restriction_enzymes_to_avoid'] = [ {'sequence' : 'TCTAGA', 'enzyme' : 'XbaI'} ]

    job = NonRepetitivePartsCalculator_CDSMakerMode(job)
    print(job['outputData'])

if __name__ == "__main__":
    testFinder()
    testMaker()
    testCDSMaker()


