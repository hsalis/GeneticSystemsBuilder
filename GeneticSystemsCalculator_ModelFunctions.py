import os, copy, sys, random, re, traceback
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def createGlobalModelFunction(inputs):
    inputs = copy.copy(inputs)
    includedModelFunctions = []
    for (rule, active) in list(inputs['activeRules'].items()):
        if active and rule in modelFunctionRules: includedModelFunctions.extend( modelFunctionRules[rule] )

    includedModelFunctions.sort(key = lambda x: x['rank'])

    inputs['upstream_sequence'] = inputs.get('upstream_sequence','').upper()
    inputs['downstream_sequence'] = inputs.get('downstream_sequence','').upper()
    inputs['model_function_type'] = 'global'

    global modelFailuresByPosition
    modelFailuresByPosition = defaultdict(int, { nt : 0 for nt in range(0,inputs['sequence_length']+1) })

    def retFcn(seq):
        seqlen = len(seq)
        for fcnInfo in includedModelFunctions:
            try:
                (Passed, traceback_pos) = fcnInfo['function'](seq, inputs)
                #if not isinstance(Passed, bool): Passed = Passed[0]  # For when a local model function is used as a global model function.

                if not Passed:
                    modelFailuresByPosition[traceback_pos] += 1
                    modelFailuresByPosition[seqlen] += 1
                    if inputs['verbosity'] > 1: print('** GLOBAL RULE {} FAILED at position {} **'.format(fcnInfo['name'], traceback_pos))
                    return False
            except:
                print('** ' + fcnInfo['name'] + ' **')
                print(traceback.format_exc())
                raise

        return True
    return lambda seq: retFcn(seq)

def createLocalModelFunction(inputs):
    inputs = copy.copy(inputs)

    includedModelFunctions = []
    for (rule, active) in list(inputs['activeRules'].items()):
        if active and rule in modelFunctionRules: includedModelFunctions.extend( modelFunctionRules[rule] )

    includedModelFunctions.sort(key = lambda x: x['rank'])

    inputs['upstream_sequence'] = inputs.get('upstream_sequence','').upper()
    inputs['downstream_sequence'] = inputs.get('downstream_sequence','').upper()
    inputs['model_function_type'] = 'local'

    global modelFailuresByPosition
    modelFailuresByPosition = defaultdict(int, { nt : 0 for nt in range(0, inputs['sequence_length']+1) })

    global modelFunctionCounter
    modelFunctionCounter = 0

    def retFcn(seq):
        global modelFunctionCounter
        seqlen = len(seq)
        for fcnInfo in includedModelFunctions:
            if inputs['verbosity'] >= 5: print('Calling Model Function {} on sequence {}'.format(fcnInfo['function'], seq))
            try:
                (Passed, pos) = fcnInfo['function'](seq, inputs)
            except:
                print('** ' + fcnInfo['name'] + ' **')
                print(traceback.format_exc())
                raise

            if not Passed:
                if inputs['verbosity'] > 1: print('** ' + fcnInfo['name'] + ' FAILED at POSITION %s **' % pos)
                modelFailuresByPosition[seqlen] += 1
                if inputs['verbosity'] > 1 and modelFailuresByPosition[seqlen] % 10 == 0: print("modelFailuresByPosition[seqlen]: ", modelFailuresByPosition[seqlen], "at seqlen: ", seqlen)
                if inputs['verbosity'] > 2: print('Top 20 positions: ', sorted( [(x, y) for (x, y) in list(modelFailuresByPosition.items())], key = lambda x: x[1], reverse = True)[0:20])
                return (False, pos)

        modelFunctionCounter += 1
        if inputs['verbosity'] > 1 and modelFunctionCounter % 100 == 0: print("Model Function Counter: %s" % modelFunctionCounter)

        if inputs['verbosity'] > 1 and seqlen % 100 == 0: print("Sequence Position: ", seqlen)
        return (True, None)

    return lambda seq: retFcn(seq)

def createGlobalModelFunctionWithInformation(inputs):
    inputs = copy.copy(inputs)
    includedModelFunctions = []
    for (rule, active) in list(inputs['activeRules'].items()):
        if active and rule in modelFunctionRules: includedModelFunctions.extend( modelFunctionRules[rule] )

    includedModelFunctions.sort(key = lambda x: x['rank'])

    inputs['upstream_sequence'] = inputs.get('upstream_sequence','').upper()
    inputs['downstream_sequence'] = inputs.get('downstream_sequence','').upper()
    inputs['model_function_type'] = 'global'

    global modelFailuresByPosition
    modelFailuresByPosition = defaultdict(int, { nt : 0 for nt in range(0,inputs['sequence_length']+1) })

    def retFcn(seq):
        seqlen = len(seq)
        for fcnInfo in includedModelFunctions:
            try:
                (Passed, traceback_pos) = fcnInfo['function'](seq, inputs)

                if not Passed:
                    modelFailuresByPosition[traceback_pos] += 1
                    modelFailuresByPosition[seqlen] += 1
                    if inputs['verbosity'] > 1: print('** GLOBAL RULE {} FAILED at position {} **'.format(fcnInfo['id'], traceback_pos))
                    return (False, traceback_pos, fcnInfo['id'])
            except:
                print('** ' + fcnInfo['name'] + ' **')
                print(traceback.format_exc())
                raise

        return (True, None, None)
    return lambda seq: retFcn(seq)

def MeltingTemperature(seq, inputs):
    Tm_settings = inputs['ruleInputs']['MeltingTemperature']['TM_SETTINGS']
    TM_bounds = inputs['ruleInputs']['MeltingTemperature']['TM_BOUNDS']  #[40.0, 75.0]

    Tm = mt.Tm_NN(seq, **Tm_settings)
    if Tm < TM_bounds[0]:
        return (False, 0)
    if Tm > TM_bounds[1]:
        return (False, 0)
    return (True, None)

modelFunctionInfo = {}
modelFunctionInfo['MeltingTemperature'] = {'id' : 'MeltingTemperature', 'name' : 'Target Melting Temperature', 'function' : MeltingTemperature, 'rank' : 2 }

modelFunctionRules = { 'rule_Tm_bounds' : [ modelFunctionInfo['MeltingTemperature'] ],
                      }