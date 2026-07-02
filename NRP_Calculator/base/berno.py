from math          import log, ceil
from scipy.special import betaincinv as ribf

import random

zero_const = -24 # Lower is Rarer # Originally: 1. / 10**10

def get_trials(prob):
    return 2 * int(ceil(-24 / log(1. - prob)))

def get_prob(trials, success):
    # probability for rare event success
    return (success + 1.) / (trials + 2.)

def get_CI(trials, success, conf_interval):
    x1 = ribf(success + 1, trials - success + 1, 0.5*(1 - conf_interval))
    x2 = ribf(success + 1, trials - success + 1, 0.5*(1 + conf_interval))
    return (x1, x2)

def test_func():
    if random.random() <= 0.002: #1. / 200:
        return True
    return False