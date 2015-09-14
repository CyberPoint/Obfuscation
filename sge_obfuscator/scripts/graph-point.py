#!/usr/bin/env python2

from __future__ import division
import os, sys

from operator import add, sub
import matplotlib.pyplot as plt
import numpy as np
import utils

def obftime(ax1):
    xs = (8, 12, 16)

    path8 = os.path.join('results', 'point.52', 'point-8.circ-52-obf-time.log')
    path12 = os.path.join('results', 'point.52', 'point-12.circ-52-obf-time.log')
    path16 = os.path.join('results', 'point.52', 'point-16.circ-52-obf-time.log')

    bar8 = utils.extract_obf_time(path8)
    bar12 = utils.extract_obf_time(path12)
    bar16 = utils.extract_obf_time(path16)

    all = zip(bar8, bar12, bar16)
    all = [(a / 60 / 60, b / 60 / 60, c / 60 / 60) for a, b, c in all]

    ind = np.arange(len(all[0]))
    width = 0.4

    encodingtime = map(add, all[6], all[7])

    total = ax1.bar(ind + width, all[8], width, color='white')
    mlm = ax1.bar(ind + width, all[4], width, color='black')
    enc = ax1.bar(ind + width, encodingtime, width, color='gray', bottom=all[4])

    ax1.set_ylabel(r'Obfuscation time (hr)')
    ax1.set_xlabel(r'Input size of point function')
    ax1.set_xticks(ind + 0.6)
    ax1.set_ylim(0, 10)
    ax1.set_xticklabels(xs)

    ax1.legend((mlm[0], enc[0]),
               ('Param gen', 'Encoding'),
               loc='upper left')

def obfsize(ax1):
    xs = (8, 12, 16)

    path8 = os.path.join('results', 'point.52', 'point-8.circ-52-obf-size.log')
    path12 = os.path.join('results', 'point.52', 'point-12.circ-52-obf-size.log')
    path16 = os.path.join('results', 'point.52', 'point-16.circ-52-obf-size.log')

    _, bar8 = utils.obfsize(path8)
    _, bar12 = utils.obfsize(path12)
    _, bar16 = utils.obfsize(path16)

    all = [bar8, bar12, bar16]
    all = [x / 2 ** 30 for x in all]

    print('obf size: %s' % all)

    ind = np.arange(len(all))
    width = 0.4

    total = ax1.bar(ind + width, all, width, color='gray')
    ax1.set_ylabel(r'Obfuscation size (GB)')
    ax1.set_xlabel(r'Input size of point function')
    ax1.set_xticks(ind + 0.6)
    ax1.set_xticklabels(xs)

def evaltime(ax1):
    xs = (8, 12, 16)

    path8 = os.path.join('results', 'point.52', 'point-8.circ-52-eval-time.log')
    path12 = os.path.join('results', 'point.52', 'point-12.circ-52-eval-time.log')
    path16 = os.path.join('results', 'point.52', 'point-16.circ-52-eval-time.log')

    _, bar8 = utils.evaltime(path8)
    _, bar12 = utils.evaltime(path12)
    _, bar16 = utils.evaltime(path16)

    all = [bar8, bar12, bar16]
    all = [x / 60 / 60 for x in all]

    print('eval time: %s' % all)

    ind = np.arange(len(all))
    width = 0.4

    total = ax1.bar(ind + width, all, width, color='gray')
    ax1.set_ylabel(r'Evaluation time (hr)')
    ax1.set_xlabel(r'Input size of point function')
    ax1.set_xticks(ind + 0.6)
    ax1.set_xticklabels(xs)

def main(argv):
    utils.init_wide()

    fig, axes = plt.subplots(nrows=1, ncols=3)
    obftime(axes.flat[0])
    obfsize(axes.flat[1])
    evaltime(axes.flat[2])

    plt.subplots_adjust(wspace=0.4, bottom=0.2)
    plt.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
