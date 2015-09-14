from __future__ import print_function

from test import test_circuit
from circuit import ParseException

import argparse, os, sys, time
import utils

__all__ = ['main']

errorstr = utils.clr_error('Error:')

def test_all(args, bpclass, obfclass, obfuscate):
    if not os.path.isdir(args.test_all):
        print("%s '%s' is not a valid directory" % (errorstr, args.test_all))
        sys.exit(1)
    ext = '.acirc' if args.zimmerman else '.circ'
    for circuit in os.listdir(args.test_all):
        path = os.path.join(args.test_all, circuit)
        if os.path.isfile(path) and path.endswith(ext):
            test_circuit(path, bpclass, obfclass, obfuscate, args)

def check_args(args):
    num_set = int(args.ananth_etal) + int(args.sahai_zhandry) + int(args.zimmerman)
    if num_set > 1:
        print('%s Only one of --ananth-etal, --sahai-zhandry, --zimmerman can be set' % errorstr)
        sys.exit(1)
    if num_set == 0:
        args.ananth_etal = True

def bp(args):
    check_args(args)

    if args.ananth_etal:
        from agis_bp import AGISBranchingProgram
        cls = AGISBranchingProgram
    if args.sahai_zhandry:
        from sz_bp import SZBranchingProgram
        cls = SZBranchingProgram
    if args.zimmerman:
        from z_obfuscator import Circuit
        cls = Circuit

    try:
        if args.test_circuit:
            test_circuit(args.test_circuit, cls, None, False, args)
        elif args.test_all:
            test_all(args, cls, None, False)
        elif args.load_circuit:
            bp = cls(args.load_circuit, verbose=args.verbose)
            if args.obliviate:
                bp.obliviate()
            if args.eval:
                r = bp.evaluate(args.eval)
                print('Output = %d' % r)
    except ParseException as e:
        print('%s %s' % (errorstr, e))
        sys.exit(1)

def obf(args):
    check_args(args)

    if args.ananth_etal:
        from agis_bp import AGISBranchingProgram
        from agis_obfuscator import AGISObfuscator
        bpclass = AGISBranchingProgram
        obfclass = AGISObfuscator
    if args.sahai_zhandry:
        from sz_bp import SZBranchingProgram
        from sz_obfuscator import SZObfuscator
        bpclass = SZBranchingProgram
        obfclass = SZObfuscator
    if args.zimmerman:
        from z_obfuscator import ZObfuscator
        bpclass = None
        obfclass = ZObfuscator

    try:
        if args.test_circuit:
            test_circuit(args.test_circuit, bpclass, obfclass, True, args)
        elif args.test_all:
            test_all(args, bpclass, obfclass, True)
        else:
            directory = None
            if args.load_obf:
                directory = args.load_obf
            elif args.load_circuit:
                start = time.time()
                obf = obfclass(verbose=args.verbose, nthreads=args.nthreads)
                directory = args.save if args.save \
                            else '%s.obf.%d' % (args.load_circuit, args.secparam)
                obf.obfuscate(args.load_circuit, args.secparam, directory,
                              obliviate=args.obliviate, nslots=args.nslots,
                              kappa=args.kappa)
                end = time.time()
                print('Obfuscation took: %f seconds' % (end - start))
            else:
                print('%s One of --load-obf, --load-circuit, or '
                      '--test-circuit must be used' % errorstr)
                sys.exit(1)

            if args.eval:
                assert directory
                obf = obfclass(verbose=args.verbose, nthreads=args.nthreads)
                r = obf.evaluate(directory, args.eval)
                print('Output = %d' % r)
    except ParseException as e:
        print('%s %s' % (errorstr, e))
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Cryptographic program obfuscator.')
    subparsers = parser.add_subparsers()

    try:
        nthreads = os.sysconf('SC_NPROCESSORS_ONLN')
    except ValueError:
        print(utils.clr_warn('Warning: Unable to count number of cores, defaulting to 1'))
        nthreads = 1
    secparam = 24

    parser_bp = subparsers.add_parser(
        'bp',
        help='commands for circuit -> branching program conversion')
    parser_bp.add_argument('--eval',
                           metavar='INPUT', action='store', type=str,
                           help='evaluate branching program on INPUT')
    parser_bp.add_argument('--load-circuit',
                           metavar='FILE', action='store', type=str,
                           help='load circuit from FILE')
    parser_bp.add_argument('--test-circuit',
                           metavar='FILE', action='store', type=str,
                           help='test BP conversion for FILE')
    parser_bp.add_argument('--test-all',
                           metavar='DIR', nargs='?', const='circuits/',
                           help='test BP conversion for all circuits in DIR (default: %(const)s)')
    parser_bp.add_argument('--secparam',
                           metavar='N', action='store', type=int, default=secparam,
                           help='security parameter (default: %(default)s)')
    parser_bp.add_argument('--obliviate',
                           action='store_true',
                           help='obliviate the branching program')
    parser_bp.add_argument('-v', '--verbose',
                           action='store_true',
                           help='be verbose')
    parser_bp.add_argument('-a', '--ananth-etal',
                            action='store_true',
                            help='use the Ananth et al. construction (default)')
    parser_bp.add_argument('-s', '--sahai-zhandry',
                           action='store_true',
                           help='use the Sahai/Zhandry construction')
    parser_bp.add_argument('-z', '--zimmerman',
                           action='store_true',
                           help='use the Zimmerman construction')
    parser_bp.set_defaults(func=bp)

    parser_obf = subparsers.add_parser(
        'obf',
        help='commands for obfuscating a circuit/branching program')
    parser_obf.add_argument('--eval',
                            metavar='INPUT', action='store', type=str,
                            help='evaluate obfuscation on INPUT')
    parser_obf.add_argument('--kappa',
                            metavar='N', action='store', type=int, default=None,
                            help='set kappa to N (for debugging)')
    parser_obf.add_argument('--load-obf',
                            metavar='DIR', action='store', type=str,
                            help='load obfuscation from DIR')
    parser_obf.add_argument('--load-circuit',
                            metavar='FILE', action='store', type=str,
                            help='load circuit from FILE and obfuscate')
    parser_obf.add_argument('--test-circuit',
                            metavar='FILE', action='store', type=str,
                            help='test circuit from FILE')
    parser_obf.add_argument('--test-all',
                            metavar='DIR', nargs='?', const='circuits/',
                            help='test obfuscation for all circuits in DIR (default: %(const)s)')
    parser_obf.add_argument('--save',
                            metavar='DIR', action='store', type=str,
                            help='save obfuscation to DIR')
    parser_obf.add_argument('--secparam',
                            metavar='N', action='store', type=int,
                            default=secparam, help='security parameter (default: %(default)s)')
    parser_obf.add_argument('--obliviate',
                            action='store_true',
                            help='obliviate the branching program')
    parser_obf.add_argument('--nthreads',
                            metavar='N', action='store', type=int, default=nthreads,
                            help='number of threads to use (default: %(default)s)')
    parser_obf.add_argument('--nslots',
                            metavar='N', action='store', type=int, default=secparam,
                            help='number of slots to fill (default: %(default)s)')
    parser_obf.add_argument('-v', '--verbose',
                            action='store_true', 
                            help='be verbose')
    parser_obf.add_argument('-a', '--ananth-etal',
                            action='store_true',
                            help='use the Ananth et al. construction (default)')
    parser_obf.add_argument('-s', '--sahai-zhandry',
                            action='store_true',
                            help='use the Sahai/Zhandry construction')
    parser_obf.add_argument('-z', '--zimmerman',
                            action='store_true',
                            help='use the Zimmerman construction')
    parser_obf.set_defaults(func=obf)

    args = parser.parse_args()
    args.func(args)
