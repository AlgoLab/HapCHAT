usage = '''

usage: python increments.py [file.log]

given a log file file.log, the output of hapchat onto stdout, obtain
the sequence of increments at the sites where the k was incremented

'''

import sys
import argparse


# print dictionary
def print_dic(steps) :

    print(*'step sequence'.split(), sep='\t')
    for step in sorted(steps) :
        print(step, end='\t')
        print(*sorted(steps[step]))


# print reverse dictionary
def print_rev_dic(steps) :
    rev_dic = {}

    for step in steps :
        key = ' '.join(str(x) for x in sorted(steps[step]))

        if key not in rev_dic :
            rev_dic[key] = []

        rev_dic[key].append(step)

    for key in sorted(rev_dic) :
        print(key, end='\t')
        print(*sorted(rev_dic[key]), end='\t')
    print()

#
# Parser
#----------------------------------------------------------------------

parser = argparse.ArgumentParser(description='''

   given a log file, the output of hapchat onto stdout, obtain the
   sequence of increments at the sites where the k was incremented

''')

parser.add_argument(
    'file',
    type = argparse.FileType('r'), default = sys.stdin,
    help = 'input log file')
parser.add_argument(
    '-r', '--reverse',
    dest = 'reverse',
    action = 'store_true',
    help = 'output a reversed dictionary')

args = parser.parse_args()

#
# Main
#----------------------------------------------------------------------

steps = {} # dictionary to store sequences of increments

for line in args.file :

    if 'INCREMENT' in line : # .. STEP n INCREMENT k from a to b

        s = line.split()
        n, a, b = (int(x) for x in (s[-7], s[-3], s[-1])) # pull these values from line

        if n not in steps :
            steps[n] = set([])

        steps[n].add(a)
        steps[n].add(b)

# print dictionary at the end
if args.reverse :
    print_rev_dic(steps)
else :
    print_dic(steps)
