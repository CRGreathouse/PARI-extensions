#!/usr/bin/env python3

import re
import subprocess
from typing import List
from array import array
import argparse

# Set these to the locations of these files on your system.
src = '../../pari/src/' # '../src'
whatnow = f'{src}/gp/whatnow.h'
paridesc = f'{src}/desc/pari.desc'
gp2c = '/usr/local/bin/gp2c'

parser = argparse.ArgumentParser(description='Show GP functions')
parser.add_argument('--format', choices=['1', 'comma'], default='1', help='1 (one function per line), comma (comma-separated list)')
parser.add_argument('--show', choices=['good', 'bad', 'removed', 'obsolete', 'current', 'all'], default='good', help='good (current functions which are not obsolete), obsolete (current functions which are slated to be removed), removed (functions no longer available), bad (removed or obsolete functions), current (all available functions, whether obsolete or not), all (all functions, including those removed)')
args = parser.parse_args()

def array_minus(a: list, b: list) -> list:
    return list(set(a) - set(b))

def unique(a: list) -> list:
    return list(set(a))


# Get raw output from gp2c.
def run_gp2c() -> List[str]:
    env = {'PATH': '/bin:/usr/bin'}
    proc = subprocess.run([gp2c, '-l'], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return proc.stdout.decode().splitlines()


# Get a list of functions from gp2c.
def get_functions() -> List[str]:
    func = []
    operator = []  # currently ignored
    for line in run_gp2c():
        if match := re.match(r'^([^ ]+) ([^ ]+) (.*)$', line, re.IGNORECASE):
            GP, PARI, prototype = match.groups()
            if re.match(r'^[a-z][a-z0-9_]*$', GP, re.IGNORECASE):
                func.append(GP)
            else:
                operator.append(GP)
        elif re.match(r'^[\s\r\n]*$', line):
            pass
        else:
            print(f"Unrecognized gp2c line: '{line}'", file=sys.stderr)
    return sorted(func)


def get_removed() -> List[str]:
    removed = []
    with open(whatnow, 'r') as f:
        for line in f:
            if match := re.match(r'^\{"[^" ]+", *_SAME\}', line):
                pass
            elif match := re.match(r'^\{"([a-zA-Z][a-zA-Z0-9_]*)", *_REMOV\}', line):
                removed.append(match.group(1))  # Explicitly marked as removed
            elif match := re.match(r'\{"([a-zA-Z][a-zA-Z0-9_]*)","([^" ]+)","[^"]*","[^"]*"\}', line):
                oldname, newname = match.groups()
                if oldname != newname:
                    removed.append(oldname)
    return sorted(removed)


def get_obsolete() -> List[str]:
    obsolete = []
    name = ''
    with open(paridesc, 'r') as f:
        for line in f:
            if match := re.match(r'^Function: ([a-zA-Z][a-zA-Z0-9_]*)$', line):
                name = match.group(1)  # Function
            elif match := re.match(r'^Function: ([^a-zA-Z]+)$', line):
                name = match.group(1)  # Operator
            elif match := re.match(r'^Function: (_[a-zA-Z0-9_]+)$', line):
                name = match.group(1)  # gp2c internal function
            elif match := re.match(r'^Function: (_\.[a-zA-Z][a-zA-Z0-9_]*)$', line):
                name = match.group(1)  # member function
            elif re.match(r'^Function:', line):  # weird, skipped
                # Currently this includes only 'O(_^_)' and '_^s'
                name = ''
            elif match := re.match(r'^Obsolete: (.+)$', line):
                obsoleteDate = match.group(1)  # currently ignored
                if name != '':
                    obsolete.append(name)
    return sorted(obsolete)


def show(out, fmt):
    if fmt == '1':
       print('\n'.join(out))
    elif fmt == 'comma':
        print(', '.join(out))

func = get_functions()
removed = get_removed()
obsolete = get_obsolete()

# If it's in both, one function was removed and another with the same name added.
removed = list(set(removed) - set(func))

bad = sorted(list(set(removed) | set(obsolete)))
good = sorted(list(set(func) - set(obsolete)))

if args.show == 'good':
    show(good, args.format)
elif args.show == 'bad':
    show(bad, args.format)
elif args.show == 'removed':
    show(removed, args.format)
elif args.show == 'obsolete':
    show(obsolete, args.format)
elif args.show == 'current':
    show(func, args.format)
elif args.show == 'all':
    if args.format == '1':
        all = sorted(list(set(good) | set(bad)))
        show(all, '1')
    else:
        for s in ['good', 'bad', 'removed', 'obsolete']:
            print(f"{s}:")
            show(s, args.format)

#print(f"In total there are {len(good)} good, {len(removed)} removed, and {len(obsolete)} obsolete functions.")
