#!/usr/bin/env python3
#
# Prints the accuracy and runtime results written by the Makefile targets as
# Markdown tables.
#
# Usage: report.py <build directory> <variant>...

import csv
import sys
from pathlib import Path

ACCURACY_FIELDS = ['correct_pct', 'mean_ulp', 'max_ulp', 'above_1ulp_pct']


def read_accuracy(build, variant):
    with open(build / f'accuracy_{variant}.csv') as f:
        return {(row['type'], int(row['n']), row['distribution'],
                 row['function']): row
                for row in csv.DictReader(f)}


def read_benchmarks(build, variant):
    # Keeps the best time over all rounds.
    times = {}
    for path in sorted(build.glob(f'benchmark_*_{variant}_*.csv')):
        compiler = path.name.split('_')[1]
        with open(path) as f:
            for row in csv.DictReader(f):
                for function in ['norm', 'rnorm']:
                    key = (compiler, row['type'], int(row['n']), row['range'],
                           function)
                    value = float(row[f'{function}_ns'])
                    times[key] = min(times.get(key, value), value)
    return times


def format_accuracy(row):
    return (f'{float(row["correct_pct"]):.2f}% / '
            f'{float(row["mean_ulp"]):.3f} / {row["max_ulp"]}')


def print_accuracy(build, variants):
    tables = {variant: read_accuracy(build, variant) for variant in variants}
    reference = tables[variants[0]]

    # Variants with identical accuracy are reported once.
    distinct = [variants[0]]
    for variant in variants[1:]:
        if all(tables[variant][key][field] == reference[key][field]
               for key in reference for field in ACCURACY_FIELDS):
            print(f'The accuracy of `{variant}` is identical to the accuracy '
                  f'of `{variants[0]}` in all configurations.\n')
        else:
            distinct.append(variant)

    for type_ in ['double', 'float']:
        for function in ['norm', 'rnorm']:
            print(f'#### {type_} {function}\n')
            print('| n | distribution | '
                  + ' | '.join(distinct) + ' |')
            print('|---|---|' + '---|' * len(distinct))
            keys = sorted((key for key in reference
                           if key[0] == type_ and key[3] == function),
                          key=lambda key: (key[1], key[2]))
            for key in keys:
                cells = [format_accuracy(tables[variant][key])
                         for variant in distinct]
                print(f'| {key[1]} | {key[2]} | ' + ' | '.join(cells) + ' |')
            print()


def print_benchmarks(build, variants):
    tables = {variant: read_benchmarks(build, variant) for variant in variants}
    reference = tables[variants[0]]
    for compiler in sorted({key[0] for key in reference}):
        print(f'#### {compiler}\n')
        print('| type | n | range | function | ' + ' | '.join(variants) + ' |')
        print('|---|---|---|---|' + '---|' * len(variants))
        keys = sorted((key for key in reference if key[0] == compiler),
                      key=lambda key: (key[1] != 'double', key[4], key[3],
                                       key[2]))
        for key in keys:
            cells = [f'{tables[variant][key]:.2f}' for variant in variants]
            print(f'| {key[1]} | {key[2]} | {key[3]} | {key[4]} | '
                  + ' | '.join(cells) + ' |')
        print()


def main():
    build = Path(sys.argv[1])
    variants = sys.argv[2:]
    print('### Accuracy\n')
    print('Each cell shows the share of correctly rounded results, the mean '
          'error, and the maximum error in ULP.\n')
    print_accuracy(build, variants)
    print('### Runtime\n')
    print('Each cell shows the best time per call in nanoseconds.\n')
    print_benchmarks(build, variants)


if __name__ == '__main__':
    main()
