#!/bin/python3
import os

def assert_with_prints(value, target, value_title):
    assert value == target, f'{value_title} is not {target}'
    print(f'success, {value_title} = {target}')

def calc_f1(orth, found, n):
    tp = orth
    fp = found - orth
    fn = n - orth
    return tp / (tp + 0.5 * (fp + fn))

def get_seeding_dir():
    return os.getenv('BLANT_SEED_SUPPL')

if __name__ == '__main__':
    assert_with_prints(5, 5, 'foo')
