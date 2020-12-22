import sys


def compare(align1, align2):
    if len(align1) != len(align2):
        return False
    for pair in align1:
        if pair not in align2:
            return False
    return True


def compareAlignemnts(f1, f2):
    curr1, curr2 = set(), set()
    for line1, line2 in zip(open(f1), open(f2)):
        if line1.startswith("seednum") and line2.startswith("seednum"):
            if not compare(curr1, curr2):
                return False
            curr1, curr2 = set(), set()
            continue
        curr1.add(tuple(line1.split()))
        curr2.add(tuple(line2.split()))
    return compare(curr1, curr2)


if __name__ == '__main__':
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    if not compareAlignemnts(f1, f2):
        print("The alignments in two files are not the same!")
        exit(1)
    print("The alignments in two files are the same!")
    exit(0)