#!/usr/bin/env python3

import csv
import random
import string

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

bases = ['A', 'C', 'G', 'T']

qual = [chr(i) for i in range(35, 75)]


# ------------------------------------------------------------------------------
# Utility functions
# ------------------------------------------------------------------------------

def read_ref(filename):
    with open(filename, newline = '') as fh:
        reader = csv.reader(fh, delimiter=',')
        return [(row[0], row[1]) for row in reader]


def read_ref1(filename):
    with open(filename, newline = '') as fh:
        return [line.strip() for line in fh.readlines()]


def seqs(ref):
    return [row[0] for row in ref]


def ids(ref):
    seen = {}
    ret = []
    for row in ref:
        id = row[1]
        if not id in seen:
            seen[id] = True
            ret.append(id)
    return ret


def read_barcodes(filename):
    with open(filename, 'rU') as fh:
        reader = csv.reader(fh, delimiter=',')
        return [row[0] for row in reader]


def random_string(alphabet, length):
    return ''.join([random.choice(alphabet) for i in range(0, length)])


def random_seq(length):
    return random_string(bases, length)


def random_length_seq(min_length, max_length):
    length = random.randint(min_length, max_length)
    return random_seq(length)


def random_qual(length):
    return random_string(qual, length)


def fastq_read_id(number):
    return '@HWUSI-EAS100R:6:23:398:3989#' + str(number)


def gaussian_qual(length, mean, var):
    '''
    returns a quality score sampled from a Gaussian distribution
    with the provided mean and standard deviation
    '''
    def trim(x):
        rounded = round(x)
        if rounded < 35:
            return 35
        elif rounded > 74:
            return 74
        else:
            return int(rounded)

    return [trim(random.gauss(mean, var)) for i in range(0, length)]


def random_barcodes(length, number):
    barcodes = {}
    while len(list(barcodes.keys())) < number:
        barcode = random_seq(length)
        barcodes[barcode] = 0
    return sorted(barcodes.keys())


def maybe_mutate_seq(seq, quals):
    chars = list(seq)
    for i in range(0, len(seq)):
        if quals[i] < 36 and random.random() < 0.7:
            chars[i] = 'N'
            return ''.join(chars)
    return seq


def single_fastq(barcode, seq, stuffer, qual_mean, qual_stddev, number):
    readid = fastq_read_id(number)
    read_seq = barcode + stuffer + seq
    quals = gaussian_qual(len(read_seq), qual_mean, qual_stddev)
    mutated_seq = maybe_mutate_seq(read_seq, quals)
    return readid, mutated_seq, quals


def print_fastq_record(file, fastq):
    (readid, seq, qual) = fastq
    file.write(readid + '\n')
    file.write(seq + '\n')
    file.write('+\n')
    file.write(''.join([chr(q) for q in qual]) + '\n')


def print_umi_fastq_files(filename1, filename2, data):
    with open(filename1, 'w') as file:
        stuffer = random_seq(18) + 'CACCG'
        for (i, (_, seq, umi)) in enumerate(data):
            stagger = random_length_seq(0, 7)
            stagger_len = len(stagger)
            fill_len = 7 - stagger_len
            fill = random_seq(fill_len)
            fastq = single_fastq('', seq + 'TATCCGT' + umi + fill, stagger + stuffer, 55, 5, i + 1)
            print_fastq_record(file, fastq)

    # print the barcode file
    with open(filename2, 'w') as file:
        for (i, (barcode, _, _)) in enumerate(data):
            fastq = single_fastq(barcode, '', '', 55, 5, i + 1)
            print_fastq_record(file, fastq)


def print_umi_counts_files(column_ref, row_ref, umi_bcs, data):
    # these store the scores in per-UMI and aggregated form
    per_umi_scores = {}
    unified_scores = {}

    # build up `per_umi_scores` and `unified_scores`
    for (c, r, u) in data:
        # start with `per_umi_scores`
        shard = per_umi_scores.get(u, {})
        count = shard.get((r, c), 0)
        count += 1
        shard[(r, c)] = count
        per_umi_scores[u] = shard

        # next work on `unified_scores`
        unified_count = unified_scores.get((r, c), 0)
        unified_count += 1
        unified_scores[(r, c)] = unified_count

    # this is a map of column IDs to a list of associated barcodes
    col_bcs_for_id = {}
    for (bc, id) in column_ref:
        bcs = col_bcs_for_id.get(id, [])
        bcs.append(bc)
        col_bcs_for_id[id] = bcs

    col_ids = ids(column_ref)
    header = ['Row Barcode', 'Row Barcode IDs'] + col_ids

    # write the UMI-barcode-specific counts files
    us = umi_bcs[:]
    us.append('TTTTT')
    for umi in us:
        fumi = umi
        if umi == 'TTTTT':
            fumi = 'UNMATCHED-UMI'
        scores = per_umi_scores.get(umi, {})
        with open('counts-' + fumi + '.txt', 'w') as file:
            print_scores_to_fh(col_bcs_for_id, row_ref, col_ids, header, scores, file)

    # write the aggregate counts file
    with open('counts.txt', 'w') as file:
        print_scores_to_fh(col_bcs_for_id, row_ref, col_ids, header, unified_scores, file)


def print_scores_to_fh(col_bcs_for_id, row_ref, col_ids, header, scores, fh):
    # write the header
    fh.write('\t'.join(header) + '\n')
    for (row_bc, row_id) in row_ref:
        counts = {}
        for (col_id, bcs) in list(col_bcs_for_id.items()):
            count = sum([scores.get((row_bc, col_bc), 0) for col_bc in bcs])
            counts[col_id] = count
        row_data = [str(counts.get(id, 0)) for id in col_ids]
        fh.write('\t'.join([row_bc, row_id] + row_data) + '\n')


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------


if __name__ == '__main__':
    column_ref = read_ref('Conditions.csv')
    row_ref = read_ref('Reference.csv')
    umi_bcs = read_ref1('Umi.csv')

    col_bcs = seqs(column_ref)
    row_bcs = seqs(row_ref)

    data = [(random.choice(col_bcs), random.choice(row_bcs), random.choice(umi_bcs)) for _ in range(0, 964)]
    data += [(random.choice(col_bcs), random.choice(row_bcs), 'TTTTT') for _ in range(0, 36)]
    random.shuffle(data)

#    print_scenario1_fastq_file('scenario1/scenario1.fastq', data)
#    print_scenario2_fastq_file('scenario2/scenario2.fastq', data)
#    print_scenario3_fastq_file('scenario3/scenario3.1.fastq', 'scenario3/scenario3.barcode_1.fastq', data)
#    print_scenario4_fastq_file('scenario4/scenario4.1.fastq', 'scenario4/scenario4.barcode_1.fastq', data)
    print_umi_fastq_files('umi1.1.fastq', 'umi1.barcode_1.fastq', data)
    print_umi_counts_files(column_ref, row_ref, umi_bcs, data)
