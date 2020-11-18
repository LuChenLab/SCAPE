from collections import Counter


def cigar_support(line, strand):
    """
    Processing the cigar information.
    [(cigar, length)], soft clip (4)
    """

    cigarinfo = line.cigar
    if cigarinfo[-1][0] == 4 and strand == '+':
        cigar_c = int(cigarinfo[-1][1])
        """
        For R2, if the length of soft clip was more than 30,
        we check the polyA[T] from 5' ot 3'
        """

        if cigar_c >= 30:
            sequence_ = line.query_sequence[-cigar_c:-cigar_c + 29]
        else:
            sequence_ = line.query_sequence[-cigar_c:]

        if cigar_c < 6:
            pa_support = 'no'
        else:
            freq = Counter(sequence_)
            if freq['A'] > len(sequence_) * 0.8:
                pa_support = 'yes'
            else:
                pa_support = 'no'

    elif cigarinfo[0][0] == 4 and strand == '-':
        cigar_c = int(cigarinfo[0][1])
        if cigar_c >= 30:
            sequence_ = line.query_sequence[cigar_c-30:cigar_c]
        else:
            sequence_ = line.query_sequence[:cigar_c]

        if cigar_c < 6:
            pa_support = 'no'
        else:
            freq = Counter(sequence_)
            if freq['T'] > len(sequence_) * 0.8:
                pa_support = 'yes'
            else:
                pa_support = 'no'

    else:
        pa_support = 'no'

    return pa_support


def collapse_intron(read):
    """
    return 1-based site, now only process the R2
    """

    end_site = read.reference_end if not read.is_reverse else read.reference_start

    cigar = read.cigar
    exon_bound = []
    for c, l in cigar:
        if c == 0:
            exon_bound.append(l)
        else:
            continue
    if read.is_reverse:
        fake_start = end_site + sum(exon_bound) - 1
    else:
        fake_start = end_site - sum(exon_bound) - 1 

    return fake_start + 1


def check_strand(read, strand):
    """
    check the read strand to remove of opposite reads,
    this is only suitable for fr-firststrand library
    """

    if strand == '+':
        if read.is_paired:
            if read.is_read1 and read.is_reverse:
                return True
            elif read.is_read2 and not read.is_reverse:
                return True
            else:
                return False
        else:
            if read.is_reverse:
                return False
            else:
                return True

    elif strand == '-':
        if read.is_paired:
            if read.is_read1 and not read.is_reverse:
                return True
            elif read.is_read2 and read.is_reverse:
                return True
            else:
                return False
        else:
            if read.is_reverse:
                return True
            else:
                return False

    else:
        raise ValueError(f'''Unrecognized strand: {strand} label, please check it.''')
