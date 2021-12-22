#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/3 下午3:18
# @Author  : Zhou Ran
# @File    : Faidx.py

import sys
import os
from collections import OrderedDict, namedtuple


class IndexRecord(
    namedtuple('IndexRecord',
               ['chrlength', 'offset', 'clen', 'blen'])):
    __slots__ = ()

    def __getitem__(self, key):
        if type(key) == str:
            return getattr(self, key)
        return tuple.__getitem__(self, key)

    def __str__(self):
        return "{chrlength:n}\t{offset:n}\t{clen:n}\t{blen:n}\n".format(**self._asdict())

    def __len__(self):
        return self.chrlength


class Fasta(object):
    """
    A class return the Fasta object
    """

    def __init__(self, chr, start, end, seq, strand, name=None):
        """

        :param chr:
        :param start:
        :param end:
        :param seq:
        :param strand:
        :param name:
        """
        assert isinstance(start, int)
        assert isinstance(end, int)
        self._chr = chr
        self._start = start
        self._end = end
        self._seq = seq
        self._strand = strand
        if name == None:
            self.name = ':'.join([self._chr,
                                  '-'.join(map(str, [self._start, self._end])),
                                  self._strand])
        else:
            self.name = name

    @property
    def seq(self):
        """
        an API to access original sequence
        :return:
        """
        return self._seq

    @property
    def realseq(self):
        """
        an API to access the 5'->3' sequence
        :return:
        """
        if self._strand == "+":
            return self.seq
        else:
            return complement(self.seq)[::-1]

    @property
    def reverse(self):
        """
        :return:
        """
        return self.__class__(self._chr, self._start, self._end, self.seq[::-1], self._strand)

    @property
    def complement(self):
        return self.__class__(self._chr, self._start, self._end, complement(self.seq), self._strand)

    def __add__(self, other):
        assert isinstance(other, Fasta), "<{}> can't add to <class: Fasta>".format(other.__class__)
        newid = ':'.join([self._chr, '-'.join(map(str, [self._start, self._end, other._start, other._end]))])
        seq = self.realseq + other.realseq
        return self.__class__(self._chr, self._start, other._end, seq, self._strand, name=newid)

    def __len__(self):
        return len(self._seq)

    def __repr__(self):
        return "<Fasta object: {}>".format(self.name)

    def __str__(self):
        seq = wrapperseq(self.seq)
        return "\n".join(['>{}'.format(self.name), seq])


class Faidx(object):

    def __init__(self, filename):
        self._filename = filename
        self._indexname = filename + '.fai'
        self.index = OrderedDict()
        try:
            self.file = open(self._filename, 'r')
        except:
            raise IOError("Can't read the Fasta files! "
                          "Please check whether the file exists or limited access")

        try:
            if os.path.exists(self._indexname):
                self.readindex()
            else:
                self.buildindex()
                self.readindex()
        except:
            raise IOError("Can't access the INDEX files! "
                          "Please check whether the file exists or limited access")

    def buildindex(self):
        """
        build a faidx file just like samtools faidx
        We need these information, site: https://www.biostars.org/p/11523/#11524
        Each line in this file contains:
        the name of the sequence
        the length of the sequence
        The byte offset where each chromosome sequence starts
        Number of bytes in each line
        Number of bases pairs appearing each line

        :return:
        """
        try:
            with open(self._filename, 'r') as Fa:
                with open(self._indexname, 'w') as Index:
                    chrname = None  # the name of the sequence
                    chrlength = 0  # chr length
                    offset = 0  # offset information
                    blen = None  # byte line length (includes '\n\r')
                    clen = None  # character line length
                    for i, line in enumerate(Fa):
                        line_blen = len(line)
                        line_clen = len(line.rstrip('\n\r'))
                        if line.startswith(">"):
                            if i != 0:
                                Index.write(
                                    '{}\t{}\t{}\t{}\t{}\n'.format(
                                        chrname, chrlength, current_offset, clen, blen))
                            chrname = line[1:].split(' ')[0]
                            chrlength = 0
                            offset += line_blen
                            current_offset = offset
                        else:
                            if not blen:
                                blen = line_blen
                            if not clen:
                                clen = line_clen

                            offset += line_blen
                            chrlength += line_clen

                    Index.write(
                        '{}\t{}\t{}\t{}\t{}\n'.format(
                            chrname, chrlength, offset, clen, blen))
        except:
            raise IOError("Can't read the files! "
                          "Please check whether the file exists or limited access")

    def readindex(self):

        """
        read index file
        :return:
        """
        with open(self._indexname) as indexfile:
            for line in indexfile.readlines():
                chrom, chrlength, offset, clen, blen = line.strip().split('\t')
                chrlength, offset, clen, blen = map(int, [chrlength, offset, clen, blen])
                self.index[chrom] = IndexRecord(chrlength, offset, clen, blen)

    def fetch(self, chr, start, end, strand):

        """
        use the 0-based coordinary for retrieving the sequence
        :param chr: chromosome name
        :param start: chromosome start site
        :param end: chromosome start end
        :return: Return the fasta sequence
        """

        assert isinstance(start, int), "The start site must be a int type, should not {}!\n" \
            .format(start.__class__)
        assert isinstance(end, int), "The end site must be a int type, should not {}!\n" \
            .format(end.__class__)

        start_ = start - 1

        try:
            chrinfo = self.index[chr]
        except:
            raise TypeError("The {} name can't find in Fasta file, pls check the chromosome or contig name.\n"
                            .format(chr))
        linebefore = int((start_) / chrinfo.clen * (chrinfo.blen - chrinfo.clen))
        seqlength = end - start_
        realstart = chrinfo.offset + start_ + linebefore * (chrinfo.blen - chrinfo.clen)

        if end >= chrinfo.chrlength:  # if the end is over the chromosome length, then return all sequence
            allseqlength = int(chrinfo.chrlength / chrinfo.clen * (chrinfo.blen - chrinfo.clen))
            seqlength = chrinfo.chrlength - start_ + allseqlength - linebefore
            self.file.seek(realstart)
            seq = self.file.read(seqlength).replace('\n', '')
            return Fasta(chr=chr, start=start, end=end, seq=seq, strand=strand)

        else:
            lineafter = int(end / chrinfo.clen * (chrinfo.blen - chrinfo.clen))
            lineinside = lineafter - linebefore

            seqlength = lineinside + seqlength
            self.file.seek(realstart)
            seq = self.file.read(seqlength).replace('\n', '')
            return Fasta(chr=chr, start=start, end=end, seq=seq, strand=strand)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.close()

    @staticmethod
    def getsplice(fai, chr, junctionlist, strand):
        """
        A function to return the junction sequence
        :param faidx:
        :param jucntionlist:
        :param strand:
        :return:
        """
        res = []
        for site in junctionlist:
            s, e = site
            res.append(fai.fetch(chr, s, e, strand))
        return res


def complement(seq):
    """
    Complement the sequence
    :param seq: "ATATCTC"
    :return: "TATGAG"
    """
    alpha = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }
    return ''.join([alpha[i] for i in list(seq)])


def wrapperseq(seq, length=60):
    """
    wrapper a sequence for beauty
    :param seq:
    :return:
    """
    final = []
    while len(seq) > length:
        final.append(seq[:length])
        seq = seq[length:]
    final.append(seq)
    return '\n'.join(final)


def main(file):
    test = Faidx(file)  # /Users/zhouran/opt/proj/2018-11-17-probe/rebuild/Probes_evaluation/testFaidx.fa
    seq = test.fetch('ENSMUST00000005815.6', 1, 6000, "+")
    print(seq.complement.reverse)


if __name__ == '__main__':
    main(sys.argv[1])
