#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/10/23 9:34 下午
__author__ = 'Zhou Ran'

import os
import gzip
import pickle

from utils.utils import AttrDict


class GTFFeature:
    def __init__(self,
                 chrom=None,
                 source=None,
                 featuretype=None,
                 start=None,
                 end=None,
                 score=None,
                 strand=None,
                 phase=None,
                 attributes=None):
        self._chrom = chrom
        self._source = source
        self._featuretype = featuretype
        self._start = start
        self._end = end
        self._score = score
        self._strand = strand
        self._attributparse(attributes)

    def _attributparse(self, attributes):
        self.attributes = AttrDict()
        if attributes:
            items = attributes.split(';')
            for i in items:
                if len(i) == 0:
                    continue
                try:
                    name, key = i.strip().split()
                except ValueError:
                    continue
                key = key.replace('"', '')
                setattr(self.attributes, name, key)
        return self.attributes

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return int(self._start)

    @property
    def end(self):
        return int(self._end)

    @property
    def feature(self):
        return self._featuretype

    @property
    def strand(self):
        return self._strand

    def __len__(self):
        length = int(self._end) - int(self._start) + 1
        if length < 1:
            raise ValueError('Zero- or negative length feature')
        return length


class GTFfile:
    GTFFeature = GTFFeature

    def __init__(self, file):
        """

        :param file:
        """
        if os.path.splitext(file)[-1] == '.gz':
            self.file = gzip.open(file)
        else:
            self.file = open(file)

    def __iter__(self):
        """
        Yield a feature for each line
        """
        f = self.file
        for line in f:
            line = line.decode('utf-8') if isinstance(line, bytes) else line
            if line.startswith('#') or len(line) == 0:
                continue
            args = line.strip().split('\t')
            yield self.__class__.GTFFeature(*args)
        f.close()

    def __repr__(self):
        return 'GFFfile: %s' % (self.file)
