from collections import defaultdict
from utils.utils import wrapdict
from utils.gtf import GTFfile

class Transcripts(GTFfile):

    def __init__(self, file):
        GTFfile.__init__(self, file)

    def structure(self):
        gene_structure = defaultdict(wrapdict)

        skip = [
            'stop_codon',
            'start_codon',
            'gene'
        ]

        for line in GTFfile.__iter__(self):

            if line.feature in set(skip):
                continue

            tid = line.attributes["transcript_id"]
            gid = line.attributes["gene_id"]

            if line.feature == 'transcript':
                try:
                    genename = line.attributes['gene_name']
                except:
                    genename = 'None'
                gene_structure[gid][tid]["chrom"] = line.chrom
                gene_structure[gid][tid]["gene_name"] = genename

                if "gene_biotype" in line.attributes:
                    gene_structure[gid][tid]["gene_biotype"] = line.attributes["gene_biotype"]
                elif "gene_type" in line.attributes:
                    gene_structure[gid][tid]["gene_type"] = line.attributes["gene_type"]


                if "transcript_biotype" in line.attributes:
                    gene_structure[gid][tid]["transcript_biotype"] = line.attributes["transcript_biotype"]
                elif "transcript_type" in line.attributes:
                    gene_structure[gid][tid]["transcript_type"] = line.attributes["transcript_type"]

                gene_structure[gid][tid]["strand"] = line.strand
                gene_structure[gid][tid][line.feature].append(
                    (line.start, line.end)
                    )

            else:
                gene_structure[gid][tid][line.feature].append(
                    (line.start, line.end)
                    )
        return gene_structure

class Iso:
    def __init__(self,
                 tid,
                 gid,
                 chrom,
                 strand,
                 region,
                 exonregion=None,
                 genename=None,
                 biotype=None,
                 cdsregion=None,
                 five_utr=None,
                 three_utr=None):
        self._tid = tid
        self._gid = gid
        self._chrom = chrom
        self._strand = strand
        self._biotype = biotype
        self._shortname = genename
        self._region = region
        self._exonregion = sorted(exonregion, key=lambda x: x[0]) if exonregion else exonregion
        self._cdsregion = sorted(cdsregion, key=lambda x: x[0]) if cdsregion else cdsregion
        self._five_utr = sorted(five_utr, key=lambda x: x[0]) if five_utr else five_utr
        self._three_utr = sorted(three_utr, key=lambda x: x[0]) if three_utr else three_utr

    @property
    def tid(self):
        return self._tid

    @property
    def chrom(self):
        return self._chrom

    @property
    def gid(self):
        return self._gid

    @property
    def strand(self):
        return self._strand

    @property
    def biotype(self):
        return self._biotype

    @property
    def region(self):
        return self._region

    @property
    def exon(self):
        return self._exonregion

    @property
    def cds(self):
        return self._cdsregion

    @property
    def five_utr(self):
        return self._five_utr

    @property
    def three_utr(self):
        return self._three_utr

    @property
    def geneshortname(self):
        return self._shortname

    @property
    def intron(self):
        intronlist = []
        for window in list(zip(self.exon[:-1], self.exon[1:])):
            intronlist.extend([[window[0][1] + 1, window[1][0] - 1]])
        return intronlist

    @property
    def exon_len(self):
        return list(map(lambda x: x[1] - x[0] + 1, self.exon))

    def truncate(self, threshold=1000):

        def decide_region(exon_list, strand, threshold):
            if strand == "+":
                exon_list = exon_list[::-1]

            num_list = list(map(lambda x: x[1] - x[0] + 1, exon_list))
            t = 0
            region = []

            for index, val in enumerate(num_list):
                t += val
                if t > threshold:
                    offset = t - threshold
                    if strand == "+":
                        offsetregion = tuple([exon_list[index][0] + offset, exon_list[index][1]])
                        region.append(offsetregion)
                    else:
                        offsetregion = tuple([exon_list[index][0], exon_list[index][1] - offset])
                        region.append(offsetregion)
                    break
                elif t == threshold:
                    region.append(exon_list[index])
                    break
                else:
                    region.append(exon_list[index])

            region = sorted(region, key=lambda x: x[0])
            return (region)

        if sum(self.exon_len) <= threshold:
            return (self.exon)

        elif len(self.exon_len) == 1:
            if self.strand == "+":
                region = [(self.exon[0][0], self.exon[0][0] + threshold - 1)]

            else:
                region = [(self.exon[0][1] - threshold + 1, self.exon[0][1])]
        else:
            region = decide_region(self.exon, self.strand, threshold)
        if region == []:
            print(self.tid, self.strand, self.exon, self.exon_len)
        return (region)

    def __repr__(self):
        return 'TrancriptId: {}, \
        GeneID: {}, \GeneName: {}, \
        ExonNum: {}, CDSNum: {}, \
        IntronNum: {}'.format(self.tid,
                              self.gid,
                              self.geneshortname,
                              len(
                                  self.exon) if self.exon else 0,
                              len(
                                  self.cds) if self.cds else 0,
                              len(
                                  self.intron) if self.intron else 0)

