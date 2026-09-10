"""Half-open interval transforms with explicit chain coverage and ambiguity."""
from __future__ import annotations
import bisect
import gzip
from collections import defaultdict


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


class ChainMap:
    def __init__(self, path, source_chrom):
        self.chains=[]
        with gzip.open(path,'rt') as handle:
            current=None
            for line in handle:
                fields=line.split()
                if not fields: continue
                if fields[0]=='chain':
                    current=None
                    if fields[2]!=source_chrom: continue
                    if fields[4]!='+': raise ValueError('Unsupported source strand')
                    current=dict(id=fields[12],score=int(fields[1]),chrom=fields[7],
                                 size=int(fields[8]),strand=fields[9],blocks=[],starts=[])
                    self.chains.append(current)
                    source=int(fields[5]); target=int(fields[10])
                elif current is not None:
                    size=int(fields[0]);current['blocks'].append((source,source+size,target))
                    current['starts'].append(source)
                    if len(fields)==3:
                        source+=size+int(fields[1]);target+=size+int(fields[2])

    def map_interval(self,start,end,min_coverage=0.95):
        if start<0 or end<=start: raise ValueError('Invalid half-open interval')
        candidates=[]
        for chain in self.chains:
            blocks=chain['blocks'];i=max(0,bisect.bisect_right(chain['starts'],start)-1)
            mapped=[];coverage=0
            for left,right,qstart in blocks[i:]:
                if left>=end:break
                a,b=max(left,start),min(right,end)
                if b<=a:continue
                q0=qstart+a-left;q1=qstart+b-left
                if chain['strand']=='-':q0,q1=chain['size']-q1,chain['size']-q0
                mapped.append((q0,q1));coverage+=b-a
            if coverage/(end-start)>=min_coverage and mapped:
                candidates.append(dict(chrom=chain['chrom'],start=min(x[0] for x in mapped),
                     end=max(x[1] for x in mapped),strand=chain['strand'],
                     coverage=coverage/(end-start),chain_id=chain['id'],score=chain['score'],
                     aligned_blocks=mapped))
        # Different chains mapping to exactly the same interval are not distinct mappings.
        unique={}
        for item in candidates:
            key=(item['chrom'],item['start'],item['end'],item['strand'])
            if key not in unique or item['score']>unique[key]['score']:unique[key]=item
        if len(unique)!=1:
            return dict(status='unmapped' if not unique else 'ambiguous',candidate_count=len(unique))
        return dict(status='mapped',**next(iter(unique.values())))


class IntervalIndex:
    """Index merged intervals for boolean overlaps; touching is not overlap."""
    def __init__(self, intervals):
        merged=[]
        for a,b in sorted(intervals):
            if b<=a:raise ValueError('Invalid indexed interval')
            if merged and a<=merged[-1][1]:merged[-1][1]=max(merged[-1][1],b)
            else:merged.append([a,b])
        self.intervals=merged;self.starts=[a for a,b in merged];self.ends=[b for a,b in merged]

    def overlaps(self,start,end):
        if end<=start:return False
        i=bisect.bisect_left(self.starts,end)-1
        return i>=0 and self.ends[i]>start


def bed_rows(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith(('#','track','browser')) or not line.strip():continue
            a=line.rstrip('\n').split('\t')
            if len(a)>=3:yield a


def read_fasta(path):
    with gzip.open(path,'rt') as f:
        return ''.join(line.strip() for line in f if not line.startswith('>')).upper()
