"""Frozen retrospective genomic classification; no causal or donor inference."""
from __future__ import annotations
import bisect
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .coordinates import IntervalIndex, bed_rows, read_fasta, reverse_complement
from .resources import ROOT, CACHE, sha256
from .human import table

BENCH=ROOT/'benchmark';RESULTS=BENCH/'results';RESULTS.mkdir(exist_ok=True)
PROTOCOL=json.loads((BENCH/'protocol.json').read_text())
BASECODE=np.full(256,4,dtype=np.int8)
for i,base in enumerate(b'ACGT'):BASECODE[base]=i


def motifs():
    out=[]
    for name in ['MA0442.3','MA0472.2','MA0090.4']:
        matrix=json.loads((CACHE/(name+'.json')).read_text())
        counts=np.array([matrix['pfm'][b] for b in 'ACGT'],dtype=float).T
        probs=(counts+0.5)/(counts.sum(axis=1,keepdims=True)+2)
        out.append(np.log2(probs/0.25))
    return out


def feature_sequence(seq,pwms):
    codes=BASECODE[np.frombuffer(seq.encode(),dtype=np.uint8)]
    if (codes>3).any():raise ValueError('Ambiguous sequence passed feature filter')
    gc=float(np.isin(codes,[1,2]).mean());cpg=seq.count('CG')/(len(seq)-1)
    values=[gc,cpg]
    for pwm in pwms:
        windows=np.lib.stride_tricks.sliding_window_view(codes,len(pwm))
        fw=pwm[np.arange(len(pwm)),windows].sum(axis=1)
        rc=pwm[::-1,::-1][np.arange(len(pwm)),windows].sum(axis=1)
        values.append(float(max(fw.max(),rc.max())))
    windows=np.lib.stride_tricks.sliding_window_view(codes,4).astype(np.int32)
    index=(windows*np.array([64,16,4,1])).sum(axis=1)
    rc_index=((3-windows[:,::-1])*np.array([64,16,4,1])).sum(axis=1)
    count=np.bincount(np.minimum(index,rc_index),minlength=256)/len(index)
    return np.r_[values,count]


def dataset():
    rng=np.random.default_rng(PROTOCOL['seed']);pwms=motifs()
    candidate=defaultdict(list);exclusions=defaultdict(list);black=defaultdict(list)
    for a in bed_rows(CACHE/'ENCFF602YVR.bed.gz'):
        candidate[a[0]].append((int(a[1]),int(a[2]),float(a[6]),int(a[1])+int(a[9])))
    for a in bed_rows(CACHE/'ENCFF632FMT.bed.gz'):
        exclusions[a[0]].append((max(0,int(a[1])-1000),int(a[2])+1000))
    for a in bed_rows(CACHE/'hg38-blacklist.v2.bed.gz'):
        black[a[0]].append((int(a[1]),int(a[2])))
    records=[];features=[];seen=set();qc=[]
    for split in ['train','validation','test']:
        for chrom in PROTOCOL[split+'_chromosomes']:
            print('Preparing',split,chrom,flush=True)
            seq=read_fasta(CACHE/f'hg38.{chrom}.fa.gz')
            bad=IntervalIndex(black[chrom]);peaks=IntervalIndex(exclusions[chrom])
            genes=json.loads((CACHE/f'gencode50_{chrom}.json').read_text())['wgEncodeGencodeBasicV50']
            tss=sorted({t['txStart'] if t['strand']=='+' else t['txEnd']-1 for t in genes})
            # One representative summit in each connected overlapping 600 bp window.
            clusters=[]
            for start,end,signal,summit in sorted(candidate[chrom],key=lambda x:x[3]):
                if summit<0:continue
                window=(summit-300,summit+300,signal,summit)
                if clusters and window[0]<clusters[-1]['end']:
                    clusters[-1]['end']=max(clusters[-1]['end'],window[1])
                    if signal>clusters[-1]['best'][2]:clusters[-1]['best']=window
                else:clusters.append(dict(end=window[1],best=window))
            options=[x['best'] for x in clusters];rng.shuffle(options)
            used=[];positives=0;counts={'ambiguous':0,'duplicate_sequence':0,'blacklist':0}
            def add(start,end,label):
                if start<0 or end>len(seq):return False
                if bad.overlaps(start,end):counts['blacklist']+=1;return False
                s=seq[start:end]
                if any(b not in 'ACGT' for b in s):counts['ambiguous']+=1;return False
                digest=hashlib.sha256(min(s,reverse_complement(s)).encode()).hexdigest()
                if digest in seen:counts['duplicate_sequence']+=1;return False
                if any(a<end and b>start for a,b in used):return False
                seen.add(digest);used.append((start,end))
                center=(start+end)//2;i=bisect.bisect_left(tss,center)
                distance=min(abs(center-t) for t in tss[max(0,i-1):i+1])
                feat=feature_sequence(s,pwms)
                records.append(dict(example_id=f'{chrom}:{start}-{end}',chrom=chrom,start=start,end=end,
                    split=split,label=label,sequence_sha256=digest,gc_fraction=feat[0],
                    distance_to_tss=distance,block=f'{chrom}:{center//1000000}'))
                features.append(np.r_[np.log1p(distance),feat]);return True
            for start,end,signal,summit in options:
                if add(start,end,1):positives+=1
                if positives>=PROTOCOL['max_positive_windows_per_chromosome']:break
            negatives=0;attempts=0
            while negatives<positives:
                attempts+=1
                if attempts>1000000:raise RuntimeError('Could not sample negative windows')
                start=int(rng.integers(1000,len(seq)-1600));end=start+600
                if peaks.overlaps(start,end):continue
                if add(start,end,0):negatives+=1
            qc.append(dict(chrom=chrom,split=split,positive=positives,negative=negatives,
                           source_rows=len(candidate[chrom]),merged_windows=len(clusters),**counts))
    frame=pd.DataFrame(records);x=np.array(features)
    # Version 1 was frozen with CRLF. Preserve those bytes on every platform.
    split_bytes=frame.to_csv(sep='\t',index=False,lineterminator='\r\n').encode('utf-8')
    split_path=BENCH/'splits.tsv'
    if split_path.exists() and split_path.read_bytes()!=split_bytes:
        raise ValueError('Frozen split changed; create a separately versioned benchmark instead of overwriting')
    if (RESULTS/'dataset_lock.json').exists():
        previous=json.loads((RESULTS/'dataset_lock.json').read_text())
        if previous['protocol_sha256']!=sha256(BENCH/'protocol.json'):
            raise ValueError('Frozen protocol changed; use a new benchmark version')
    split_path.write_bytes(split_bytes)
    table(RESULTS/'sampling_qc.tsv',qc)
    np.savez_compressed(CACHE/'benchmark_features.npz',x=x)
    stamp=dict(protocol_sha256=sha256(BENCH/'protocol.json'),split_sha256=sha256(BENCH/'splits.tsv'),
               examples=len(frame),feature_count=x.shape[1],seed=PROTOCOL['seed'])
    (RESULTS/'dataset_lock.json').write_text(json.dumps(stamp,indent=2)+'\n')
    return frame,x


def metrics(y,p):
    if len(set(y))<2:return dict(average_precision=None,AUROC=None,Brier_score=None)
    return dict(average_precision=float(average_precision_score(y,p)),AUROC=float(roc_auc_score(y,p)),
                Brier_score=float(brier_score_loss(y,p)))


def bootstrap(frame,preds):
    rng=np.random.default_rng(PROTOCOL['seed']+1)
    blocks=frame['block'].unique();indices={b:np.where(frame['block'].to_numpy()==b)[0] for b in blocks}
    scores=defaultdict(list);deltas=defaultdict(list);y=frame.label.to_numpy()
    for _ in range(500):
        ix=np.concatenate([indices[b] for b in rng.choice(blocks,len(blocks),replace=True)])
        if len(np.unique(y[ix]))<2:continue
        ap={name:average_precision_score(y[ix],p[ix]) for name,p in preds.items()}
        for name,value in ap.items():scores[name].append(value);deltas[name].append(value-ap['GC_CpG'])
    return [dict(model=name,genomic_blocks=len(blocks),bootstrap_draws=len(values),
        AP_low=float(np.quantile(values,.025)),AP_high=float(np.quantile(values,.975)),
        AP_delta_vs_GC_low=float(np.quantile(deltas[name],.025)),AP_delta_vs_GC_high=float(np.quantile(deltas[name],.975)))
        for name,values in scores.items()]


def main():
    frame,x=dataset();train=frame.split.eq('train').to_numpy();val=frame.split.eq('validation').to_numpy()
    test=frame.split.eq('test').to_numpy();y=frame.label.to_numpy()
    feature_sets={'GC_CpG':[1,2],'distance_to_annotated_TSS':[0],
                  'SOX10_EGR2_TEAD1_motifs':[3,4,5],'4mer_sequence_logistic':list(range(1,x.shape[1]))}
    preds={'prevalence':np.full(test.sum(),y[train].mean())};fitlog=[];models={};validation=[]
    for name,columns in feature_sets.items():
        best=None
        for c in PROTOCOL['regularization_C_grid']:
            model=make_pipeline(StandardScaler(),LogisticRegression(C=c,max_iter=2000,random_state=PROTOCOL['seed'],solver='lbfgs'))
            model.fit(x[train][:,columns],y[train]);score=average_precision_score(y[val],model.predict_proba(x[val][:,columns])[:,1])
            validation.append(dict(model=name,C=c,validation_AP=score))
            if best is None or score>best[0]:best=(score,c,model)
        score,c,model=best;models[name]=(columns,model)
        preds[name]=model.predict_proba(x[test][:,columns])[:,1]
        fitlog.append(dict(model=name,C=c,validation_AP=score))
        print('Fit',name,'C',c,'validation AP',round(score,4),flush=True)
    # Label-permutation sanity control, with its own validation-only C selection.
    perm=np.random.default_rng(PROTOCOL['seed']+2).permutation(y[train]);columns=feature_sets['4mer_sequence_logistic'];best=None
    for c in PROTOCOL['regularization_C_grid']:
        m=make_pipeline(StandardScaler(),LogisticRegression(C=c,max_iter=2000,random_state=PROTOCOL['seed']))
        m.fit(x[train][:,columns],perm);score=average_precision_score(y[val],m.predict_proba(x[val][:,columns])[:,1])
        if best is None or score>best[0]:best=(score,c,m)
    preds['permuted_training_labels']=best[2].predict_proba(x[test][:,columns])[:,1]
    testframe=frame.loc[test].copy().reset_index(drop=True);metricrows=[];calibration=[]
    for name,p in preds.items():
        testframe[name]=p;metricrows.append(dict(model=name,subset='all_test',n=len(p),**metrics(y[test],p)))
        for chrom in PROTOCOL['test_chromosomes']:
            mask=testframe.chrom.eq(chrom).to_numpy();metricrows.append(dict(model=name,subset=chrom,n=int(mask.sum()),**metrics(y[test][mask],p[mask])))
        mask=~(testframe.chrom.eq('chr17')&testframe.start.between(14229778,16265360)).to_numpy()
        metricrows.append(dict(model=name,subset='outside_PMP22_neighborhood',n=int(mask.sum()),**metrics(y[test][mask],p[mask])))
        for b in range(10):
            mask=(p>=b/10)&(p<(b+1)/10 if b<9 else p<=1)
            if mask.any():calibration.append(dict(model=name,bin_start=b/10,n=int(mask.sum()),mean_prediction=float(p[mask].mean()),observed_fraction=float(y[test][mask].mean())))
    # GC-matched assessment only: match positives and negatives within 0.02 GC bins per chromosome.
    matched=[];rng=np.random.default_rng(PROTOCOL['seed']+3)
    testframe['gc_bin']=(testframe.gc_fraction/0.02).astype(int)
    for (_,group) in testframe.groupby(['chrom','gc_bin']):
        pos=group.index[group.label.eq(1)].to_numpy();neg=group.index[group.label.eq(0)].to_numpy();n=min(len(pos),len(neg))
        if n:matched.extend(rng.choice(pos,n,replace=False));matched.extend(rng.choice(neg,n,replace=False))
    matched=np.array(sorted(matched))
    for name,p in preds.items():metricrows.append(dict(model=name,subset='GC_matched_test',n=len(matched),**metrics(y[test][matched],p[matched])))
    table(RESULTS/'metrics.tsv',metricrows);table(RESULTS/'validation.tsv',validation);table(RESULTS/'selected_models.tsv',fitlog)
    table(RESULTS/'calibration.tsv',calibration);table(RESULTS/'bootstrap.tsv',bootstrap(testframe,preds))
    testframe['in_gc_matched_subset']=testframe.index.isin(matched)
    testframe.to_csv(RESULTS/'predictions.tsv',sep='\t',index=False)
    # Store coefficients and scaling, avoiding executable pickle serialization.
    params={}
    for name,(cols,model) in models.items():
        scaler,lr=model.steps[0][1],model.steps[1][1]
        params[name]=dict(columns=cols,mean=scaler.mean_.tolist(),scale=scaler.scale_.tolist(),
                          coefficients=lr.coef_[0].tolist(),intercept=float(lr.intercept_[0]))
    (RESULTS/'model_parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    print(pd.DataFrame(metricrows).query("subset == 'all_test'").to_string(index=False))

if __name__=='__main__':main()
