"""Post-test diagnostics declared in benchmark/amendments.md; no model re-selection."""
from concurrent.futures import ProcessPoolExecutor
import gzip
import json

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score

from .benchmark import BENCH,RESULTS,PROTOCOL,bootstrap,metrics,feature_sequence,motifs
from .coordinates import read_fasta
from .resources import CACHE,sha256
from .human import table


def conservation_chrom(chrom):
    frame=pd.read_csv(BENCH/'splits.tsv',sep='\t')
    windows=frame[frame.chrom.eq(chrom)].sort_values('start')
    intervals=list(windows[['example_id','start','end']].itertuples(index=False,name=None))
    result={name:[0.,0] for name,_,_ in intervals};i=0;pos=0
    with gzip.open(CACHE/(chrom+'.phastCons100way.wigFix.gz'),'rt') as handle:
        for line in handle:
            if line.startswith('fixedStep'):
                header=dict(x.split('=',1) for x in line.split()[1:])
                if header.get('step')!='1' or header.get('span','1')!='1':raise ValueError(header)
                if header['chrom']!=chrom:raise ValueError(header)
                pos=int(header['start'])-1
                continue
            while i<len(intervals) and intervals[i][2]<=pos:i+=1
            if i==len(intervals):break
            name,start,end=intervals[i]
            if start<=pos<end:
                result[name][0]+=float(line);result[name][1]+=1
            pos+=1
    print('Conservation extracted',chrom,len(result),flush=True)
    return [dict(example_id=name,phastcons_mean=total/n if n else 0.,
                 covered_fraction=n/600,covered_bases=n) for name,(total,n) in result.items()]


def matched_model_contrast(predictions):
    """The strongest non-sequence comparator matters, not only the GC baseline."""
    part=predictions[predictions.in_gc_matched_subset].reset_index(drop=True)
    blocks=part.block.unique();byblock={b:np.where(part.block.eq(b))[0] for b in blocks}
    rng=np.random.default_rng(PROTOCOL['seed']+1);deltas=[]
    for _ in range(500):
        ix=np.concatenate([byblock[b] for b in rng.choice(blocks,len(blocks),replace=True)])
        p=part.iloc[ix]
        if p.label.nunique()<2:continue
        deltas.append(average_precision_score(p.label,p['4mer_sequence_logistic'])-
                      average_precision_score(p.label,p['distance_to_annotated_TSS']))
    table(RESULTS/'matched_sequence_vs_distance.tsv',[dict(subset='GC_matched_test',
        contrast='sequence_composite_AP_minus_distance_AP',blocks=len(blocks),draws=len(deltas),
        AP_difference=average_precision_score(part.label,part['4mer_sequence_logistic'])-
                      average_precision_score(part.label,part['distance_to_annotated_TSS']),
        low=float(np.quantile(deltas,.025)),high=float(np.quantile(deltas,.975)),
        interpretation='post-test diagnostic; no new model selection')])


def main():
    frame=pd.read_csv(BENCH/'splits.tsv',sep='\t')
    path=RESULTS/'conservation_features.tsv'
    with ProcessPoolExecutor(max_workers=3) as pool:
        rows=[row for group in pool.map(conservation_chrom,frame.chrom.unique()) for row in group]
    table(path,rows)
    features=frame[['example_id']].merge(pd.DataFrame(rows),on='example_id',validate='one_to_one')
    x=features[['phastcons_mean','covered_fraction']].to_numpy()
    train=frame.split.eq('train');val=frame.split.eq('validation');test=frame.split.eq('test');y=frame.label.to_numpy()
    best=None;fitlog=[]
    for c in PROTOCOL['regularization_C_grid']:
        model=make_pipeline(StandardScaler(),LogisticRegression(C=c,max_iter=2000,random_state=PROTOCOL['seed']))
        model.fit(x[train],y[train]);score=average_precision_score(y[val],model.predict_proba(x[val])[:,1])
        fitlog.append(dict(model='conservation',C=c,validation_AP=score))
        if best is None or score>best[0]:best=(score,c,model)
    predictions=pd.read_csv(RESULTS/'predictions.tsv',sep='\t',float_precision='round_trip')
    test_ids=frame.loc[test,'example_id'].tolist()
    if test_ids!=predictions.example_id.tolist():raise ValueError('Test order mismatch')
    p=best[2].predict_proba(x[test])[:,1]
    predictions['conservation']=p
    matched_model_contrast(predictions)
    names=['prevalence','GC_CpG','distance_to_annotated_TSS','SOX10_EGR2_TEAD1_motifs',
           '4mer_sequence_logistic','permuted_training_labels','conservation']
    preds={name:predictions[name].to_numpy() for name in names}
    subset=predictions.in_gc_matched_subset.to_numpy(dtype=bool)
    metricrows=[]
    for name,mask in [('all_test',np.ones(len(p),dtype=bool)),('GC_matched_test',subset)]:
        metricrows.append(dict(model='conservation',subset=name,n=int(mask.sum()),**metrics(y[test][mask],p[mask])))
    matched=predictions.loc[subset].reset_index(drop=True)
    table(RESULTS/'GC_matched_bootstrap.tsv',bootstrap(matched,{k:v[subset] for k,v in preds.items()}))
    table(RESULTS/'supplemental_bootstrap.tsv',bootstrap(predictions,{'conservation':p,'GC_CpG':preds['GC_CpG']}))
    table(RESULTS/'conservation_validation.tsv',fitlog)
    predictions[['example_id','conservation']].to_csv(RESULTS/'conservation_predictions.tsv',sep='\t',index=False)
    scaler,lr=best[2].steps[0][1],best[2].steps[1][1]
    (RESULTS/'conservation_model.json').write_text(json.dumps(dict(C=best[1],validation_AP=best[0],
        feature_names=['phastcons_mean','covered_fraction'],mean=scaler.mean_.tolist(),scale=scaler.scale_.tolist(),
        coefficients=lr.coef_[0].tolist(),intercept=float(lr.intercept_[0]),features_sha256=sha256(path)),indent=2)+'\n')
    params=json.loads((RESULTS/'model_parameters.json').read_text())['4mer_sequence_logistic'];pwms=motifs()
    shiftrows=[]
    for chrom in PROTOCOL['test_chromosomes']:
        sequence=read_fasta(CACHE/f'hg38.{chrom}.fa.gz')
        for shift in [-100,100]:
            for row in predictions[predictions.chrom.eq(chrom)].itertuples(index=False):
                s=sequence[row.start+shift:row.end+shift]
                if len(s)!=600 or any(b not in 'ACGT' for b in s):
                    shiftrows.append(dict(example_id=row.example_id,shift_bp=shift,label=row.label,prediction='',status='ambiguous_or_out_of_bounds'))
                    continue
                vector=np.r_[0.,feature_sequence(s,pwms)][params['columns']]
                z=(vector-np.array(params['mean']))/np.array(params['scale'])
                probability=float(expit(z@np.array(params['coefficients'])+params['intercept']))
                shiftrows.append(dict(example_id=row.example_id,shift_bp=shift,label=row.label,prediction=probability,status='scored'))
    shifts=pd.DataFrame(shiftrows);table(RESULTS/'shift_predictions.tsv',shiftrows)
    for shift in [-100,100]:
        part=shifts[shifts.shift_bp.eq(shift)&shifts.status.eq('scored')]
        metricrows.append(dict(model='4mer_sequence_logistic',subset=f'shift_{shift}_bp',n=len(part),
            **metrics(part.label.to_numpy(),part.prediction.to_numpy(dtype=float))))
    table(RESULTS/'supplemental_metrics.tsv',metricrows)
    print(pd.DataFrame(metricrows).to_string(index=False))

if __name__=='__main__':main()
