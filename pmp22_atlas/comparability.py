"""Post-review test-population diagnostics; frozen models and labels are unchanged."""
import hashlib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from .benchmark import RESULTS,metrics
from .coordinates import read_fasta,reverse_complement
from .human import table
from .resources import CACHE

MODELS=['prevalence','GC_CpG','distance_to_annotated_TSS','SOX10_EGR2_TEAD1_motifs',
        '4mer_sequence_logistic','permuted_training_labels']

def main():
    frame=pd.read_csv(RESULTS/'predictions.tsv',sep='\t',float_precision='round_trip')
    frame['cpg_fraction']=np.nan
    for chrom,group in frame.groupby('chrom'):
        sequence=read_fasta(CACHE/f'hg38.{chrom}.fa.gz')
        for i,row in group.iterrows():
            s=sequence[row.start:row.end]
            digest=hashlib.sha256(min(s,reverse_complement(s)).encode()).hexdigest()
            if digest!=row.sequence_sha256:raise ValueError(f'Reference sequence mismatch: {row.example_id}')
            frame.loc[i,'cpg_fraction']=s.count('CG')/(len(s)-1)
    frame['within_2kb_TSS']=frame.distance_to_tss.le(2000)
    composition=[];score_rows=[];strata=[]
    for subset,mask in [('all_test',np.ones(len(frame),bool)),('GC_matched_test',frame.in_gc_matched_subset)]:
        part=frame.loc[mask]
        for label,g in part.groupby('label'):
            composition.append(dict(subset=subset,label=label,n=len(g),fraction_of_original_class=len(g)/int(frame.label.eq(label).sum()),
                within_2kb_TSS=int(g.within_2kb_TSS.sum()),within_2kb_TSS_fraction=float(g.within_2kb_TSS.mean()),
                median_TSS_distance=float(g.distance_to_tss.median()),median_GC=float(g.gc_fraction.median()),
                median_CpG=float(g.cpg_fraction.median())))
        for name,col in [('GC_fraction_only','gc_fraction'),('CpG_fraction_only','cpg_fraction')]:
            score_rows.append(dict(subset=subset,ranking_feature=name,n=len(part),positive=int(part.label.sum()),
                prevalence=float(part.label.mean()),average_precision=float(average_precision_score(part.label,part[col])),
                interpretation='exploratory feature ranking, not a fitted biological model'))
    subsets=[('all_test',np.ones(len(frame),bool)),('GC_matched_test',frame.in_gc_matched_subset),
             ('within_2kb_TSS',frame.within_2kb_TSS),('beyond_2kb_TSS',~frame.within_2kb_TSS),
             ('GC_matched_beyond_2kb_TSS',frame.in_gc_matched_subset & ~frame.within_2kb_TSS)]
    for subset,mask in subsets:
        part=frame.loc[mask]
        for model in MODELS:
            strata.append(dict(subset=subset,model=model,n=len(part),positive=int(part.label.sum()),
                prevalence=float(part.label.mean()),**metrics(part.label.to_numpy(),part[model].to_numpy()),
                interpretation='post-review exploratory stratum; AP across strata changes population/prevalence'))
    table(RESULTS/'test_composition.tsv',composition)
    table(RESULTS/'composition_only_scores.tsv',score_rows)
    table(RESULTS/'stratified_metrics.tsv',strata)
    frame[['example_id','label','gc_fraction','cpg_fraction','distance_to_tss','within_2kb_TSS','in_gc_matched_subset']].to_csv(
        RESULTS/'comparability_examples.tsv',sep='\t',index=False)
    print(pd.DataFrame(composition).to_string(index=False))
    print(pd.DataFrame(strata).query("subset == 'GC_matched_beyond_2kb_TSS'").to_string(index=False))

if __name__=='__main__':main()
