"""Regression checks for lost measurement units, shared evidence and shifted populations."""
import unittest
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score,roc_auc_score

ROOT=Path(__file__).resolve().parents[1]

class ReviewCorrections(unittest.TestCase):
    def test_low_coverage_overlap_is_audited_without_acceptance(self):
        audit=pd.read_csv(ROOT/'data/human/rat_mapping_overlap_audit.tsv',sep='\t')
        p=audit[(audit.source_interval_id=='rn5:chr10:49331700-49331999') &
                (audit.region_id=='PMP22_INTRONIC_TRANSGENE')].iloc[0]
        self.assertEqual(138,p.aligned_overlap_bp)
        self.assertFalse(p.accepted)
        self.assertAlmostEqual(138/299,p.mapped_fraction)
        support=pd.read_csv(ROOT/'data/human/region_rat_injury_support.tsv',sep='\t')
        r=support[(support.region_id==p.region_id)&(support.rat_track==p.rat_track)].iloc[0]
        self.assertGreater(r.filtered_candidate_peak_count,0)
        self.assertEqual('partial_mapping_filtered',r.evaluation_status)

    def test_shared_source_peak_retains_one_identity_across_regions(self):
        matrix=pd.read_csv(ROOT/'data/release/evidence_matrix.tsv',sep='\t').fillna('')
        identity='ENCFF703WKZ:chr17:15247227-15273130'
        rows=matrix[matrix.source_interval_ids.str.split(';').apply(lambda ids:identity in ids)]
        self.assertEqual({'PMP22_INTRONIC','PMP22_INTRONIC_TRANSGENE','PMP22_P2_REPORTER','PMP22_P1','PMP22_P2'},set(rows.region_id))
        self.assertEqual(1,rows.evidence_family_id.nunique())
        self.assertEqual(1,rows.donor_ids.nunique())
        self.assertTrue(rows.donor_ids.str.len().gt(0).all())

    def test_stratified_scores_preserve_actual_prevalence(self):
        pred=pd.read_csv(ROOT/'benchmark/results/predictions.tsv',sep='\t',float_precision='round_trip')
        part=pred[pred.in_gc_matched_subset & pred.distance_to_tss.gt(2000)]
        scores=pd.read_csv(ROOT/'benchmark/results/stratified_metrics.tsv',sep='\t')
        scores=scores[scores.subset.eq('GC_matched_beyond_2kb_TSS')].set_index('model')
        self.assertEqual(139,int(part.label.sum()))
        self.assertEqual(386,len(part))
        self.assertAlmostEqual(part.label.mean(),scores.loc['prevalence','average_precision'])
        for name in ['4mer_sequence_logistic','GC_CpG','distance_to_annotated_TSS']:
            self.assertAlmostEqual(average_precision_score(part.label,part[name]),scores.loc[name,'average_precision'])
            self.assertAlmostEqual(roc_auc_score(part.label,part[name]),scores.loc[name,'AUROC'])

    def test_accessible_controls_have_distinct_real_peak_support(self):
        controls=pd.read_csv(ROOT/'data/human/candidate_accessible_controls.tsv',sep='\t')
        peaks=pd.read_csv(ROOT/'data/human/accessible_control_source_peaks.tsv',sep='\t').set_index('source_interval_id')
        self.assertEqual(2,len(controls))
        self.assertEqual(2,controls.source_interval_id.nunique())
        for r in controls.itertuples():
            self.assertIn(r.source_interval_id,peaks.index)
            source=peaks.loc[r.source_interval_id]
            self.assertEqual(source.summit,r.start+204)
            self.assertLess(r.start,source.end)
            self.assertGreater(r.end,source.start)
            self.assertEqual(409,r.end-r.start)
            self.assertLessEqual(r.gc_difference,.03)
            self.assertGreater(r.distance_to_annotated_TSS,2000)
            self.assertFalse(r.start<16265360 and r.end>14229778)

if __name__=='__main__':unittest.main()
