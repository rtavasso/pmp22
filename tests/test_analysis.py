"""Tests for errors that could change scientific interpretation, not just execution."""
import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score,roc_auc_score,brier_score_loss

from pmp22_atlas.coordinates import ChainMap,IntervalIndex,reverse_complement
from pmp22_atlas.human import count_ctss
from pmp22_atlas.benchmark import feature_sequence
from pmp22_atlas import resources

ROOT=Path(__file__).resolve().parents[1]

class CoordinateTests(unittest.TestCase):
    def chain(self,text):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        path=Path(tmp.name)/'test.chain.gz'
        with gzip.open(path,'wt') as f:f.write(text)
        return ChainMap(path,'chrS')

    def test_forward_and_reverse_chain_orientation(self):
        for strand,expected in [('+',(510,530)),('-',(1470,1490))]:
            chain=self.chain(f'chain 100 chrS 1000 + 100 200 chrT 2000 {strand} 500 600 1\n100\n')
            hit=chain.map_interval(110,130,1)
            self.assertEqual((hit['start'],hit['end']),expected)
            self.assertEqual(hit['strand'],strand)
            self.assertEqual(hit['coverage'],1)

    def test_gaps_are_not_mapped_bases(self):
        chain=self.chain('chain 100 chrS 1000 + 100 200 chrT 2000 + 500 605 1\n40 5 10\n55\n')
        self.assertEqual(chain.map_interval(135,155,.95)['status'],'unmapped')
        hit=chain.map_interval(135,155,.75)
        self.assertEqual(hit['coverage'],.75)
        self.assertEqual(hit['aligned_blocks'],[(535,540),(550,560)])
        self.assertFalse(IntervalIndex(hit['aligned_blocks']).overlaps(542,548))

    def test_multimapping_is_not_silently_ranked_by_score(self):
        chain=self.chain('chain 100 chrS 1000 + 100 200 chrT 2000 + 500 600 1\n100\n\n'
                         'chain 50 chrS 1000 + 100 200 chrT 2000 + 700 800 2\n100\n')
        self.assertEqual(chain.map_interval(120,130)['status'],'ambiguous')

    def test_duplicate_equivalent_chains_do_not_create_false_ambiguity(self):
        chain=self.chain('chain 100 chrS 1000 + 100 200 chrT 2000 + 500 600 1\n100\n\n'
                         'chain 50 chrS 1000 + 100 200 chrT 2000 + 500 600 2\n100\n')
        self.assertEqual(chain.map_interval(120,130)['status'],'mapped')

    def test_half_open_touching_empty_and_invalid_intervals(self):
        index=IntervalIndex([(10,20),(20,25)])
        self.assertFalse(index.overlaps(25,30));self.assertFalse(index.overlaps(5,10))
        self.assertFalse(index.overlaps(15,15));self.assertTrue(index.overlaps(24,26))
        chain=self.chain('chain 100 chrS 1000 + 100 200 chrT 2000 + 500 600 1\n100\n')
        with self.assertRaises(ValueError):chain.map_interval(100,100)
        with self.assertRaises(ValueError):chain.map_interval(-1,10)

    def test_ctss_counts_strand_and_boundaries(self):
        sites=[(9,10,'-',2),(10,11,'-',3),(19,20,'-',5),(20,21,'-',7),(12,13,'+',100)]
        self.assertEqual(count_ctss(sites,10,20,'-'),8)
        self.assertEqual(count_ctss(sites,10,20,'+'),100)

    def test_sequence_features_are_reverse_complement_invariant(self):
        seq=('AACGGTCTAGCGATATCGGGCCAATCTT'*24)[:600]
        pwm=np.array([[1.,2.,-1.,-2.],[.2,-.3,1.,.4],[-1.,.5,.8,.2]])
        np.testing.assert_allclose(feature_sequence(seq,[pwm]),feature_sequence(reverse_complement(seq),[pwm]),atol=1e-12)


class ReleaseIntegrityTests(unittest.TestCase):
    def test_test_examples_have_no_overlap_or_sequence_duplicates(self):
        frame=pd.read_csv(ROOT/'benchmark/splits.tsv',sep='\t')
        self.assertFalse(frame.sequence_sha256.duplicated().any())
        expected={'train':{'chr1','chr2'},'validation':{'chr3'},'test':{'chr17','chr22'}}
        for split,chroms in expected.items():self.assertEqual(set(frame[frame.split.eq(split)].chrom),chroms)
        for _,part in frame.groupby('chrom'):
            part=part.sort_values('start')
            self.assertTrue((part.start.to_numpy()[1:]>=part.end.to_numpy()[:-1]).all())
            self.assertEqual(part.label.value_counts().loc[0],part.label.value_counts().loc[1])
        lock=json.loads((ROOT/'benchmark/results/dataset_lock.json').read_text())
        for key,name in [('split_sha256','splits.tsv'),('protocol_sha256','protocol.json')]:
            self.assertEqual(lock[key],resources.sha256(ROOT/'benchmark'/name))

    def test_reported_scores_recompute_from_exported_predictions(self):
        p=pd.read_csv(ROOT/'benchmark/results/predictions.tsv',sep='\t',float_precision='round_trip')
        metrics=pd.read_csv(ROOT/'benchmark/results/metrics.tsv',sep='\t',float_precision='round_trip')
        for row in metrics[metrics.subset.isin(['all_test','GC_matched_test'])].itertuples(index=False):
            q=p if row.subset=='all_test' else p[p.in_gc_matched_subset]
            self.assertEqual(row.n,len(q))
            self.assertAlmostEqual(row.average_precision,average_precision_score(q.label,q[row.model]),places=12)
            self.assertAlmostEqual(row.AUROC,roc_auc_score(q.label,q[row.model]),places=12)
            self.assertAlmostEqual(row.Brier_score,brier_score_loss(q.label,q[row.model]),places=12)

    def test_model_selection_uses_validation_maximum(self):
        val=pd.read_csv(ROOT/'benchmark/results/validation.tsv',sep='\t')
        chosen=pd.read_csv(ROOT/'benchmark/results/selected_models.tsv',sep='\t')
        for row in chosen.itertuples(index=False):
            part=val[val.model.eq(row.model)]
            self.assertAlmostEqual(row.validation_AP,part.validation_AP.max(),places=12)
            self.assertEqual(row.C,part.loc[part.validation_AP.idxmax(),'C'])

    def test_effects_keep_species_units_and_missing_uncertainty(self):
        effects=pd.read_csv(ROOT/'data/assay_effects.tsv',sep='\t',keep_default_na=False)
        self.assertTrue(effects.species.astype(bool).all())
        self.assertTrue(effects.locator.astype(bool).all())
        self.assertTrue(effects.control.astype(bool).all())
        self.assertFalse((effects.value_origin=='digitized').any())
        qualitative=effects.value_origin.eq('author_qualitative')
        self.assertTrue(effects.loc[qualitative,'effect_value'].eq('').all())
        native=effects.evidence_layer.isin(['native_RNA','protein','function'])
        self.assertTrue(set(effects.loc[native,'species'])<= {'rat','mouse'})

    def test_human_reference_sequences_have_correct_lengths_and_orientation_labels(self):
        regions=pd.read_csv(ROOT/'data/regions.tsv',sep='\t').set_index('region_id')
        text=(ROOT/'data/human/regions.grch38.fa').read_text()
        for record in text.split('>')[1:]:
            lines=record.splitlines();rid=lines[0].split()[0];seq=''.join(lines[1:])
            self.assertEqual(len(seq),regions.loc[rid,'end']-regions.loc[rid,'start'])
            self.assertIn('strand=-',lines[0]);self.assertIn('reference_only',lines[0])
        audit=pd.read_csv(ROOT/'data/human/mapping_audit.tsv',sep='\t')
        self.assertTrue(audit.reciprocal_exact.all());self.assertTrue(audit.mapped_fraction.eq(1).all())

    def test_proposed_controls_are_distinct_and_distant_from_known_intervals(self):
        controls=pd.read_csv(ROOT/'data/human/candidate_negative_controls.tsv',sep='\t')
        known=pd.read_csv(ROOT/'data/regions.tsv',sep='\t')
        known=known[~known.region_id.eq('PMP22_DISTAL')]
        self.assertEqual(len(controls),2);self.assertTrue(controls.gc_difference.le(.03).all())
        for r in controls.itertuples(index=False):
            self.assertEqual(r.end-r.start,409)
            self.assertFalse(((known.start-10000<r.end)&(known.end+10000>r.start)).any())
            self.assertIn('untested',r.status)

    def test_frozen_public_metadata_snapshots_match_acquisition_hashes(self):
        lock=json.loads((ROOT/'data/resource_lock.json').read_text())
        for path in (ROOT/'data/snapshots').glob('*.gz'):
            original=path.name[:-3]
            self.assertIn(original,lock)
            self.assertEqual(hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest(),lock[original]['sha256'])


class ResourceLockTests(unittest.TestCase):
    def fixture(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup)
        root=Path(tmp.name);cache=root/'data/cache';cache.mkdir(parents=True)
        lock=root/'data/resource_lock.json';payload=b'{"frozen":true}\n'
        entry=dict(url='https://example.org/public.json',sha256=hashlib.sha256(payload).hexdigest())
        lock.write_text(json.dumps({'public.json':entry}))
        return root,cache,lock,payload

    def test_fresh_checkout_restores_metadata_without_network(self):
        root,cache,lock,payload=self.fixture();snap=root/'data/snapshots';snap.mkdir()
        (snap/'public.json.gz').write_bytes(gzip.compress(payload))
        with patch.object(resources,'ROOT',root),patch.object(resources,'CACHE',cache),patch.object(resources,'LOCK',lock),patch.object(resources.requests,'Session') as network:
            path=resources.fetch('public.json','https://example.org/public.json',source='test')
            self.assertEqual(path.read_bytes(),payload);network.assert_not_called()

    def test_corrupt_cache_and_changed_url_are_rejected(self):
        root,cache,lock,payload=self.fixture();path=cache/'public.json'
        with patch.object(resources,'ROOT',root),patch.object(resources,'CACHE',cache),patch.object(resources,'LOCK',lock):
            path.write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError,'Cached bytes changed'):
                resources.fetch('public.json','https://example.org/public.json',source='test')
            path.write_bytes(payload)
            with self.assertRaisesRegex(ValueError,'Source URL changed'):
                resources.fetch('public.json','https://example.org/other.json',source='test')


if __name__=='__main__':unittest.main()
