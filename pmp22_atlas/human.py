"""Human interval atlas, donor-resolved CAGE and separate rodent state evidence."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
from pathlib import Path

from .coordinates import ChainMap, bed_rows, read_fasta, reverse_complement
from .resources import CACHE, ROOT

DATA=ROOT/'data'; OUT=DATA/'human'; OUT.mkdir(exist_ok=True)


def table(path, rows, fields=None):
    rows=list(rows)
    if not rows and fields is None:raise ValueError(f'Empty table without schema: {path}')
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def mapping():
    forward=ChainMap(CACHE/'hg18ToHg38.over.chain.gz','chr17')
    reverse=ChainMap(CACHE/'hg38ToHg18.over.chain.gz','chr17')
    sequence=read_fasta(CACHE/'hg38.chr17.fa.gz')
    # Coordinates transcribed from Methods, not reconstructed from approximate distances.
    source=[
        ('PMP22_INTRONIC',15091959,15092201,'intronic_enhancer','JONES2011','Methods: +11 kb luciferase construct'),
        ('PMP22_INTRONIC_TRANSGENE',15090965,15092611,'intronic_transgene','JONES2011','Methods: mouse transgenic human insert'),
        ('PMP22_P2_REPORTER',15106498,15106933,'promoter_reporter','JONES2011','Methods: -2 kb P2 promoter construct'),
        ('PMP22_DISTAL_A',15253855,15254150,'distal_enhancer','JONES2012','Methods: region A human reporter'),
        ('PMP22_DISTAL_B',15250013,15250409,'distal_enhancer','JONES2012','Methods: region B human reporter'),
        ('PMP22_DISTAL_C',15221688,15222096,'distal_enhancer','JONES2012','Methods: region C human reporter'),
    ]
    rows=[];audit=[]
    for rid,start,end,kind,source_id,locator in source:
        result=forward.map_interval(start-1,end,1.0)
        if result['status']!='mapped':raise ValueError((rid,result))
        back=reverse.map_interval(result['start'],result['end'],1.0)
        exact=back['status']=='mapped' and back['start']==start-1 and back['end']==end
        if not exact:raise ValueError(f'Failed reciprocal mapping: {rid}')
        rows.append(dict(region_id=rid,assembly='GRCh38',chrom=result['chrom'],
                         start=result['start'],end=result['end'],strand='-' if result['strand']=='+' else '+',
                         element_class=kind,source_id=source_id,evidence_status='measured',
                         orthology_status='human_reference_remapped_reporter_evidence'))
        audit.append(dict(region_id=rid,source_assembly='hg18',source_chrom='chr17',
            published_start=start,published_end=end,interpretation='UCSC displayed 1-based closed; start minus one for BED',
            source_start_0based=start-1,source_end_exclusive=end,
            target_assembly='GRCh38',target_start=result['start'],target_end=result['end'],
            chain_strand=result['strand'],mapped_fraction=result['coverage'],
            reciprocal_exact=exact,chain_id=result['chain_id'],locator=locator,
            sequence_status='reference sequence; cloned insert endpoints require plasmid/supplement confirmation to resolve any 1 bp convention ambiguity'))
    # Promoter-associated first exon families include the longer GENCODE 50 5-prime boundaries.
    gencode=json.loads((CACHE/'gencode50_locus.json').read_text())['wgEncodeGencodeCompV50']
    pmp=[t for t in gencode if t['name2']=='PMP22']
    for rid,exon_start in [('PMP22_P1',15265153),('PMP22_P2',15262428)]:
        family=[t for t in pmp if int(t['exonStarts'].strip(',').split(',')[-1])==exon_start and t['cdsStart']<t['cdsEnd']]
        end=max(t['txEnd'] for t in family)
        rows.append(dict(region_id=rid,assembly='GRCh38',chrom='chr17',start=exon_start,end=end,strand='-',
            element_class='promoter_first_exon',source_id='GENCODE_R50',evidence_status='measured',
            orthology_status='annotated_human_first_exon_family'))
    distal=[r for r in rows if r['region_id'].startswith('PMP22_DISTAL_')]
    rows.append(dict(region_id='PMP22_DISTAL',assembly='GRCh38',chrom='chr17',
        start=min(r['start'] for r in distal),end=max(r['end'] for r in distal),strand='-',
        element_class='distal_cluster_envelope',source_id='JONES2012',evidence_status='proposed',
        orthology_status='envelope_of_human_reporter_intervals_not_a_tested_human_deletion'))
    table(OUT/'mapping_audit.tsv',audit)
    table(DATA/'regions.tsv',rows)
    table(OUT/'gencode50_pmp22_transcripts.tsv',pmp)
    with (OUT/'regions.grch38.bed').open('w') as f, (OUT/'regions.grch38.fa').open('w') as fa:
        for row in rows:
            f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['region_id']}\t0\t{row['strand']}\n")
            seq=sequence[row['start']:row['end']]
            if row['strand']=='-':seq=reverse_complement(seq)
            fa.write(f">{row['region_id']} GRCh38 {row['chrom']}:{row['start']}-{row['end']} strand={row['strand']} reference_only\n")
            fa.write('\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
    return rows


def human_activity(regions):
    peaks=list(bed_rows(CACHE/'ENCFF602YVR.bed.gz'))
    unique={tuple(r[:3]):r for r in peaks}
    peaks17=[r for r in unique.values() if r[0]=='chr17' and int(r[1])<15500000 and int(r[2])>15000000]
    peak_export=[dict(chrom=r[0],start=int(r[1]),end=int(r[2]),name=r[3],signal=float(r[6]),summit_offset=int(r[9]),source='ENCFF602YVR') for r in peaks17]
    table(OUT/'schwann_atac_locus.tsv',peak_export)
    relaxed={tuple(r[:3]):r for r in bed_rows(CACHE/'ENCFF632FMT.bed.gz')
             if r[0]=='chr17' and int(r[1])<15500000 and int(r[2])>15000000}
    table(OUT/'schwann_atac_broad_locus.tsv',[
        dict(chrom=r[0],start=int(r[1]),end=int(r[2]),name=r[3],signal=float(r[6]),source='ENCFF632FMT')
        for r in relaxed.values()])
    cage={};qc=[];cage_locus=[]
    to19=ChainMap(CACHE/'hg38ToHg19.over.chain.gz','chr17')
    to38=ChainMap(CACHE/'hg19ToHg38.over.chain.gz','chr17')
    for donor in ['CNhs12073','CNhs12345','CNhs12621']:
        total=0;sites=0;local=[]
        for fields in bed_rows(CACHE/f'{donor}.hg19.ctss.bed.gz'):
            count=float(fields[4]);total+=count;sites+=1
            if fields[0]=='chr17' and 14900000<int(fields[1])<15500000:
                local.append((int(fields[1]),int(fields[2]),fields[5],count))
                mapped=to38.map_interval(int(fields[1]),int(fields[2]),1)
                if mapped['status']=='mapped':
                    strand=fields[5] if mapped['strand']=='+' else ('-' if fields[5]=='+' else '+')
                    cage_locus.append(dict(donor=donor,chrom=mapped['chrom'],start=mapped['start'],end=mapped['end'],strand=strand,tags=count))
        cage[donor]=(total,local)
        qc.append(dict(donor=donor,total_tag_weight=total,ctss_records=sites,
                       source_assembly='hg19',cell_context='cultured primary Schwann; myelinating state not established'))
    table(OUT/'cage_qc.tsv',qc);table(OUT/'cage_locus.grch38.tsv',cage_locus)
    measurements=[];promoter=[]
    conservation=json.loads((CACHE/'phastCons_locus.json').read_text())['phastCons100way']
    for r in regions:
        overlap=[p for p in peak_export if p['start']<r['end'] and p['end']>r['start']]
        con=[p for p in conservation if p['start']<r['end'] and p['end']>r['start']]
        covered=sum(min(p['end'],r['end'])-max(p['start'],r['start']) for p in con)
        weighted=sum((min(p['end'],r['end'])-max(p['start'],r['start']))*p['value'] for p in con)
        measurements.append(dict(region_id=r['region_id'],schwann_atac_unique_peak_intervals=len(overlap),
            atac_peak_signal_max=max((p['signal'] for p in overlap),default=''),
            conservation_covered_bp=covered,phastcons100way_mean=round(weighted/covered,4) if covered else '',
            inference='observed_accessibility_or_conservation_not_causal_assignment'))
        if r['region_id'] not in {'PMP22_P1','PMP22_P2'}:continue
        # Count in exon-family windows plus 100 bp downstream and a declared upstream flank.
        for flank in [100,250,500]:
            a,b=r['start']-100,r['end']+flank
            mapped=to19.map_interval(a,b,1)
            if mapped['status']!='mapped':raise ValueError(mapped)
            for donor,(total,sites) in cage.items():
                count=count_ctss(sites,mapped['start'],mapped['end'],'-')
                promoter.append(dict(donor=donor,promoter=r['region_id'],upstream_flank_bp=flank,
                    grch38_start=a,grch38_end=b,hg19_start=mapped['start'],hg19_end=mapped['end'],
                    minus_strand_tags=count,tags_per_million=round(count/total*1e6,5),
                    interpretation='promoter-family CAGE initiation; not mature-RNA isoform abundance'))
    for region,measurement in zip(regions,measurements):
        measurement['broad_atac_unique_peak_intervals']=sum(
            int(p[1])<region['end'] and int(p[2])>region['start'] for p in relaxed.values())
        measurement['atac_threshold_status']=('IDR_supported' if measurement['schwann_atac_unique_peak_intervals']
            else 'broad_only' if measurement['broad_atac_unique_peak_intervals'] else 'not_detected_in_released_peaks')
    table(OUT/'region_activity.tsv',measurements);table(OUT/'promoter_cage.tsv',promoter)
    return measurements,qc,promoter


def count_ctss(sites,start,end,strand):
    """Sum tag weights in a half-open interval on the requested genomic strand."""
    return sum(n for x,y,s,n in sites if s==strand and x<end and y>start)


def source_provenance():
    rows=[]
    for acc in ['ENCSR726QTF','ENCSR205TUH']:
        data=json.loads((CACHE/(acc+'.json')).read_text())
        for rep in data['replicates']:
            sample=rep['library']['biosample'];donor=sample['donor']
            rows.append(dict(experiment=acc,biosample=sample['accession'],donor=donor['accession'],
                age=donor.get('age'),age_units=donor.get('age_units'),sex=donor.get('sex'),
                biological_replicate=rep['biological_replicate_number'],technical_replicate=rep['technical_replicate_number'],
                pooled_product='ENCSR301DUJ / ENCFF439MKP',
                note='two donor source experiments; released Schwann peak product is pooled and cannot estimate donor variance'))
    table(OUT/'atac_provenance.tsv',rows)
    return rows


def injury(regions):
    rows=[];summary=[];rejected=[]
    for name in ['InjuryDB_peaks','ShamDB_peaks','peaks_Injury','peaks_Sham']:
        path=CACHE/f'GSE63103_H3K27ac_{name}.bed.gz';n=0;local=[]
        for a in bed_rows(path):
            n+=1
            if any('e' in v.lower() for v in a[1:3]):
                rejected.append(dict(track=name,chrom=a[0],original_start=a[1],original_end=a[2],reason='scientific_notation_coordinate_precision_ambiguous'))
                continue
            if a[0]=='chr10' and int(a[1])<50346994 and int(a[2])>48316999:
                row=dict(track=name,assembly='rn5',chrom=a[0],start=int(a[1]),end=int(a[2]),value=a[3] if len(a)>3 else '')
                rows.append(row);local.append(row)
        summary.append(dict(track=name,genome_records=n,locus_records=len(local),
            rejected_coordinates=sum(r['track']==name for r in rejected),
            reported_test='original authors DBChIP' if 'DB' in name else 'original authors peak calls',
            new_test='none; reuse of submitted interval labels'))
    table(OUT/'rat_injury_locus.rn5.tsv',rows);table(OUT/'rat_injury_summary.tsv',summary)
    table(OUT/'rat_injury_rejected_coordinates.tsv',rejected)
    # Report orthology as an annotation, never as a human injury measurement.
    lift=ChainMap(CACHE/'rn5ToHg38.over.chain.gz','chr10');mapped=[]
    for row in rows:
        result=lift.map_interval(row['start'],row['end'],0.5)
        mapped.append(dict(source_track=row['track'],rn5_start=row['start'],rn5_end=row['end'],
            status=result['status'],grch38_chrom=result.get('chrom',''),grch38_start=result.get('start',''),
            grch38_end=result.get('end',''),mapped_fraction=result.get('coverage',''),
            strand=result.get('strand',''),aligned_blocks=json.dumps(result.get('aligned_blocks',[])),
            interpretation='orthology annotation; tissue response remains rat evidence'))
    table(OUT/'rat_injury_orthology.tsv',mapped)
    # Preserve lower-coverage candidate evidence rather than hiding it in a zero count.
    audits=[]
    for row in rows:
        accepted=next(p for p in mapped if p['source_track']==row['track']
                      and p['rn5_start']==row['start'] and p['rn5_end']==row['end'])
        for candidate in lift.candidates(row['start'],row['end']):
            for region in regions:
                if candidate['chrom']!=region['chrom']:continue
                bp=sum(max(0,min(b,region['end'])-max(a,region['start']))
                       for a,b in candidate['aligned_blocks'])
                if not bp:continue
                use=accepted['status']=='mapped' and candidate['coverage']>=.5
                audits.append(dict(region_id=region['region_id'],rat_track=row['track'],
                    source_interval_id=f"rn5:chr10:{row['start']}-{row['end']}",chain_id=candidate['chain_id'],
                    mapped_fraction=candidate['coverage'],aligned_overlap_bp=bp,
                    accepted=use,exclusion_reason='' if use else 'below_50pct_source_span' if candidate['coverage']<.5 else accepted['status']))
    table(OUT/'rat_mapping_overlap_audit.tsv',audits)
    overlaps=[]
    for region in regions:
        for track in ['InjuryDB_peaks','ShamDB_peaks','peaks_Injury','peaks_Sham']:
            hits=[p for p in mapped if p['source_track']==track and p['status']=='mapped'
                  and p['grch38_chrom']==region['chrom'] and any(
                      a<region['end'] and b>region['start'] for a,b in json.loads(p['aligned_blocks']))]
            filtered={p['source_interval_id'] for p in audits if p['region_id']==region['region_id']
                      and p['rat_track']==track and not p['accepted']}
            overlaps.append(dict(region_id=region['region_id'],rat_track=track,
                orthologous_peak_count=len(hits),minimum_mapped_fraction=min((p['mapped_fraction'] for p in hits),default=''),
                filtered_candidate_peak_count=len(filtered),
                evaluation_status='partial_mapping_filtered' if filtered else 'accepted_overlap_present' if hits else 'no_accepted_overlap_not_a_functional_negative',
                interpretation='Counts conditional on source-span mapping and peak calling; not independent experiments, a state effect, or complete activity coverage'))
    table(OUT/'region_rat_injury_support.tsv',overlaps)
    return summary


def mouse_deletion():
    # Pantera 2020 Fig. 1 states mm10 chr11:63001867-63042405.
    # Preserve the display convention uncertainty, as for historical human constructs.
    result=ChainMap(CACHE/'mm10ToHg38.over.chain.gz','chr11').map_interval(63001866,63042405,.5)
    table(OUT/'mouse_se_orthology.tsv',[dict(source='PANTERA2020',source_assembly='mm10',
        published_interval='chr11:63001867-63042405',source_start_0based=63001866,
        source_end_exclusive=63042405,source_convention='displayed interval interpreted 1-based closed; endpoints need junction confirmation',
        status=result['status'],grch38_chrom=result.get('chrom',''),grch38_start=result.get('start',''),
        grch38_end=result.get('end',''),mapped_fraction=result.get('coverage',''),strand=result.get('strand',''),
        aligned_blocks=json.dumps(result.get('aligned_blocks',[])),
        interpretation='mapped envelope contains alignment gaps; mouse deletion effect is not a human interval effect')])


def model_exposure(regions):
    spans=list(bed_rows(CACHE/'borzoi_sequences_human.bed.gz'))
    rows=[]
    for region in regions:
        hits=[p for p in spans if p[0]==region['chrom'] and int(p[1])<region['end'] and int(p[2])>region['start']]
        folds=sorted({p[3] for p in hits})
        center=(region['start']+region['end'])//2
        context=[p for p in spans if p[0]==region['chrom'] and int(p[1])<center+262144 and int(p[2])>center-262144]
        rows.append(dict(region_id=region['region_id'],borzoi_folds=';'.join(folds),
                         centered_input_bp=524288,input_context_start=center-262144,input_context_end=center+262144,
                         context_folds=';'.join(sorted({p[3] for p in context})),
                         overlapping_windows=len(hits),standard_checkpoint_role='training' if folds and not set(folds)&{'fold3','fold4'} else 'review',
                         eligibility='ineligible_for_independent_PMP22_genomic_benchmark'))
    table(OUT/'borzoi_locus_exposure.tsv',rows)


def main():
    regions=mapping();activity,qc,promoter=human_activity(regions)
    provenance=source_provenance();state=injury(regions);mouse_deletion();model_exposure(regions)
    summary=dict(regions=len(regions),source_donors=len({r['donor'] for r in provenance}),
                 human_cage_donors=len(qc),activity=activity,promoter_cage=promoter,rat_injury=state)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
