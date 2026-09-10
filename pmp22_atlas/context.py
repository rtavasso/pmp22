"""Bulk human nerve context, source QC, transcript sensitivity and control intervals."""
import csv
import json
from collections import defaultdict

from .coordinates import bed_rows,read_fasta,IntervalIndex
from .human import OUT,DATA,table,count_ctss
from .resources import CACHE


def read(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    regions=read(DATA/'regions.tsv')
    genes=json.loads((CACHE/'gencode50_locus.json').read_text())['wgEncodeGencodeCompV50']
    tx_to_symbol={r['name'].split('.')[0]:r['name2'] for r in genes}
    grouped=defaultdict(list)
    for row in genes:grouped[row['name2']].append(row)
    table(OUT/'neighbor_genes.tsv',[dict(symbol=name,chrom='chr17',start=min(t['txStart'] for t in ts),
        end=max(t['txEnd'] for t in ts),strand=ts[0]['strand'],transcripts=len(ts),
        note='GENCODE50 locus annotation; not an enhancer-to-gene link') for name,ts in sorted(grouped.items())])
    provenance=[];rna=[];locus=[];overlap=[]
    for selected in read(OUT/'bulk_selected_files.tsv'):
        acc=selected['experiment'];d=json.loads((CACHE/(acc+'.json')).read_text())
        target=d.get('target',{}).get('label','')
        for rep in d['replicates']:
            library=rep['library'];bio=library['biosample'];donor=bio.get('donor',{})
            provenance.append(dict(experiment=acc,assay=d['assay_title'],target=target,
                library=library['accession'],biosample=bio['accession'],donor=donor.get('accession',''),
                age=donor.get('age',''),sex=donor.get('sex',''),bio_replicate=rep['biological_replicate_number'],
                controls=';'.join(c['accession'] if isinstance(c,dict) else c for c in d.get('possible_controls',[])),
                selected_file=selected['file'],annotation=selected['genome_annotation'],
                limitation='bulk tissue; one biological sample per experiment; shared donors across assays'))
        if selected['output_type']=='gene quantifications':
            for row in read(CACHE/selected['filename']):
                symbols={tx_to_symbol[t.split('.')[0]] for t in row['transcript_id(s)'].split(',') if t.split('.')[0] in tx_to_symbol}
                if len(symbols)!=1:continue
                symbol=next(iter(symbols))
                rna.append(dict(experiment=acc,donor=provenance[-1]['donor'],file=selected['file'],symbol=symbol,
                    gene_id=row['gene_id'],annotation=selected['genome_annotation'],TPM=float(row['TPM']),
                    primary_PMP22_gene=row['gene_id'].split('.')[0]=='ENSG00000109099',
                    symbol_assignment='GENCODE50 transcript-ID match; source V29 gene IDs retained separately',
                    expected_count=float(row['expected_count']),source_posterior_TPM_low=row.get('TPM_ci_lower_bound',''),
                    source_posterior_TPM_high=row.get('TPM_ci_upper_bound',''),
                    interpretation='bulk tissue expression; RSEM posterior interval is not between-donor biological uncertainty'))
        else:
            peaks={tuple(r[:3]):r for r in bed_rows(CACHE/selected['filename']) if r[0]=='chr17' and int(r[1])<15500000 and int(r[2])>15000000}
            for p in peaks.values():
                locus.append(dict(experiment=acc,file=selected['file'],target=target or d['assay_title'],chrom=p[0],
                    start=int(p[1]),end=int(p[2]),output_type=selected['output_type']))
            for r in regions:
                hits=sum(int(p[1])<int(r['end']) and int(p[2])>int(r['start']) for p in peaks.values())
                overlap.append(dict(region_id=r['region_id'],experiment=acc,target=target or d['assay_title'],
                    unique_peak_intervals=hits,file=selected['file'],context='bulk tibial nerve; not purified Schwann'))
    table(OUT/'bulk_provenance.tsv',provenance);table(OUT/'bulk_locus_rna.tsv',rna)
    table(OUT/'bulk_locus_peaks.tsv',locus);table(OUT/'region_bulk_support.tsv',overlap)
    qc=[];audit=[]
    for acc in ['ENCFF029IGD','ENCFF540FVR','ENCFF439MKP','ENCFF602YVR','ENCFF632FMT']:
        d=json.loads((CACHE/(acc+'.json')).read_text())
        for item in d.get('quality_metrics',[]):
            if not isinstance(item,dict):continue
            for key,value in item.items():
                if isinstance(value,(int,float)):
                    qc.append(dict(file=acc,metric_type=item['@type'][0],metric=key,value=value,
                        scope='source bulk snATAC library; not Schwann-subset QC or biological replication'))
        for level,items in d.get('audit',{}).items():
            for item in items:
                audit.append(dict(file=acc,level=level,category=item['category'],detail=item['detail'],
                    handling='matching peak MD5 is ENCSR221DFU leg Schwann annotation: same bytes, not independent; internal-status flags retained'))
    table(OUT/'atac_source_qc.tsv',qc);table(OUT/'encode_audit.tsv',audit)
    # Compare promoter-family annotation with reference-transcript first-exon boundaries.
    ref=json.loads((CACHE/'refseq_locus.json').read_text())['ncbiRefSeq']
    ref=[r for r in ref if r['name2']=='PMP22' and r['name'].startswith('NM_')]
    promoter_exons=sorted({(int(r['exonStarts'].strip(',').split(',')[-1]),int(r['exonEnds'].strip(',').split(',')[-1])) for r in ref})
    cages=read(OUT/'cage_locus.grch38.tsv');sensitivity=[]
    for a,b in promoter_exons:
        if a not in {15265153,15262428}:continue
        for donor in ['CNhs12073','CNhs12345','CNhs12621']:
            sites=[(int(r['start']),int(r['end']),r['strand'],float(r['tags'])) for r in cages if r['donor']==donor]
            sensitivity.append(dict(annotation='RefSeq coding transcript first exon',promoter='PMP22_P1' if a==15265153 else 'PMP22_P2',
                exon_start=a,exon_end=b,donor=donor,window_start=a-100,window_end=b+100,
                tags=count_ctss(sites,a-100,b+100,'-'),interpretation='annotation sensitivity; no new state contrast'))
    table(OUT/'promoter_annotation_sensitivity.tsv',sensitivity)
    # Candidate negative intervals are hypotheses, never experimentally negative labels.
    seq=read_fasta(CACHE/'hg38.chr17.fa.gz')
    exclusions=IntervalIndex([(int(r[1])-1000,int(r[2])+1000) for r in bed_rows(CACHE/'ENCFF632FMT.bed.gz') if r[0]=='chr17'])
    blacklist=IntervalIndex([(int(r[1]),int(r[2])) for r in bed_rows(CACHE/'hg38-blacklist.v2.bed.gz') if r[0]=='chr17'])
    positive=next(r for r in regions if r['region_id']=='PMP22_DISTAL_C')
    length=int(positive['end'])-int(positive['start'])
    pseq=seq[int(positive['start']):int(positive['end'])]
    gc=lambda s:(s.count('G')+s.count('C'))/len(s)
    targetgc=gc(pseq);candidates=[]
    for start in range(15300000,15500000,100):
        end=start+length
        if exclusions.overlaps(start,end) or blacklist.overlaps(start,end):continue
        if any(int(r['start'])-10000<end and int(r['end'])+10000>start
               for r in regions if r['region_id']!='PMP22_DISTAL'):continue
        if any(t['txStart']<end and t['txEnd']>start for t in genes):continue
        s=seq[start:end]
        if len(s)!=length or any(b not in 'ACGT' for b in s):continue
        candidates.append((abs(gc(s)-targetgc),abs(start-int(positive['start'])),start,end,gc(s)))
    if not candidates:raise ValueError('No eligible negative interval')
    candidates.sort();chosen=[]
    for delta,distance,start,end,gcvalue in candidates:
        if delta>.03:continue
        if any(abs(start-r['start'])<5000 for r in chosen):continue
        chosen.append(dict(control_id=f'NEG_C_{len(chosen)+1}',chrom='chr17',start=start,end=end,
            matched_region='PMP22_DISTAL_C',gc_fraction=gcvalue,positive_gc=targetgc,gc_difference=delta,
            rule='length matched; GC difference <=0.03; no broad ATAC within 1kb; at least 10kb from known reporter/promoter intervals; no transcript or blacklist',
            status='proposed control; experimentally untested; no guarantee of no regulatory function'))
        if len(chosen)==2:break
    if len(chosen)!=2:raise ValueError('Insufficient independent control windows')
    table(OUT/'candidate_negative_controls.tsv',chosen)
    with (OUT/'negative_controls.grch38.bed').open('w') as f:
        for r in chosen:f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['control_id']}\n")
    with (OUT/'negative_controls.grch38.fa').open('w') as f:
        for r in chosen:f.write(f">{r['control_id']} GRCh38 chr17:{r['start']}-{r['end']} forward_reference\n{seq[r['start']:r['end']]}\n")
    print('Bulk donors:',len({r['donor'] for r in provenance}))
    print('PMP22 bulk TPM:',[(r['experiment'],r['TPM']) for r in rna if r['primary_PMP22_gene']])
    print('Controls:',[(r['control_id'],r['start'],r['end']) for r in chosen])
    # Accessible comparison sites address recruitment/activity mismatches of closed controls.
    annotation=json.loads((CACHE/'gencode50_chr17.json').read_text())['wgEncodeGencodeBasicV50']
    transcripts=IntervalIndex([(t['txStart'],t['txEnd']) for t in annotation])
    tss=[t['txStart'] if t['strand']=='+' else t['txEnd']-1 for t in annotation]
    options=[];seen=set()
    for p in bed_rows(CACHE/'ENCFF602YVR.bed.gz'):
        if p[0]!='chr17' or int(p[9])<0:continue
        center=int(p[1])+int(p[9]);start=center-length//2;end=start+length
        if (start,end) in seen:continue
        seen.add((start,end))
        if start<16265360 and end>14229778:continue
        if blacklist.overlaps(start,end) or transcripts.overlaps(start,end):continue
        distance=min(abs(center-t) for t in tss)
        if distance<=2000:continue
        s=seq[start:end]
        if len(s)!=length or any(b not in 'ACGT' for b in s):continue
        delta=abs(gc(s)-targetgc)
        if delta>.03:continue
        options.append((delta,abs(center-int(positive['start'])),start,end,gc(s),distance,
                        f'ENCFF602YVR:chr17:{p[1]}-{p[2]}'))
    accessible=[]
    for delta,_,start,end,gcvalue,distance,source in sorted(options):
        if any(abs(start-r['start'])<5000 for r in accessible):continue
        accessible.append(dict(control_id=f'ACCESSIBLE_C_{len(accessible)+1}',chrom='chr17',start=start,end=end,
            matched_region='PMP22_DISTAL_C',gc_fraction=gcvalue,positive_gc=targetgc,gc_difference=delta,
            distance_to_annotated_TSS=distance,source_interval_id=source,
            rule='409 bp summit-centered stringent ATAC; GC difference <=0.03; outside PMP22 +/-1Mb neighborhood, GENCODE50 basic transcripts, 2kb TSS windows, blacklist',
            status='proposed accessible comparison; PMP22 non-effect unproven; confirm recruitment and local modulation in experimental state'))
        if len(accessible)==2:break
    if len(accessible)!=2:raise ValueError('Insufficient accessible comparison sites')
    table(OUT/'candidate_accessible_controls.tsv',accessible)
    selected_ids={r['source_interval_id'] for r in accessible}
    selected_peaks={}
    for p in bed_rows(CACHE/'ENCFF602YVR.bed.gz'):
        identity=f'ENCFF602YVR:{p[0]}:{p[1]}-{p[2]}'
        if identity in selected_ids:
            selected_peaks[identity]=dict(source_interval_id=identity,chrom=p[0],start=int(p[1]),end=int(p[2]),
                summit=int(p[1])+int(p[9]),source_file='ENCFF602YVR')
    table(OUT/'accessible_control_source_peaks.tsv',selected_peaks.values())
    with (OUT/'accessible_controls.grch38.bed').open('w') as bed, (OUT/'accessible_controls.grch38.fa').open('w') as fa:
        for r in accessible:
            bed.write(f"chr17\t{r['start']}\t{r['end']}\t{r['control_id']}\n")
            fa.write(f">{r['control_id']} GRCh38 chr17:{r['start']}-{r['end']} forward_reference\n{seq[r['start']:r['end']]}\n")
    print('Accessible comparison sites:',[(r['control_id'],r['start'],r['end']) for r in accessible])

if __name__=='__main__':main()
