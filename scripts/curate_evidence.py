"""Explicit paper transcription and metadata reconciliation, never figure digitization.

Run after human/context analysis. Each effect retains its original experimental
unit, endpoint, control and locator. Approximate author statements remain approximate.
"""
import csv
import gzip
import importlib.metadata
import json
import platform
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import ROOT,CACHE
from pmp22_atlas.human import table
from pmp22_atlas.cli import SCHEMAS

DATA=ROOT/'data'

def read(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def main():
    effects=[]
    def effect(rid,source,comparison,outcome,value,unit,n,layer,origin,species,context,interval,control,uncertainty,locator,limitations):
        values=[f'AE{len(effects)+1:03}',rid,source,comparison,outcome,value,unit,n,layer,origin,species,context,interval,control,uncertainty,locator,limitations]
        effects.append(dict(zip(SCHEMAS['assay_effects.tsv'],values)))
    for region,mutation,reduction,panel in [('A','SOX10',85,'4D'),('B','EGR2_1',45,'4E'),('C','SOX10_1',65,'4F')]:
        effect('PMP22_DISTAL_'+region,'JONES2012',mutation+'_mutant_vs_WT','luciferase_activity_reduction',
            reduction,'percent_reduction_approximate','independent luciferase n not stated in figure legend',
            'reporter','author_text_approximate','human insert in rat cells','S16 Schwann cell line',
            'historical human reporter; see mapping_audit.tsv','same-region wild-type reporter',
            'SD plotted; numerical SD and individual values unavailable','Results and Fig. '+panel,
            'Not native human RNA; C text estimate retained separately from visually approximate bar; no digitization')
    for region,mutation in [('B','SOX10'),('B','EGR2_2'),('C','SOX10_2'),('C','SOX10_3')]:
        effect('PMP22_DISTAL_'+region,'JONES2012',mutation+'_mutant_vs_WT','no significant reduction reported',
            '','not_numerically_extracted','independent luciferase n not stated in figure legend','reporter',
            'author_qualitative','human insert in rat cells','S16 Schwann cell line','same historical reporter',
            'same-region wild-type reporter','null result is not equivalence; precision unknown','Fig. 4E-F and Results',
            'Retained negative mutation; blank is not zero')
    for region,positive,total in [('A',73,78),('B',91,94),('C',79,86)]:
        effect('PMP22_DISTAL_'+region,'JONES2012','motor_nerve_positive_among_EGFP_positive_embryos',
            'embryos_with_motor_nerve_EGFP',positive,'positive_embryos',f'{total} EGFP-positive embryos',
            'in_vivo_reporter','author_exact_counts','human insert in zebrafish','transient transgenic embryos',
            'same historical reporter','conditional denominator: EGFP-positive embryos, not all injected embryos',
            'no new interval calculated; nesting by injection batch unavailable','Results preceding Fig. 5',
            'Counts measure reporter localization; not treatment-control expression or population replication')
    effect('PMP22_INTRONIC','JONES2011','EGR2_site4_mutant_vs_WT','EGR2-dependent_reporter_reduction',
        50,'percent_reduction_approximate','n not resolved for this numerical statement','reporter',
        'author_text_approximate','human insert in mouse cells','B16/F10 with EGR2 expression plasmid',
        'hg18 +11 kb reporter; see mapping_audit.tsv','wild-type enhancer with same EGR2 transfection',
        'numerical uncertainty unavailable','Results / Fig. 3B','Not native Schwann RNA; site 2 had the largest effect')
    for site in [1,3]:
        effect('PMP22_INTRONIC','JONES2011',f'EGR2_site{site}_mutant_vs_WT','activity relatively unimpaired',
            '','not_numerically_extracted','n not resolved','reporter','author_qualitative','human insert in mouse cells',
            'B16/F10 with EGR2 expression plasmid','same +11 kb reporter','wild-type reporter',
            'no equivalence margin or numerical interval','Results / Fig. 3B','Retained negative mutation')
    effect('PMP22_DISTAL','PANTERA2018','C6-C8_two_WT_allele_deletions_vs_C1-C5','WT_allele_Pmp22_RNA',
        .5,'fold_of_control_approximate','3 deletion clones vs 5 deletion-negative clones','native_RNA',
        'derived_from_author_approximate_twofold_reduction','rat','S16 F2sN reporter line; approximately three chr10 copies',
        'approximately 38 kb rat SE deletion; clone-specific junctions require supplement','deletion-negative clones; Actb normalization',
        'SD displayed per clone; raw values unavailable; no new CI','Fig. 3B-C and Results',
        'Two WT alleles deleted; reporter allele intact. Mpz normalization supports direction. Human envelope is only a region-class link')
    effect('PMP22_DISTAL','PANTERA2018','C6-C8_vs_C1-C5','reporter_allele_Pmp22_RNA_unchanged',
        '','not_numerically_extracted','3 deletion clones vs 5 negative clones','native_RNA','author_qualitative','rat',
        'same triploid reporter line','reporter allele retains SE','same clones; Actb and Mpz normalization',
        'not an equivalence test','Fig. 3B-C','Internal allelic control, not independent validation')
    for clone,endpoint,value in [('C11','WT_allele_RNA',.75),('C12','reporter_allele_RNA',.5)]:
        effect('PMP22_DISTAL','PANTERA2018',clone+'_single_allele_deletion_vs_controls',endpoint,
            value,'fold_of_control_approximate','one clone per configuration','native_RNA',
            'derived_from_author_approximate_percent_reduction','rat','S16 F2sN reporter line',
            'alternative rat SE deletion; junctions in supplement','parent line and C10 negative clone',
            'single-clone observation; no population CI','Results referring to Supplementary Fig. S3',
            'Configuration-specific support, not a calibrated gene-dosage curve')
    for age,n,outcome in [('P0','WT=3; heterozygous=6; homozygous=4','total RNA not significantly changed; homozygous P1 reduced'),
                           ('P10','WT=5; heterozygous=5; homozygous=3','total RNA reduced in both deletion genotypes'),
                           ('P56','WT=9; heterozygous=7; homozygous=4','total RNA reduced in both deletion genotypes')]:
        effect('PMP22_DISTAL','PANTERA2020',age+'_SE_deletion_vs_WT',outcome,'','direction_only',n,
            'native_RNA','author_qualitative','mouse',age+' sciatic nerve','mm10 chr11:63001867-63042405; see mouse_se_orthology.tsv',
            'wild-type littermates; Actb normalization','no raw values/CI extracted; author ANOVA only','Fig. 2A-C and Results',
            'Age-specific native mouse deletion; no human fold-change inferred; Mag/Mpz unchanged and Tvp23b increased at P0')
    effect('PMP22_DISTAL','PANTERA2020','P56_SE_deletion_vs_WT','PMP22 protein reduced','','direction_only',
        'western blot: individual nerves shown; numerical n not transcribed','protein','author_qualitative','mouse',
        'P56 sciatic nerve','same mouse SE deletion','wild-type littermates; western blot control',
        'no numerical bar digitization or uncertainty','Fig. 2D','Separate endpoint from RNA; no human protein estimate')
    effect('PMP22_DISTAL','PANTERA2020','SE_deletion_vs_WT','shorter time to 60 percent conduction block','','direction_only',
        '12 mice per genotype','function','author_qualitative','mouse','3-5 month sciatic nerve under mechanical compression',
        'same mouse SE deletion','wild-type mice under same compression','original P<0.05; no effect-size CI extracted',
        'Fig. 5A','60 percent is a challenge endpoint, not treatment efficacy; deletion creates vulnerability, not activation rescue')
    table(DATA/'assay_effects.tsv',effects,SCHEMAS['assay_effects.tsv'])

    priors=[
        ['JONES2011','Jones et al. 2011; doi:10.1523/JNEUROSCI.5893-10.2011; PMC3100536','human/rodent','reporter_and_binding','intronic regulatory mechanism','PDF inspected; intervals and selected effects extracted','Known element; no novelty claim'],
        ['JONES2012','Jones et al. 2012; doi:10.1093/hmg/ddr595; PMC3298281','human/rat/zebrafish','reporter_and_binding','A/B/C enhancers and positive/negative mutations','PDF inspected; exact printed intervals and effect text extracted','Known elements; human native function remains separate'],
        ['PANTERA2018','Pantera et al. 2018; doi:10.1093/hmg/ddy191; PMC6077802','rat','native_SE_deletion','allele-resolved endogenous RNA in triploid S16 reporter line','PDF inspected; approximate effects and controls extracted','Established causal rodent prior; clone junction supplement unresolved'],
        ['PANTERA2020','Pantera et al. 2020; doi:10.1093/hmg/ddaa082; PMC7322568','mouse','native_SE_deletion','age-specific RNA/protein and nerve vulnerability','PDF inspected; interval and genotype n extracted','Does not test activation rescue'],
        ['HUNG2015','Hung et al. 2015; doi:10.1074/jbc.M114.622878; PMC4358118','rat','injury_H3K27ac','3 day injury versus sham','PDF and GSE63103 processed files analyzed','Publication and deposit are one evidence source'],
        ['GSE63103','GEO GSE63103; Hung 2015','rat','H3K27ac_ChIP','injury and sham peaks and differential peak labels','four BED files and full SOFT downloaded','Reuse author differential labels; no new differential test'],
        ['GSE64703','GEO GSE64703 / SRP051729; Lopez-Anido et al. 2015','rat','SOX10_ChIP','PNS/CNS descriptive occupancy','12-input pilot checksums reverified','Not Ma et al.; original publication and deposit share experiments'],
        ['GSE64971','GEO GSE64971 / SRP052233; Lopez-Anido et al. 2015','rat','H3K27ac_ChIP','P15 sciatic activity','processed files and metadata analyzed','One IP and input; no biological variance estimate'],
        ['FANTOM5','FANTOM5 primary-cell hCAGE and SDRF','human','CAGE','three donor-labeled cultured Schwann samples','all three CTSS files downloaded and counted','Model donor label/path mismatch retained'],
        ['ENCSR301DUJ','ENCODE ENCSR301DUJ; ENCFF602YVR and ENCFF632FMT','human','snATAC_pseudobulk','adult tibial Schwann pooled accessibility','two original donor paths traced; both peak thresholds analyzed','Cannot provide donor-generalization or matched state effects'],
        ['GENCODE_R50','GENCODE release 50 via UCSC wgEncodeGencodeCompV50','human','annotation','minus-strand transcript and promoter families','locus and full benchmark chromosomes pinned','Annotation is not measured function'],
        ['YOSHIOKA2023','Yoshioka et al. 2023; doi:10.1038/s43856-023-00400-y; PMC10684506','human','PMP22_copy_deletion','CMT1A iPSC-derived Schwann correction','PDF inspected; excluded from enhancer effect labels','Coding-copy removal is not noncoding-only regulation'],
        ['JUN2022','PMC9410756','rodent','JUN/injury regulatory study','repair-associated regulatory context','PDF inspected; context only','Not independent human state validation'],
    ]
    table(DATA/'prior_study_inventory.tsv',[dict(zip(SCHEMAS['prior_study_inventory.tsv'],r)) for r in priors])
    samples=[r for r in read(DATA/'samples.tsv') if r['source_id']=='FANTOM5']
    for r in samples:
        if r['source_id']=='FANTOM5':
            r['download_status']='CTSS_and_SDRF_verified'
            r['permission']='public FANTOM distribution; source attribution; see data/permissions.md'
    geo=[]
    for series,study,folder in [('GSE64703','SRP051729',DATA/'raw'),('GSE64971','SRP052233',DATA/'raw'),
                                ('GSE63103','SRP049618',CACHE)]:
        runs=read(CACHE/(study+'_read_runs.tsv'))
        with gzip.open(folder/(series+'_family.soft.gz'),'rt') as f:blocks=f.read().split('^SAMPLE = ')[1:]
        for block in blocks:
            lines=block.splitlines();accession=lines[0]
            title=next(x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_title = '))
            relation=next(x.split('?term=',1)[1] for x in lines if x.startswith('!Sample_relation = SRA:'))
            linked=[r for r in runs if r['experiment_accession']==relation]
            geo.append(dict(source=series,sample_id=accession,title=title,experiment=relation,
                runs=';'.join(r['run_accession'] for r in linked),read_count=sum(int(r['read_count']) for r in linked),
                library_layout=';'.join(sorted({r['library_layout'] for r in linked})),animals_per_pool='not resolved from deposited sample metadata',
                control='input' if 'input' in title.lower() else 'IP',raw_reprocessing='not performed; archive read depth is not effective mapped depth'))
            if series=='GSE63103':continue
            tissue='spinal_cord' if 'spinal' in title.lower() else 'sciatic_nerve'
            samples.append(dict(zip(SCHEMAS['samples.tsv'],[accession,series,'rat',tissue,'P15',
                'input' if 'input' in title.lower() else 'SOX10_ChIP' if series=='GSE64703' else 'H3K27ac_ChIP',
                'control' if 'input' in title.lower() else 'IP','pooled tissue; '+title,
                'processed_SOFT_and_run_metadata_verified','public GEO'])))
    table(DATA/'human'/'geo_library_provenance.tsv',geo)
    table(DATA/'human'/'geo_run_metadata.tsv',[dict(source=study,**r) for study in ['SRP051729','SRP052233','SRP049618'] for r in read(CACHE/(study+'_read_runs.tsv'))])
    for r in read(DATA/'human'/'atac_provenance.tsv'):
        samples.append(dict(zip(SCHEMAS['samples.tsv'],[r['experiment'],'ENCSR301DUJ','human','tibial_nerve_Schwann',
            'adult; '+r['age']+' years; '+r['sex'],'snATAC','source_donor_of_pooled_product',r['donor'],
            'metadata_and_pooled_peaks_downloaded','public ENCODE'])))
    # Parse GEO sample identifiers and roles directly; do not turn IP/input into donor n.
    with gzip.open(CACHE/'GSE63103_family.soft.gz','rt') as f:
        blocks=f.read().split('^SAMPLE = ')[1:]
    for block in blocks:
        accession=block.splitlines()[0]; lines=block.splitlines()
        title=next(x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_title = '))
        samples.append(dict(zip(SCHEMAS['samples.tsv'],[accession,'GSE63103','rat','sciatic_nerve',
            'sham' if 'sham' in title.lower() else '3_days_after_transection','input' if 'input' in title.lower() else 'H3K27ac_ChIP',
            'control' if 'input' in title.lower() else 'IP','pooled_tissue; '+title,
            'processed_and_SOFT_downloaded','public GEO'])))
    if (DATA/'human'/'bulk_provenance.tsv').exists():
        for r in read(DATA/'human'/'bulk_provenance.tsv'):
            samples.append(dict(zip(SCHEMAS['samples.tsv'],[r['experiment'],r['experiment'],'human','bulk_tibial_nerve',
                r['age']+' years; '+r['sex'],r['assay'],'measurement',r['donor'],
                'metadata_and_selected_processed_file_downloaded','public ENCODE'])))
    table(DATA/'samples.tsv',samples,SCHEMAS['samples.tsv'])
    overlaps=[
        ['FANTOM5','Borzoi_targets_806_807','training_target_source','not_independent','CNhs12621 label / CNhs12073 path mismatch unresolved'],
        ['CATlas','Borzoi_kai199','training_target_source','not_independent','kai199 and ENCSR301DUJ equivalence not established'],
        ['GSE64703','Lopez_Anido_2015','deposit_for_publication','not_independent','Multiple derived tracks remain one source experiment family'],
        ['GSE64971','Lopez_Anido_2015','deposit_for_publication','not_independent','One IP plus input is not n=2'],
        ['GSE63103','HUNG2015','deposit_for_publication','not_independent','Author differential peak labels reused'],
        ['ENCSR726QTF','ENCSR301DUJ','source_alignment_ENCFF029IGD','not_independent','Female donor ENCDO793LXB; age 53'],
        ['ENCSR205TUH','ENCSR301DUJ','source_alignment_ENCFF540FVR','not_independent','Female donor ENCDO271OUW; age 51'],
        ['ENCSR301DUJ','ENCFF602YVR','IDR_peak_derivative','not_independent','One pooled Schwann product'],
        ['ENCSR301DUJ','ENCFF632FMT','broad_peak_derivative','not_independent','Threshold sensitivity, not independent replication'],
        ['ENCSR301DUJ','ENCSR221DFU','identical_peak_MD5','not_independent','Leg Schwann alternate annotation; ENCFF414RAQ/ENCFF251PBC duplicate the selected peak bytes'],
        ['ENCFF602YVR','new_sequence_benchmark','fixed_external_labels','not_independent','Chromosome holdout only; retrospective within-profile classification'],
    ]
    if (DATA/'human'/'bulk_provenance.tsv').exists():
        for r in read(DATA/'human'/'bulk_provenance.tsv'):
            overlaps.append([r['donor'],r['experiment'],'same_donor_tissue_source','not_independent',
                             'Different assays/subsections cannot be counted as independent human donors'])
    table(DATA/'source_overlap.tsv',[dict(zip(SCHEMAS['source_overlap.tsv'],r)) for r in overlaps])
    claims=[
        ['C001','Six historical human constructs map uniquely and reciprocally','reference interval','coordinate correspondence','JONES2011;JONES2012;UCSC chains','same reference resources','1 bp printed convention uncertainty; reference is not sequenced clone','Traceable human reference intervals'],
        ['C002','C and intronic element overlap stringent human Schwann ATAC peaks','pooled peak interval','observed overlap','ENCSR301DUJ','one pooled product from two donors','Threshold-sensitive A/B; no donor-specific peaks or perturbations','Human accessibility support, not enhancer-to-gene causality'],
        ['C003','P2 CAGE predominates in two informative cultured donors','donor promoter-family tag count','P1/P2 initiation tag distribution','FANTOM5','three labeled donors; donor2 sparse','Promoter window/strand rules fixed; cultured state; no mature-RNA measurement','Descriptive donor-resolved promoter initiation'],
        ['C004','Rat injury has fewer submitted locus peaks than sham','submitted peak set','window overlap count','GSE63103;HUNG2015','same experiment/deposit','Widths and thresholds differ; no new differential test','Separate rodent state context only'],
        ['C005','Reference-sequence features classify held-out ATAC windows','600 bp genomic example','average precision on fixed case-control sample','ENCSR301DUJ','one source; 119 test genomic blocks','Composition/background effects; retrospective labels; no donor holdout','Within-profile genomic classification'],
        ['C006','Standard Borzoi models saw the PMP22 locus','checkpoint genomic window','training fold overlap','Borzoi public split BED and README','four replicas share splits','Donor/path mismatch also unresolved','Ineligible for independent PMP22 genomic test'],
        ['C007','SE deletion reduces rodent native Pmp22','clone or mouse','within-assay deletion/control outcome','PANTERA2018;PANTERA2020','distinct studies; shared biological prior','Triploid line or mouse development; not human CRISPRi/a','Causal rodent regulatory prior'],
        ['C008','Selective native human PMP22 regulation is established','native human perturbation','RNA/protein change with preserved identity','NONE','missing','No qualifying perturbation result generated','NOT ALLOWED'],
    ]
    table(DATA/'claim_ledger.tsv',[dict(zip(SCHEMAS['claim_ledger.tsv'],r)) for r in claims])
    exposure=[
        ['Borzoi','standard replicas 0-3; shared split','FANTOM5;CATlas;tibial tracks','PMP22 fold7 training; centered input context audited','known or unresolved donor overlap','RNA/CAGE/ATAC','ineligible_PMP22_holdout','No checkpoint inference run; see human/borzoi_locus_exposure.tsv'],
        ['new_logistic_models','frozen coefficients and dataset lock','ENCFF602YVR','train chr1/2; val chr3; test chr17/22','same pooled product','binary accessibility','executed','Genomic holdout only; no donor/state/perturbation generalization'],
        ['AlphaGenome','not selected; no checkpoint run','unknown','not certified','not certified','multimodal','not_evaluated','No independent matched perturbation labels; no head-to-head claim'],
        ['ChromBPNet','not trained','ENCFF602YVR','not applicable','one pooled product','counts/profile','not_evaluated','Peak-only release analysis does not fit read-level bias/count/profile models'],
    ]
    table(DATA/'model_exposure.tsv',[dict(zip(SCHEMAS['model_exposure.tsv'],r)) for r in exposure])
    versions=[dict(component='Python',version=platform.python_version(),reference='CPython',purpose='executed analysis')]
    for package in ['numpy','scipy','pandas','scikit-learn','matplotlib','requests','pyliftover','pybigtools','beautifulsoup4']:
        versions.append(dict(component=package,version=importlib.metadata.version(package),reference='requirements-analysis.txt',purpose='locked analysis environment'))
    for c,v,p in [('reference_genome','UCSC hg38 primary chromosomes','FASTA source bytes SHA256 pinned; no claim that patch contigs were used'),
                  ('annotation','GENCODE 50 / bulk quantification GENCODE29','locus and promoter definitions / original RNA quantification respectively'),
                  ('Borzoi_source',json.loads((CACHE/'borzoi_commit.json').read_text())['sha'],'training exposure audit only')]:
        versions.append(dict(component=c,version=v,reference='data/resource_lock.json',purpose=p))
    table(DATA/'software_lock.tsv',versions)
    print(f'Curated {len(effects)} assay effects, {len(samples)} sample rows and {len(claims)} bounded claims.')

if __name__=='__main__':main()
