"""Build a standalone, offline research report from committed analysis outputs."""
import base64
import csv
import html
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

from .resources import ROOT

DATA=ROOT/'data';HUMAN=DATA/'human';BENCH=ROOT/'benchmark/results';FIG=ROOT/'figures'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans','Arial','Liberation Sans','sans-serif'],'font.size':10,'axes.spines.top':False,
    'axes.spines.right':False,'svg.fonttype':'none','svg.hashsalt':'pmp22-release-v1'})
INK='#173a45';TEAL='#147d82';GOLD='#c08b35';BLUE='#4676ad';GRAY='#75858d'


def read(path):
    return pd.read_csv(path,sep='\t',keep_default_na=False)


def save(fig,name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight',metadata={'Date':None})
    fig.savefig(FIG/(name+'.png'),dpi=160,bbox_inches='tight')
    plt.close(fig)
    svg=(FIG/(name+'.svg')).read_text(encoding='utf-8')
    svg='\n'.join(line.rstrip() for line in svg.splitlines())+'\n'
    (FIG/(name+'.svg')).write_text(svg,encoding='utf-8',newline='\n')
    return svg[svg.index('<svg '):]


def locus_plot():
    regions=read(DATA/'regions.tsv');strict=read(HUMAN/'schwann_atac_locus.tsv')
    broad=read(HUMAN/'schwann_atac_broad_locus.tsv');bulk=read(HUMAN/'bulk_locus_peaks.tsv')
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),gridspec_kw={'width_ratios':[1,1.25]})
    for ax,lo,hi,title in [(axes[0],15246000,15267000,'Promoters and intronic element'),
                           (axes[1],15368000,15421000,'Distal regulatory neighborhood')]:
        tracks=[(strict,3,'Schwann ATAC · IDR',TEAL),(broad,2,'Schwann ATAC · broad',BLUE),
                (bulk[bulk.target.eq('H3K27ac')],1,'Bulk nerve · H3K27ac',GOLD),
                (bulk[bulk.target.eq('EP300')],0,'Bulk nerve · EP300','#a55b76')]
        for frame,y,label,color in tracks:
            for row in frame.itertuples(index=False):
                if row.start<hi and row.end>lo:
                    ax.broken_barh([(max(row.start,lo),min(row.end,hi)-max(row.start,lo))],(y-.17,.34),facecolors=color)
        names={'PMP22_INTRONIC':'Intronic','PMP22_P1':'P1','PMP22_P2':'P2',
               'PMP22_DISTAL_A':'A','PMP22_DISTAL_B':'B','PMP22_DISTAL_C':'C'}
        for r in regions.itertuples(index=False):
            if r.region_id not in names or not (lo<r.start<hi):continue
            ax.broken_barh([(r.start,r.end-r.start)],(3.95,.22),facecolors=INK)
            y=4.36 if r.region_id!='PMP22_P2' else 4.12
            ax.annotate(names[r.region_id],((r.start+r.end)/2,4.15),xytext=((r.start+r.end)/2,y),
                ha='center',fontsize=10,color=INK)
        ax.set(xlim=(lo,hi),ylim=(-.6,4.8),yticks=[4,3,2,1,0],
               yticklabels=['Reference intervals','ATAC · IDR','ATAC · broad','Bulk H3K27ac','Bulk EP300'])
        ax.set_title(title,loc='left',weight='bold',pad=15,color=INK)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x/1e6:.3f}'))
        ax.set_xlabel('GRCh38 chr17 position (Mb) →')
        ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
        ax.grid(axis='x',alpha=.12)
    fig.subplots_adjust(wspace=.48,bottom=.18)
    return save(fig,'human_locus')


def cage_plot():
    cage=read(HUMAN/'promoter_cage.tsv');cage=cage[cage.upstream_flank_bp.eq(100)]
    donors=['CNhs12073','CNhs12345','CNhs12621']
    p1=np.array([float(cage[cage.donor.eq(d)&cage.promoter.eq('PMP22_P1')].minus_strand_tags.iloc[0]) for d in donors])
    p2=np.array([float(cage[cage.donor.eq(d)&cage.promoter.eq('PMP22_P2')].minus_strand_tags.iloc[0]) for d in donors])
    fig,ax=plt.subplots(figsize=(9,3.3))
    y=np.arange(3);ax.barh(y,p2,color=TEAL,label='P2 tags');ax.barh(y,p1,left=p2,color=GOLD,label='P1 tags')
    for i,(a,b) in enumerate(zip(p1,p2)):
        ax.text(a+b+5,i,f'P1 {a:.0f}  /  P2 {b:.0f}'+('   · sparse' if i==1 else ''),va='center',fontsize=10)
    ax.set(yticks=y,yticklabels=['Donor 1 · CNhs12073','Donor 2 · CNhs12345','Donor 3 · CNhs12621'],xlim=(0,430),xlabel='Observed minus-strand CAGE tags')
    ax.invert_yaxis();ax.legend(frameon=False,loc='lower right');ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    return save(fig,'promoter_cage')


def benchmark_plot():
    metric=pd.concat([read(BENCH/'metrics.tsv'),read(BENCH/'supplemental_metrics.tsv')])
    ci=read(BENCH/'bootstrap.tsv');cim=read(BENCH/'GC_matched_bootstrap.tsv')
    cic=read(BENCH/'supplemental_bootstrap.tsv')
    names=['prevalence','GC_CpG','SOX10_EGR2_TEAD1_motifs','distance_to_annotated_TSS','conservation','4mer_sequence_logistic','permuted_training_labels']
    labels=['Prevalence','GC / CpG','Motifs','TSS distance','Conservation*','Sequence composite','Permuted labels']
    fig,axes=plt.subplots(2,1,figsize=(11,7.4));x=np.arange(len(names))
    for ax,(subset,color,title,interval) in zip(axes,[('all_test',BLUE,'Original test · n=2,120 · 78.3% of positives within 2 kb of TSS',ci),
                                               ('GC_matched_test',TEAL,'GC-matched subset · n=600 · 53.7% of positives within 2 kb of TSS',cim)]):
        values=[];lower=[];upper=[]
        for name in names:
            value=float(metric[metric.model.eq(name)&metric.subset.eq(subset)].average_precision.iloc[0])
            c=interval[interval.model.eq(name)]
            if c.empty:c=cic[cic.model.eq(name)]
            values.append(value);lower.append(value-float(c.AP_low.iloc[0]));upper.append(float(c.AP_high.iloc[0])-value)
        ax.bar(x,values,.64,color=color,yerr=[lower,upper],capsize=2,error_kw={'elinewidth':.8})
        ax.axhline(.5,color=GRAY,linestyle='--',lw=.8)
        ax.set(xticks=x,xticklabels=labels,ylim=(0,1.05),ylabel='Average precision',title=title)
        ax.tick_params(axis='x',labelrotation=12);ax.grid(axis='y',alpha=.12)
    fig.suptitle('Compare models within each panel; panels use different populations',fontsize=12)
    fig.subplots_adjust(bottom=.11,hspace=.58,top=.89)
    return save(fig,'benchmark')


def htable(frame,columns=None,rename=None):
    if columns is not None:frame=frame[columns]
    if rename:frame=frame.rename(columns=rename)
    return '<div class="table-scroll">'+frame.to_html(index=False,escape=True,border=0,classes='data-table')+'</div>'


def download(path,label):
    path=Path(path)
    mime='text/tab-separated-values' if path.suffix=='.tsv' else 'application/octet-stream'
    payload=base64.b64encode(path.read_bytes()).decode()
    return f'<a class="download" download="{html.escape(path.name)}" href="data:{mime};base64,{payload}">{html.escape(label)} ↧</a>'


def simple_markdown(text):
    """Small safe renderer for the repository's dossier/notes subset."""
    out=[];paragraph=[]
    def flush():
        if paragraph:out.append('<p>'+html.escape(' '.join(paragraph))+'</p>');paragraph.clear()
    for line in text.splitlines():
        if line.startswith('#'):
            flush();level=min(len(line)-len(line.lstrip('#'))+1,4)
            out.append(f'<h{level}>'+html.escape(line.lstrip('# ').strip())+f'</h{level}>')
        elif line.startswith('- **'):
            flush();out.append('<p class="meta">'+html.escape(line.replace('**','')[2:])+'</p>')
        elif not line.strip():flush()
        else:paragraph.append(line)
    flush();return ''.join(out)


def main():
    regions=read(DATA/'regions.tsv');activity=read(HUMAN/'region_activity.tsv')
    evidence=regions.merge(activity,on='region_id',validate='one_to_one')
    evidence=evidence[evidence.region_id.isin(['PMP22_INTRONIC','PMP22_P1','PMP22_P2','PMP22_DISTAL_A','PMP22_DISTAL_B','PMP22_DISTAL_C'])]
    evidence['Interval (GRCh38; BED)']=evidence.apply(lambda r:f"chr17:{r.start:,}–{r.end:,}",axis=1)
    evidence['Element']=evidence.region_id.str.replace('PMP22_','',regex=False)
    evidence=evidence[['Element','Interval (GRCh38; BED)','schwann_atac_unique_peak_intervals','broad_atac_unique_peak_intervals','phastcons100way_mean']].rename(
        columns={'schwann_atac_unique_peak_intervals':'IDR overlaps','broad_atac_unique_peak_intervals':'Broad overlaps','phastcons100way_mean':'Mean conservation'})
    cage=read(HUMAN/'cage_qc.tsv');bulk=read(HUMAN/'bulk_locus_rna.tsv');bulk=bulk[bulk.primary_PMP22_gene.eq(True)]
    injury=read(HUMAN/'rat_injury_summary.tsv');effects=read(DATA/'assay_effects.tsv')
    candidates=read(DATA/'release/candidate_panel.tsv');gates=read(DATA/'acceptance_audit.tsv')
    metrics=pd.concat([read(BENCH/'metrics.tsv'),read(BENCH/'supplemental_metrics.tsv')])
    metric_table=metrics[metrics.subset.isin(['all_test','GC_matched_test'])][['model','subset','n','average_precision','AUROC','Brier_score']].copy()
    for c in ['average_precision','AUROC','Brier_score']:metric_table[c]=metric_table[c].astype(float).round(3)
    composition=read(BENCH/'test_composition.tsv').copy()
    composition['label']=composition.label.map({0:'background',1:'positive'})
    composition=composition[['subset','label','n','within_2kb_TSS','median_TSS_distance','median_CpG']]
    composition['median_CpG']=composition.median_CpG.round(5)
    distal=read(BENCH/'stratified_metrics.tsv')
    distal=distal[distal.subset.eq('GC_matched_beyond_2kb_TSS')][['model','n','positive','prevalence','average_precision','AUROC']].copy()
    for c in ['prevalence','average_precision','AUROC']:distal[c]=distal[c].astype(float).round(3)
    figures=[locus_plot(),cage_plot(),benchmark_plot()]
    dossiers=''.join('<details class="dossier"><summary>'+html.escape(path.read_text(encoding='utf-8').splitlines()[0].lstrip('# '))+'</summary>'+
        simple_markdown(path.read_text(encoding='utf-8'))+'</details>' for path in sorted((ROOT/'candidate_dossiers').glob('*.md')))
    downloads=''.join(download(path,label) for path,label in [
        (DATA/'regions.tsv','Regions'),(HUMAN/'mapping_audit.tsv','Mapping audit'),
        (HUMAN/'regions.grch38.bed','Human BED'),(HUMAN/'regions.grch38.fa','Reference FASTA'),
        (DATA/'release/evidence_matrix.tsv','Evidence matrix'),(DATA/'assay_effects.tsv','Published effects'),
        (HUMAN/'promoter_cage.tsv','CAGE counts'),(DATA/'samples.tsv','Samples'),
        (DATA/'source_overlap.tsv','Source relationships'),(ROOT/'benchmark/splits.tsv','Frozen split'),
        (BENCH/'predictions.tsv','Test predictions'),(BENCH/'metrics.tsv','Primary metrics'),
        (BENCH/'supplemental_metrics.tsv','Supplemental metrics'),(HUMAN/'candidate_negative_controls.tsv','Control intervals'),
        (HUMAN/'candidate_accessible_controls.tsv','Accessible comparison sites'),
        (HUMAN/'rat_mapping_overlap_audit.tsv','Rat mapping and overlap audit'),
        (BENCH/'test_composition.tsv','Test population composition'),(BENCH/'stratified_metrics.tsv','Stratified performance'),
        (BENCH/'composition_only_scores.tsv','GC and CpG ranking checks'),
        (DATA/'acceptance_audit.tsv','Acceptance audit'),(DATA/'resource_lock.json','Resource lock')])
    sources=[
        ('Jones 2011 · intronic reporter','https://pmc.ncbi.nlm.nih.gov/articles/PMC3100536/'),
        ('Jones 2012 · distal reporters A/B/C','https://pmc.ncbi.nlm.nih.gov/articles/PMC3298281/'),
        ('Pantera 2018 · rat native deletion','https://pmc.ncbi.nlm.nih.gov/articles/PMC6077802/'),
        ('Pantera 2020 · mouse deletion and phenotype','https://pmc.ncbi.nlm.nih.gov/articles/PMC7322568/'),
        ('Hung 2015 · nerve injury','https://pmc.ncbi.nlm.nih.gov/articles/PMC4358118/'),
        ('GSE63103 · injury deposit','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63103'),
        ('Lopez-Anido 2015 · differential SOX10 occupancy','https://pmc.ncbi.nlm.nih.gov/articles/PMC4644515/'),
        ('ENCODE · pooled Schwann annotation','https://www.encodeproject.org/annotations/ENCSR301DUJ/'),
        ('FANTOM · primary-cell source files','https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.primary_cell.hCAGE/'),
        ('Borzoi · official source and split documentation','https://github.com/calico/borzoi'),
        ('UCSC · hg38 reference','https://hgdownload.soe.ucsc.edu/goldenPath/hg38/'),
        ('JASPAR · motif source','https://jaspar.elixir.no/'),
    ]
    source_html=''.join(f'<li><a href="{url}">{html.escape(label)}</a></li>' for label,url in sources)
    document=f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Project A — PMP22 atlas results</title><style>
:root{{--ink:{INK};--teal:{TEAL};--muted:#5b6d73;--line:#dce5e7;--paper:#f4f7f7;--gold:{GOLD}}}
*{{box-sizing:border-box}}html{{scroll-behavior:smooth;scroll-padding-top:76px}}
body{{margin:0;background:var(--paper);color:var(--ink);font:16px/1.62 system-ui,-apple-system,"Segoe UI",sans-serif}}
a{{color:#106e79}}a:hover{{text-decoration-thickness:2px}}main{{max-width:1220px;margin:auto;padding:30px 30px 80px}}
header{{padding:48px 0 28px;max-width:980px}}.eyebrow{{font-size:12px;letter-spacing:.13em;text-transform:uppercase;font-weight:750;color:var(--teal)}}
h1{{font-size:clamp(32px,5vw,54px);line-height:1.1;letter-spacing:-.035em;margin:14px 0 20px}}
.lead{{font-size:21px;line-height:1.5;max-width:900px}}.subtle{{color:var(--muted);font-size:14px}}
nav{{position:sticky;top:0;z-index:3;background:#f4f7f7f5;border-bottom:1px solid var(--line);padding:14px 0;display:flex;gap:24px;flex-wrap:wrap;font-size:14px}}
nav a{{text-decoration:none;font-weight:650}}section{{padding:35px 0;border-bottom:1px solid var(--line)}}h2{{font-size:29px;letter-spacing:-.02em;line-height:1.2}}h3{{font-size:20px}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:24px 0}}.card{{padding:21px;background:white;border:1px solid var(--line);border-radius:12px}}
.card strong{{display:block;font-size:29px;font-weight:720}}.card span{{font-size:14px;color:var(--muted)}}
.note{{padding:18px 22px;border-left:4px solid var(--gold);background:#faf3e7;border-radius:4px;margin:22px 0}}
.finding{{padding:20px 24px;background:#e4f1ef;border-radius:10px;margin:22px 0}}.two{{display:grid;grid-template-columns:1fr 1fr;gap:26px}}
figure{{margin:24px 0;background:white;padding:20px 16px;border-radius:12px;border:1px solid var(--line)}}figure svg{{display:block;width:100%;height:auto}}
figcaption{{font-size:14px;color:var(--muted);padding:8px 8px 0;line-height:1.5}}
.table-scroll{{overflow:auto;background:white;border:1px solid var(--line);border-radius:8px;margin:18px 0}}table{{border-collapse:collapse;width:100%;font-size:14px;line-height:1.5}}th,td{{padding:11px 14px;text-align:left!important;vertical-align:top;border-bottom:1px solid var(--line)}}th{{background:#eaf0f1;font-size:12px;letter-spacing:.025em}}tr:last-child td{{border-bottom:0}}
details{{background:white;border:1px solid var(--line);border-radius:8px;margin:12px 0;padding:15px 20px}}summary{{cursor:pointer;font-weight:650}}details p{{max-width:970px}}details h2{{font-size:22px}}
.downloads{{display:flex;flex-wrap:wrap;gap:9px}}.download{{padding:8px 12px;background:white;border:1px solid var(--line);border-radius:6px;text-decoration:none;font-size:13px}}
.meta{{font-size:14px;color:var(--muted)}}code{{font-size:.9em}}.source-list{{columns:2;padding-left:22px}}.source-list li{{break-inside:avoid;margin:8px 0}}
footer{{padding-top:28px;color:var(--muted);font-size:13px}}button{{font:inherit;border:1px solid var(--line);border-radius:7px;background:white;padding:7px 12px;cursor:pointer}}
@media(max-width:760px){{main{{padding:18px}}header{{padding-top:25px}}.cards{{grid-template-columns:1fr 1fr}}.two{{grid-template-columns:1fr}}nav{{gap:12px;padding:10px 0}}.source-list{{columns:1}}.lead{{font-size:18px}}figure{{padding:10px 3px}}}}
@media print{{body{{background:white;font-size:11pt}}main{{max-width:none;padding:0}}nav,.downloads,button{{display:none}}section{{break-before:auto}}figure,.card{{break-inside:avoid}}h2,h3{{break-after:avoid}}details{{break-inside:auto}}.cards{{grid-template-columns:repeat(4,1fr)}}}}
</style></head><body><main>
<header><div class="eyebrow">Project A · public-data computational release · 10 September 2026</div>
<h1>A clearer map of PMP22 regulation.<br>Precise limits on what it proves.</h1>
<p class="lead">Human accessibility and promoter measurements now connect to published regulatory intervals.
The strongest causal evidence is still rodent. The benchmark tests genomic accessibility, and the next
human experiments are specified in seven dossiers.</p>
<p class="subtle">Standalone report · figures and downloadable tables work offline · source links open online</p></header>
<nav aria-label="Report sections"><a href="#findings">Findings</a><a href="#atlas">Human atlas</a><a href="#promoters">Promoters</a><a href="#effects">Prior effects</a><a href="#benchmark">Benchmark</a><a href="#dossiers">Dossiers</a><a href="#audit">Audit & files</a></nav>
<section id="findings"><h2>What completing the analysis gives you</h2>
<div class="cards"><div class="card"><strong>9</strong><span>human reference intervals</span></div><div class="card"><strong>3</strong><span>separate CAGE donor files</span></div><div class="card"><strong>22</strong><span>typed published assay records</span></div><div class="card"><strong>7</strong><span>experimental dossiers</span></div></div>
<div class="finding"><strong>Accessibility-supported human nominations:</strong> distal C and the intronic element
overlap stringent Schwann accessibility peaks. A and B remain plausible, with support at the broader
threshold. This prioritizes experiments; accessibility alone does not assign a causal PMP22 effect.</div>
<div class="two"><div><h3>Measured human evidence</h3><p>Two adult female donors contribute to the pooled
Schwann ATAC product. Three cultured donors have separate CAGE files, with two informative at PMP22.
Eleven bulk-nerve assays trace to four donors and add tissue context.</p></div>
<div><h3>The model result is conditional</h3><p>The sequence classifier looks strong against random
background. GC matching changes the test population and leaves CpG differences. Within the matched set,
its advantage over a distance-to-gene baseline is uncertain.
The result supports a modest within-profile benchmark, not a disease or perturbation prediction.</p></div></div>
<div class="note"><strong>Completion boundary.</strong> This release completes the public-data analysis and
experimental handoff. The original scientific contract still requires native human noncoding perturbations,
a replicated human-state design and stronger target-assignment evidence. Those results were not generated.</div></section>
<section id="atlas"><div class="eyebrow">Observed accessibility · annotated intervals</div><h2>The human locus, with threshold sensitivity visible</h2>
<p>Six published hg18 constructs map uniquely to GRCh38 with full coverage and exact reciprocal mapping.
P1/P2 first-exon families come from GENCODE50. Human PMP22 runs on the minus strand;
the plots below show increasing chromosome coordinates, with a different scale in each panel.</p>
<figure>{figures[0]}<figcaption>Intervals are 0-based, half-open. Track bars mark peak locations, not comparable signal
heights or replicate counts. IDR and broad tracks are two thresholds of the same pooled Schwann data.
Bulk H3K27ac/EP300 add tissue-level evidence. The distal envelope is not a tested human deletion.</figcaption></figure>
{htable(evidence)}
<p><strong>A/B are threshold-sensitive, not proven inactive.</strong> Both have bulk H3K27ac support.
C and the intronic interval also overlap EP300 in bulk nerve. Neither chromatin overlap nor conservation
demonstrates a selective native effect on PMP22.</p>
<p>One 25,903 bp H3K27ac peak overlaps five annotated regions. The evidence matrix now retains
shared source-interval IDs, evidence families and donor IDs so repeated entries cannot be mistaken
for independent support. Candidate priority is an experimental choice, not measured enhancer strength.</p>
<details><summary>Coordinate and source caveats</summary><p>Printed coordinates retain a possible 1 bp convention
ambiguity until cloned inserts are sequenced or verified against original plasmids/supplements. Exported
FASTA is reference sequence. Mouse-to-human SE mapping covers only 54.5% of the mouse interval; gaps
are retained, and rodent overlap uses aligned blocks. Human contact loops were not measured.</p>
{htable(read(HUMAN/'mapping_audit.tsv'),['region_id','published_start','published_end','target_start','target_end','mapped_fraction','reciprocal_exact'])}</details>
</section>
<section id="promoters"><div class="eyebrow">Transcription initiation · cultured primary Schwann cells</div><h2>P2 dominates in the two informative CAGE donors</h2>
<figure>{figures[1]}<figcaption>Primary promoter-family windows include 100 bp flanks. Counts are observed CAGE
initiation tags, not mature-RNA isoform ratios. Donor2 has only one promoter tag and cannot support a precise
promoter comparison. No biological error bars are inferred from sequencing tags.</figcaption></figure>
<p>Increasing the upstream flank to 250 or 500 bp does not change the counts. RefSeq first-exon boundaries
give donor1 P2=284 instead of 285; all other counts agree. Cultured-source metadata do not establish
adult myelinating identity, so these data cannot measure a mature-versus-repair state effect.</p>
<div class="two"><div><h3>Library depth matters</h3>{htable(cage,['donor','total_tag_weight','ctss_records'])}</div>
<div><h3>Bulk RNA adds context</h3>{htable(bulk,['experiment','donor','TPM'])}
<p class="subtle">PMP22 gene ENSG00000109099, GENCODE29 quantification. Bulk RNA TPM is a different
measurement from CAGE tags per million; these values do not form a cell-state comparison.</p></div></div>
<details><summary>Human donor provenance and assay dependence</summary>
<p>Schwann ATAC source donors are ENCDO793LXB (female, 53) and ENCDO271OUW (female, 51).
Several bulk assays use the same donor tissues. An ENCODE matching-MD5 warning resolves to the alternate
leg-Schwann annotation ENCSR221DFU, whose peak files are identical. Do not count it as validation.
Whole-library QC is available, but Schwann-subset FRiP, cell count and TSS enrichment remain unverified.</p>
{htable(read(HUMAN/'bulk_provenance.tsv'),['experiment','assay','target','donor','age','sex','controls','selected_file'])}</details>
</section>
<section id="effects"><div class="eyebrow">Published experiments · separated by assay and species</div><h2>What the causal prior actually says</h2>
<p>Jones 2012 reports approximate reporter-activity losses of <strong>85% for A</strong> after a SOX10-site
mutation, <strong>45% for B</strong> after EGR2_1 mutation, and <strong>65% for C</strong> after SOX10_1 mutation.
These are human DNA constructs tested in rat S16 cells. They are not predicted human RNA or protein changes.
Other tested B/C mutations with no significant reduction remain in the effect table.</p>
<div class="note"><strong>Corrected intronic endpoint:</strong> the approximately 50% site-4 effect in Jones 2011
is a reduction in EGR2 fold induction. Each construct is divided by its own no-EGR2 baseline.
It is not a directly comparable 50% loss of induced reporter output. Figure 3 reports n=6;
the independence of those measurements remains unresolved. Distinct mutations and native RNA assays
retain their own outcomes and denominators. Mouse RNA statistical provenance no longer asserts an
unsupported ANOVA-only procedure.</div>
<div class="finding"><strong>Native rat deletion:</strong> Pantera 2018 reports approximately half the
WT-allele Pmp22 RNA after deleting the SE on two alleles in a roughly three-copy reporter line.
The comparison is three deletion clones versus five controls; the intact reporter allele is an internal
control. Pantera 2020 adds native mouse developmental RNA/protein and nerve-vulnerability evidence.</div>
<p>The mouse effects are retained as directions where numerical values were unavailable. No figure was
digitized and no uncertainty was invented. Human whole-copy PMP22 deletion is excluded from the
noncoding-only regulatory labels. Deletion vulnerability does not establish activation rescue.</p>
<details><summary>All 22 extracted effects, including negative results</summary>
{htable(effects,['effect_id','region_id','source_id','comparison','outcome','effect_value','effect_unit','replicates','value_origin','locator'])}</details>
<h3>Rat injury is a separate state context</h3>
{htable(injury,['track','genome_records','locus_records','rejected_coordinates'])}
<p>The 2 Mb rn5 Pmp22 window contains 83 sham and 38 injury peak records, with 25 sham-enriched and
3 injury-enriched records from the original differential files. Counts depend on peak widths and thresholds;
this is not a new differential test. Eight imprecise coordinates were excluded outside the locus.
No exact C differential peak is present. Two partially aligned sham overlaps do not demonstrate
a C-specific injury response. These whole-nerve measurements also leave cell composition unresolved.
The new mapping audit retains excluded candidate matches; accepted counts are not complete activity coverage.
The original tissue pilot reproduces <strong>35</strong> PNS versus <strong>3</strong> CNS SOX10 peaks and
28 P15 H3K27ac peaks. GSE64703 library-depth imbalance prevents treating those counts as an expression effect.</p>
</section>
<section id="benchmark"><div class="eyebrow">Retrospective genomic classification · one pooled source</div><h2>A strong raw score, a narrower biological conclusion</h2>
<p>New logistic models were trained on chromosomes 1/2, selected on chromosome 3, and tested on
chromosomes 17/22. The test has 2,120 balanced positive/background windows in 119 genomic blocks.
It is a chromosome holdout within one source, not an independent donor or intervention test.</p>
<figure>{figures[2]}<figcaption>Bars show AP; error bars are 95% intervals from 500 genomic-block bootstrap
draws, conditional on this source and sampling. The dashed line is the constructed 50% prevalence.
The GC-matched subset retains 28.3% of the original examples. The score change across panels cannot
be attributed to GC alone. *Conservation and additional diagnostics were added
after the first test run and are explicitly supplemental.</figcaption></figure>
<div class="note"><strong>The comparison that changes the interpretation:</strong> on GC-matched windows,
sequence AP is 0.784 versus 0.771 for TSS distance. Their paired AP difference is 0.013
(95% interval −0.036 to 0.062). The benchmark does not establish superiority over this simple comparator.</div>
<h3>Which examples are being compared?</h3>
{htable(composition)}
<p>Matching GC percentage changes promoter proximity and does not match CpG frequency.
In the matched subset, ranking by GC alone gives AP 0.506; CpG alone gives AP 0.679.
The original positives are predominantly promoter-proximal; random backgrounds avoid broad ATAC peaks
plus 1 kb. This is not a benchmark of enhancer function or Schwann specificity.</p>
<details><summary>Exploratory reanalysis beyond 2 kb from annotated TSSs</summary>
{htable(distal)}
<p>The frozen sequence model gives AP 0.580 and AUROC 0.746 on the matched subset beyond 2 kb from
annotated TSSs. Here the positive prevalence, and therefore the constant-score AP reference, is 0.360.
These are post-review diagnostics. Different strata change prevalence and case mix; compare models
within a stratum. Distance from a TSS does not certify enhancer identity. No models were reselected.</p></details>
<p>The sequence composite uses GC/CpG, SOX10/EGR2/TEAD1 motif maxima and canonical 4-mers. AP is a
ranking measure, not “96.6% accuracy.” Calibration applies to this balanced sample, not arbitrary genomic
loci. A ±100 bp window shift gives AP 0.964 in both directions without refitting. The permutation control
is near the prevalence baseline, not a formal significance test.</p>
<details><summary>All primary and GC-matched metrics</summary>{htable(metric_table)}</details>
<details><summary>Information tiers, leakage audit and limits</summary><p>Standardization and parameter fitting
use training data only; regularization is selected on validation AP. Identical reverse-complement sequence
hashes and overlapping selected windows are excluded. The original peak labels were produced genome-wide
by the source pipeline: this is fixed-label retrospective classification, not a train-only raw-read
peak/bias pipeline. Conservation uses a multispecies alignment and distance uses gene annotation.</p>
<p>All nine PMP22 intervals and their centered 524,288 bp contexts overlap Borzoi fold7, training in the
standard four replicas. Test fold3/validation fold4 do not include this locus. No Borzoi inference or
foundation-model head-to-head comparison is represented as independent validation. The classifier
does not predict enhancer-to-gene links, CRISPRi/a magnitude, PMP22 protein or clinical efficacy.</p>
{htable(read(BENCH/'sampling_qc.tsv'),['chrom','split','positive','negative','source_rows','merged_windows'])}</details>
</section>
<section id="dossiers"><div class="eyebrow">Falsifiable hypotheses · no model-derived human effect sizes</div><h2>Seven experiments worth distinguishing</h2>
<p>Priority is an ordinal evidence decision. Each dossier defines the interval, direction, decisive
comparison, controls and uncertainty. Dossiers include repression, a crossed-state test, retained-normal-locus
activation for deletion-HNPP, and two length/GC-matched intervals that remain experimentally untested.</p>
{dossiers}
</section>
<section id="audit"><div class="eyebrow">Reproducibility and scientific acceptance</div><h2>What is delivered, and what is still evidence-dependent</h2>
{htable(gates)}
<details><summary>Methods, replay instructions and remaining gaps</summary>
{simple_markdown((ROOT/'METHODS.md').read_text(encoding='utf-8'))}
{simple_markdown((ROOT/'data/open_gaps.md').read_text(encoding='utf-8'))}</details>
<h3>Download the underlying records</h3><div class="downloads">{downloads}</div>
<p class="subtle">These downloads are embedded in the HTML and work offline. The repository includes
all additional tables, model parameters, source snapshots, figures, tests and reproduction commands.</p>
<h3>Original sources</h3><ul class="source-list">{source_html}</ul>
<p>Full acquisition URLs, timestamps, byte counts and checksums are in the resource lock. Public API
metadata snapshots preserve the exact responses. Large raw files and published PDFs are not redistributed
in the Git repository. This integration does not certify publication novelty or native human causality.</p>
</section><footer>Project A · source-derived analysis and proposed experimental handoff ·
<a href="https://github.com/rtavasso/pmp22/tree/project-a-complete">project-a-complete branch</a> ·
<button type="button" onclick="window.print()">Print report</button></footer>
</main></body></html>"""
    (ROOT/'report.html').write_text(document,encoding='utf-8',newline='\n')
    print('Built standalone report.html and three scientific figures (SVG + PNG).')

if __name__=='__main__':main()
