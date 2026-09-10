"""Inspect every bulk tibial lead named in the original contract."""
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import encode, encode_file, fetch
from pmp22_atlas.human import table, OUT

ACCESSIONS=['ENCSR771YJT','ENCSR580XUK','ENCSR401ESD','ENCSR469POZ',
    'ENCSR154GUK','ENCSR157ZVP','ENCSR479GPU','ENCSR272UNO','ENCSR648OSR','ENCSR796HLX','ENCSR858QEL']

def acquire(acc):
    data=encode(acc)
    files=[f for f in data.get('files',[]) if isinstance(f,dict)
           and f.get('status')=='released' and f.get('assembly')=='GRCh38']
    # Record all candidates; choose one frozen released peak or gene-quantification product.
    priority={'gene quantifications':0,'IDR thresholded peaks':0,'optimal IDR thresholded peaks':0,
        'replicated peaks':1,'pseudoreplicated peaks':2,'pseudoreplicated IDR thresholded peaks':2,'peaks':3}
    eligible=[f for f in files if f.get('file_format') in ['bed','tsv'] and f.get('output_type') in priority]
    # Selection uses output class, annotation release and release timestamp; never locus signal.
    eligible.sort(key=lambda f:(f.get('date_created',''),f['accession']),reverse=True)
    eligible.sort(key=lambda f:(priority[f['output_type']],-int(f.get('genome_annotation','V0').lstrip('V'))))
    if not eligible:raise ValueError('No suitable processed product: '+acc)
    chosen=eligible[0];encode_file(chosen['accession'])
    print('Selected',acc,chosen['accession'],chosen['output_type'],flush=True)
    return dict(experiment=acc,file=chosen['accession'],filename=chosen['href'].rsplit('/',1)[-1],
        output_type=chosen['output_type'],genome_annotation=chosen.get('genome_annotation',''),
        assembly=chosen['assembly'],selection='output class then newest annotation/release; no locus-based selection')

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows=list(pool.map(acquire,ACCESSIONS))
    table(OUT/'bulk_selected_files.tsv',rows)
