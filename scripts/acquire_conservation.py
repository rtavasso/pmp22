import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import fetch

def acquire(chrom):
    name=chrom+'.phastCons100way.wigFix.gz'
    return fetch(name,'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/'+name,source='UCSC 100-way conservation')

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(acquire,['chr1','chr2','chr3','chr17','chr22']))
