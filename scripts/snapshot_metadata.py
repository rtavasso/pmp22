"""Preserve acquired public text/API bytes; exclude papers and large assay files."""
import gzip
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import ROOT,CACHE,LOCK,sha256

def main():
    dest=ROOT/'data/snapshots';dest.mkdir(exist_ok=True);count=0
    for name,entry in json.loads(LOCK.read_text()).items():
        if Path(name).suffix not in {'.json','.tsv','.txt','.md'}:continue
        path=CACHE/name
        if sha256(path)!=entry['sha256']:raise ValueError(name)
        (dest/(name+'.gz')).write_bytes(gzip.compress(path.read_bytes(),mtime=0))
        count+=1
    print('Preserved',count,'public metadata snapshots with deterministic gzip headers.')

if __name__=='__main__':main()
