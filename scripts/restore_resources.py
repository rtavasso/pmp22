"""Restore every locked public resource, refusing changed bytes or paths."""
import argparse
import json
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pmp22_atlas.resources import LOCK, fetch

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args()
    lock=json.loads(LOCK.read_text())
    def restore(item):
        name,entry=item
        if Path(name).name!=name:raise ValueError('Non-local cache filename')
        return fetch(name,entry['url'],source=entry['source'],md5=entry.get('publisher_md5'))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        paths=list(pool.map(restore,lock.items()))
    print(f'Verified {len(paths)} locked resources.')

if __name__=='__main__':main()
