"""Acquire frozen public inputs; data/cache is deliberately excluded from Git."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from concurrent.futures import ThreadPoolExecutor, as_completed
from pmp22_atlas.resources import fetch, encode, encode_file


def main():
    jobs = []
    with ThreadPoolExecutor(max_workers=5) as pool:
        for acc in ['ENCSR301DUJ', 'ENCSR726QTF', 'ENCSR205TUH', 'ENCFF029IGD', 'ENCFF540FVR']:
            jobs.append(pool.submit(encode, acc))
        for acc in ['ENCFF632FMT', 'ENCFF602YVR']:
            jobs.append(pool.submit(encode_file, acc))
        for sample in ['Schwann%20Cells%2c%20donor1.CNhs12073.11498-119F4',
                       'Schwann%20Cells%2c%20donor2.CNhs12345.11578-120F3',
                       'Schwann%20Cells%2c%20donor3.CNhs12621.11659-122F3']:
            name = sample.split('.')[1] + '.hg19.ctss.bed.gz'
            # The server's filenames literally contain percent-escaped text;
            # its directory links therefore escape percent signs a second time.
            url = 'https://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.primary_cell.hCAGE/' + sample.replace('%', '%25') + '.hg19.ctss.bed.gz'
            jobs.append(pool.submit(fetch, name, url, source='FANTOM5'))
        for source, target in [('hg18','Hg38'),('hg38','Hg18'),('hg19','Hg38'),('hg38','Hg19'),('rn5','Hg38')]:
            name = source + 'To' + target + '.over.chain.gz'
            jobs.append(pool.submit(fetch, name, f'https://hgdownload.soe.ucsc.edu/goldenPath/{source}/liftOver/{name}', source='UCSC liftOver'))
        for chrom in ['chr1','chr2','chr3','chr17','chr22']:
            jobs.append(pool.submit(fetch, 'hg38.'+chrom+'.fa.gz', f'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz', source='UCSC GRCh38'))
        for name,url in {
            'gencode50_locus.json': 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=wgEncodeGencodeCompV50;chrom=chr17;start=15000000;end=15600000',
            'refseq_locus.json': 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeq;chrom=chr17;start=15000000;end=15600000',
            'hg38-blacklist.v2.bed.gz': 'https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz',
            'borzoi_targets_human.txt':'https://raw.githubusercontent.com/calico/borzoi/main/examples/targets_human.txt',
            'borzoi_README.md':'https://raw.githubusercontent.com/calico/borzoi/main/README.md',
            'borzoi_commit.json':'https://api.github.com/repos/calico/borzoi/commits/main',
            'borzoi_sequences_human.bed.gz':'https://raw.githubusercontent.com/calico/borzoi/main/data/sequences_human.bed.gz',
        }.items():
            jobs.append(pool.submit(fetch, name, url, source=url.split('/')[2]))
        failures=[]
        for job in as_completed(jobs):
            try: job.result()
            except Exception as exc: failures.append(str(exc)); print('FAILED',exc,flush=True)
        if failures: raise SystemExit('\n'.join(failures))

if __name__ == '__main__': main()
