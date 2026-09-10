"""Download public resources with an immutable content lock and request ledger."""
from __future__ import annotations

import hashlib
import gzip
import json
import shutil
import threading
from datetime import datetime, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / 'data' / 'cache'
LOCK = ROOT / 'data' / 'resource_lock.json'
GUARD = threading.Lock()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fetch(name, url, *, source, md5=None):
    """First successful acquisition pins bytes. Subsequent runs verify the pin."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / name
    with GUARD:
        lock = json.loads(LOCK.read_text()) if LOCK.exists() else {}
        previous = lock.get(name)
    if path.exists() and previous:
        if sha256(path) != previous['sha256']:
            raise ValueError(f'Cached bytes changed: {name}')
        if previous['url'] != url:
            raise ValueError(f'Source URL changed: {name}')
        return path
    # API responses contain changing timestamps or annotations. Preserve the
    # acquired public metadata as compressed Git snapshots for a fresh checkout.
    snapshot=ROOT/'data'/'snapshots'/(name+'.gz')
    if previous and snapshot.exists():
        if previous['url']!=url:raise ValueError(f'Source URL changed: {name}')
        temporary=path.with_name(path.name+'.part')
        with gzip.open(snapshot,'rb') as source_handle,temporary.open('wb') as destination:
            shutil.copyfileobj(source_handle,destination)
        if sha256(temporary)!=previous['sha256']:
            raise ValueError(f'Snapshot content changed: {name}')
        temporary.replace(path)
        return path
    session = requests.Session()
    session.mount('https://', HTTPAdapter(max_retries=Retry(total=3, backoff_factor=1,
                  status_forcelist=[429, 500, 502, 503, 504])))
    response = session.get(url, timeout=(20, 90), stream=True)
    response.raise_for_status()
    temporary = path.with_name(path.name + '.part')
    with temporary.open('wb') as handle:
        for chunk in response.iter_content(1024 * 1024):
            handle.write(chunk)
    digest = sha256(temporary)
    if previous and digest != previous['sha256']:
        raise ValueError(f'Remote content changed: {name}; retained .part for audit')
    if md5:
        h = hashlib.md5()
        with temporary.open('rb') as f:
            for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
        if h.hexdigest() != md5: raise ValueError(f'Publisher MD5 mismatch: {name}')
    temporary.replace(path)
    entry = dict(url=url, resolved_url=response.url, source=source,
                 retrieved_at=datetime.now(timezone.utc).isoformat(),
                 bytes=path.stat().st_size, sha256=digest, publisher_md5=md5,
                 etag=response.headers.get('ETag'))
    with GUARD:
        lock = json.loads(LOCK.read_text()) if LOCK.exists() else {}
        lock[name] = previous or entry
        LOCK.write_text(json.dumps(lock, indent=2, sort_keys=True) + '\n')
    print(f'Acquired {name}: {path.stat().st_size:,} bytes', flush=True)
    return path


def encode(accession):
    kind = 'files' if accession.startswith('ENCFF') else 'experiments'
    path = fetch(accession + '.json', f'https://www.encodeproject.org/{kind}/{accession}/?format=json', source='ENCODE')
    return json.loads(path.read_text())


def encode_file(accession):
    metadata = encode(accession)
    filename = metadata['href'].rsplit('/', 1)[-1]
    return fetch(filename, metadata.get('cloud_metadata', {}).get('url') or
                 'https://www.encodeproject.org' + metadata['href'],
                 source=accession, md5=metadata['md5sum'])
