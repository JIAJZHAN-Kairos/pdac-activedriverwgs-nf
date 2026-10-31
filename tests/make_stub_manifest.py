import hashlib
import json
from pathlib import Path

p=Path('stub_inputs');p.mkdir(exist_ok=True)
for name in ['mutations.rds','elements.txt','ids.txt','universe.txt','membership.txt','seed.txt','expected.txt','counts.txt']:
    (p/name).write_text('synthetic fixture only\n')
def sha(name):return hashlib.sha256((p/name).read_bytes()).hexdigest()
manifest=dict(schema_version=1,max_vcpus=64,
    universe=str((p/'universe.txt').resolve()),membership=str((p/'membership.txt').resolve()),
    seed_results=str((p/'seed.txt').resolve()),expected_ids=str((p/'expected.txt').resolve()),independent_counts=str((p/'counts.txt').resolve()),
    blocks=[dict(block_id='synthetic_block',mutations=str((p/'mutations.rds').resolve()),
        elements=str((p/'elements.txt').resolve()),ids=str((p/'ids.txt').resolve()),
        mutations_sha256=sha('mutations.rds'),elements_sha256=sha('elements.txt'),ids_sha256=sha('ids.txt'))])
(p/'manifest.json').write_text(json.dumps(manifest))
