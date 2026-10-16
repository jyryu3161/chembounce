"""Root assessment: installed REINVENT sources/data and CPU-only API imports."""
import argparse,base64,datetime,hashlib,importlib,importlib.metadata as metadata,json,os,pathlib,sys
p=argparse.ArgumentParser();p.add_argument('--upstream',type=pathlib.Path,required=True);p.add_argument('--pins',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);args=p.parse_args()
assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
assert all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'])
assert not args.output.exists()
sha=lambda q:hashlib.sha256(pathlib.Path(q).read_bytes()).hexdigest()
pins=json.loads(args.pins.read_text());source=pins['sources']['mol2mol_upstream'];upstream=args.upstream.resolve();prefix=pathlib.Path(sys.prefix).resolve()
assert upstream not in [pathlib.Path(s).resolve() for s in sys.path if s]
dist=metadata.distribution('reinvent');assert dist.version=='4.4.22';assert pathlib.Path(dist._path).resolve().is_relative_to(prefix)
expected={k:v for k,v in source['expected_files_sha256'].items() if k.startswith(('reinvent/','reinvent_plugins/')) and k.endswith('.py')}
expected.update({k:v['sha256'] for k,v in source['additional_files'].items() if k.endswith('.pkl.gz')})
files={str(x):x for x in dist.files if str(x).startswith(('reinvent/','reinvent_plugins/')) and not str(x).endswith('.pyc')}
assert set(files)==set(expected),(set(files)-set(expected),set(expected)-set(files))
actual={}
for rel,entry in sorted(files.items()):
 path=pathlib.Path(dist.locate_file(entry)).resolve();assert path.is_relative_to(prefix)
 assert sha(path)==expected[rel]==sha(upstream/rel),rel
 assert entry.hash and entry.hash.mode=='sha256' and entry.size==path.stat().st_size
 assert base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest()).rstrip(b'=').decode()==entry.hash.value
 actual[rel]={'sha256':sha(path),'bytes':path.stat().st_size}
# Class/API imports only; no checkpoint, model object, GPU API, scorer or sampling.
from reinvent.models.transformer.mol2mol.mol2mol import Mol2MolModel
from reinvent.models.transformer.core.dataset.dataset import Dataset
assert callable(Mol2MolModel.load_from_file) and callable(Dataset.collate_fn)
origins={}
for name in ['reinvent','reinvent.models.transformer.mol2mol.mol2mol','reinvent.models.transformer.core.dataset.dataset']:
 module=importlib.import_module(name);path=pathlib.Path(module.__file__).resolve();assert path.is_relative_to(prefix) and not path.is_relative_to(upstream);origins[name]={'path':str(path),'sha256':sha(path)}
# Imports must not mutate any installed/upstream scientific file.
assert all(sha(dist.locate_file(rel))==v==sha(upstream/rel) for rel,v in expected.items())
report={'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'passed':True,'profile':'mol2mol','python':sys.version,'prefix':str(prefix),'upstream':str(upstream),'upstream_commit':source['commit'],'distribution_version':dist.version,'pins_sha256':sha(args.pins),'helper_sha256':sha(__file__),'scientific_package_files':actual,'record_hashes_verified':len(actual),'all_273_python_files_and_sascore_data_match_pinned_upstream':len(actual)==274,'installed_api_origins':origins,'source_files_unchanged_after_imports':True,'scope':'Installed non-editable REINVENT package membership, bytes, RECORD and pinned upstream parity; fixed CPU-only API imports. SAScore data is hashed without unpickling/executing. No checkpoint, model object, GPU API, sampling, campaign, output ranking or general numerical equivalence.'}
assert report['all_273_python_files_and_sascore_data_match_pinned_upstream'];args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({'report':str(args.output),'sha256':sha(args.output),'passed':True,'scientific_package_files':len(actual)}))
