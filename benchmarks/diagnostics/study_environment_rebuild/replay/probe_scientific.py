"""Parent assessment only: fixed CPU smoke cases, no scientific campaign/model."""
import argparse, dataclasses, hashlib, importlib, json, math, os, pathlib, sys

p=argparse.ArgumentParser();p.add_argument('--profile',choices=['research','forward-cpu','crem','mol2mol'],required=True);p.add_argument('--source',required=True);p.add_argument('--upstream');p.add_argument('--output',required=True);args=p.parse_args()
assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
assert all(os.environ.get(k)=='1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'])
root=pathlib.Path(args.source).resolve();sys.path.insert(0,str(root))
if args.upstream:sys.path.insert(0,str(pathlib.Path(args.upstream).resolve()))
import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem, rdFingerprintGenerator, rdChemReactions
from chembounce_synthesis.identity import chemical_identity
from chembounce_synthesis.common import canonical, FP
from chemistry import count_fingerprint

SMILES=['c1ccccc1','[13cH]1ccccc1','c1ccncc1','c1ccc(-c2ccccc2)cc1',
        'N[C@H](C)c1ccccc1','N[C@@H](C)c1ccccc1','Oc1ncccc1','O=c1cccc[nH]1',
        'CC(=O)O','CC(=O)[O-]','CN','C[NH3+]','C[N+](C)(C)C','CCOc1ccccc1']
identities=[dataclasses.asdict(chemical_identity(s)) for s in SMILES]
for a,b in [(0,1),(4,5),(6,7),(8,9),(10,11)]:
    assert identities[a]['family_smiles']==identities[b]['family_smiles']
    assert identities[a]['exact_smiles']!=identities[b]['exact_smiles']
assert identities[0]['scaffold_family_smiles']!=identities[2]['scaffold_family_smiles']
assert identities[0]['ring_skeleton_smiles']==identities[2]['ring_skeleton_smiles']
assert '+' in identities[12]['family_smiles']
invalid={}
for s in ['CC.O','*c1ccccc1','']:
    try:canonical(s)
    except ValueError as exc:invalid[s]=str(exc)
    else:raise AssertionError(('invalid accepted',s))
mols=[Chem.MolFromSmiles(s) for s in SMILES]
bits=[FP.GetFingerprint(m) for m in mols];counts=[count_fingerprint(m) for m in mols]
fingerprints={'bit_on_bits':[list(x.GetOnBits()) for x in bits],
              'count_sparse':[dict(x.GetNonzeroElements()) for x in counts],
              'query_bit_similarities':[DataStructs.TanimotoSimilarity(bits[0],x) for x in bits],
              'query_count_similarities':[DataStructs.TanimotoSimilarity(counts[0],x) for x in counts]}
rxn=rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])[Cl:3].[N;H2:4]>>[C:1](=[O:2])[N:4]')
products=sorted({Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(product))) for ps in rxn.RunReactants((Chem.MolFromSmiles('CC(=O)Cl'),Chem.MolFromSmiles('CN'))) for product in ps})
assert products==['CNC(C)=O']
mol=Chem.AddHs(Chem.MolFromSmiles('CCOc1ccccc1'));param=AllChem.ETKDGv3();param.randomSeed=20261016;param.numThreads=1
assert AllChem.EmbedMolecule(mol,param)==0
properties=AllChem.MMFFGetMoleculeProperties(mol);force=AllChem.MMFFGetMoleculeForceField(mol,properties);status=force.Minimize(maxIts=100)
energy=float(force.CalcEnergy());assert math.isfinite(energy)
coords=mol.GetConformer().GetPositions();assert np.isfinite(coords).all()
arithmetic=np.asarray([[1,2],[3,4]],dtype=np.float64)@np.asarray([[5,6],[7,8]],dtype=np.float64)
assert arithmetic.tolist()==[[19.,22.],[43.,50.]]
extras={};modules=['numpy','rdkit','chembounce_synthesis.identity','chembounce_synthesis.common','chemistry']
if args.profile in ['research','forward-cpu','mol2mol']:
    import torch
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    actual=(torch.tensor([[1.,2.],[3.,4.]],device='cpu')@torch.tensor([[5.,6.],[7.,8.]],device='cpu')).tolist()
    assert actual==[[19.,22.],[43.,50.]]
    extras['torch']={'version':torch.__version__,'compiled_cuda':torch.version.cuda,'cpu_matmul':actual}
    modules+=['torch']
if args.profile in ['research','forward-cpu']:
    from transformers import T5ForConditionalGeneration, AutoTokenizer
    from scipy.optimize import linear_sum_assignment
    rows,cols=linear_sum_assignment([[4,1,3],[2,0,5],[3,2,2]])
    assert cols.tolist()==[1,0,2]
    extras['linear_assignment']=[rows.tolist(),cols.tolist()]
    extras['transformers_classes']=[T5ForConditionalGeneration.__module__,AutoTokenizer.__module__]
    modules+=['transformers','scipy']
if args.profile=='crem':
    import crem.crem
    assert callable(crem.crem.mutate_mol)
    extras['crem_mutate_api']=True;modules+=['crem','crem.crem']
if args.profile=='mol2mol':
    from reinvent.models.transformer.mol2mol.mol2mol import Mol2MolModel
    from reinvent.models.transformer.core.dataset.dataset import Dataset
    assert callable(Mol2MolModel.load_from_file) and callable(Dataset.collate_fn)
    extras['official_api_classes']=[Mol2MolModel.__module__,Dataset.__module__]
    modules+=['reinvent','reinvent.models.transformer.mol2mol.mol2mol','reinvent.models.transformer.core.dataset.dataset']
origins={}
for name in modules:
    module=importlib.import_module(name);path=pathlib.Path(module.__file__).resolve()
    origins[name]={'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
for name in ['chembounce_synthesis.identity','chembounce_synthesis.common','chemistry']:
    assert pathlib.Path(origins[name]['path']).is_relative_to(root)
if args.upstream:
    for name in modules:
        if name.split('.')[0] in ['reinvent','crem']:assert pathlib.Path(origins[name]['path']).is_relative_to(pathlib.Path(args.upstream).resolve())
result={'profile':args.profile,'python':sys.version,'rdkit_version':rdBase.rdkitVersion,'identities':identities,'invalid_inputs_rejected':invalid,'fingerprints':fingerprints,'reaction_products':products,'conformer':{'seed':20261016,'minimize_status':status,'energy':energy,'positions':coords.tolist()},'numpy_arithmetic':arithmetic.tolist(),'extra_api_checks':extras,'module_origins':origins,'scope':'Fixed same-host CPU smoke assessment; no checkpoint/DB access, learned inference, GPU API call, campaign execution, or comparative endpoint.'}
destination=pathlib.Path(args.output);assert not destination.exists();destination.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
print(json.dumps({'output':str(destination),'sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),'status':'fixed_cpu_checks_passed'}))
