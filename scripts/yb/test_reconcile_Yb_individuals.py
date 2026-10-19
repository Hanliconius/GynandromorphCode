import importlib.util
from pathlib import Path
s=importlib.util.spec_from_file_location('ledger',Path(__file__).with_name('reconcile_Yb_individuals.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
assert m.resolve([dict(interval='d',value=1)],'Heliconius melpomene ssp. unknown','v')[0] is None
assert m.resolve([], 'Heliconius cydno ssp. weymeri f. gustavi','d')[0]==1
assert m.resolve([], 'Heliconius cydno ssp. weymeri f. gustavi','v')[0]==1
assert m.resolve([], 'Heliconius cydno ssp. chioneus','d')[0]==0
assert m.resolve([], 'Heliconius pachinus','v')[0]==1
assert m.resolve([], 'Heliconius melpomene ssp. vulcanus','d')[0]==0
assert m.resolve([], 'Heliconius melpomene ssp. vulcanus','v')[0]==1
r=m.resolve([dict(interval='v',value=0)],'Heliconius pachinus','v');assert r[0]==1 and r[2]
assert m.resolve([dict(interval='v',value=1),dict(interval='v',value=0)],'unresolved','v')[0] is None
print('PASS: documented cydno and surface-discordant states; conflicts retained; no dorsal-to-ventral imputation')
