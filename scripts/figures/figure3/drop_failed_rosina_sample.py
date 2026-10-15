#!/usr/bin/env python3
"""Exclude one failed comparative sample from a wholly pending Rosina callset.
Preserves BAMs; archives the old frozen code/config before cohort changes.
Refuses active calling, already-produced calls, or any other missing BAM.
"""
import argparse
import ast
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def command(args):
    return subprocess.check_output(args,text=True).strip()


def write_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')


class CohortCounts(ast.NodeTransformer):
    def visit_Constant(self,node):
        if type(node.value) is int and node.value in (43,47):
            key='comparative' if node.value==43 else 'samples'
            return ast.copy_location(ast.Call(func=ast.Name(id='len',ctx=ast.Load()),
                args=[ast.Subscript(value=ast.Name(id='cfg',ctx=ast.Load()),
                    slice=ast.Constant(key),ctx=ast.Load())],keywords=[]),node)
        if isinstance(node.value,str):
            return ast.copy_location(ast.Constant(node.value.replace('43comparative','42comparative')),node)
        return node


def update_engine(text):
    tree=ast.parse(text)
    wanted={'validate','call_window','combine_genome','publish'}
    found=set()
    for i,node in enumerate(tree.body):
        if isinstance(node,ast.FunctionDef) and node.name in wanted:
            tree.body[i]=CohortCounts().visit(node);found.add(node.name)
    if found!=wanted:raise RuntimeError('Unexpected frozen worker implementation')
    ast.fix_missing_locations(tree)
    result=ast.unparse(tree)+'\n'
    compile(result,'updated frozen worker','exec')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--exclude',default='ERR1307105')
    a=p.parse_args();root=a.run.resolve(strict=True)
    cfg=json.loads((root/'run.json').read_text())
    jobs=json.loads((root/'jobs.json').read_text())
    row=[r for r in cfg['samples'] if r['sample']==a.exclude]
    if len(row)!=1 or row[0]['role']!='comparative':
        raise RuntimeError('Exclusion is absent, duplicated or focal; no changes made')
    if len(cfg['comparative'])!=43 or len(cfg['samples'])!=47:
        raise RuntimeError('Expected original 43 comparative + four focal; no changes made')
    if a.exclude not in cfg['comparative']:raise RuntimeError('Cohort mismatch')
    required={'validate','call','chrom','combine','extract','focal','publish'}
    if set(jobs)!=required:raise RuntimeError('Unexpected jobs.json stages')
    ids=[]
    for value in jobs.values():ids.extend([value] if isinstance(value,str) else value)
    ids=list(dict.fromkeys(str(x) for x in ids))
    for job in ids:
        states=command(['squeue','-h','-j',job,'-o','%T']).splitlines()
        if not states or any(s.strip()!='PENDING' for s in states):
            raise RuntimeError('All downstream jobs must still be PENDING: '+job)
    if (root/'complete.json').exists() or any((root/'windows').glob('*.complete.json')):
        raise RuntimeError('Calling has already produced results; refusing mixed cohorts')
    keep=[r for r in cfg['samples'] if r['sample']!=a.exclude]
    missing=[]
    for r in keep:
        bam=Path(r['bam'])
        if not bam.is_file() or not Path(str(bam)+'.bai').is_file():missing.append(r['sample'])
    if missing:raise RuntimeError('Other BAMs/indexes missing; no changes made: '+','.join(missing))
    engine=root/'rosina_whole_genome_bcftools.py'
    original=engine.read_text();patched=update_engine(original)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    backup=root/('before_exclude_'+a.exclude+'_'+stamp)
    backup.mkdir()
    paths=[engine,root/'run.json',root/'bam_manifest.tsv',root/'comparative_bams.list']
    for path in paths:
        if not path.is_file():raise RuntimeError('Missing frozen input '+str(path))
        shutil.copy2(path,backup/path.name)
    # Hold all downstream stages while editing. On any error they stay held;
    # files in backup allow recovery. No job is cancelled or resubmitted.
    subprocess.run(['scontrol','hold',','.join(ids)],check=True)
    try:
        cfg['samples']=keep
        cfg['comparative']=[s for s in cfg['comparative'] if s!=a.exclude]
        cfg['verified_reuse']=[s for s in cfg.get('verified_reuse',[]) if s!=a.exclude]
        cfg['pending_alignment_ids']=[int(r['array_task']) for r in keep if r['origin']=='new_oct14']
        cfg.setdefault('exclusions',[]).append(dict(sample=a.exclude,
            reason='User requested continuation without failed download',timestamp=stamp))
        engine.write_text(patched)
        cfg['frozen_script_sha256']=hashlib.sha256(engine.read_bytes()).hexdigest()
        write_json(root/'run.json',cfg)
        with (root/'bam_manifest.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(keep[0]),delimiter='\t');w.writeheader();w.writerows(keep)
        by_name={r['sample']:r for r in keep}
        (root/'comparative_bams.list').write_text(''.join(by_name[s]['bam']+'\n' for s in cfg['comparative']))
        write_json(root/'cohort_exclusion.json',dict(sample=a.exclude,comparative=42,focal=4,
            backup=str(backup),worker_sha256=cfg['frozen_script_sha256'],
            note='Frozen worker revised before any calling; old code/config preserved in backup',
            original_alignment_dependency_removed=True))
        validation=jobs['validate']
        validation=[validation] if isinstance(validation,str) else validation
        for job in validation:
            subprocess.run(['scontrol','update','JobId='+str(job),'Dependency=none'],check=True)
        subprocess.run(['scontrol','release',','.join(ids)],check=True)
    except Exception:
        print('UPDATE FAILED: downstream jobs remain held where possible. Backup: '+str(backup),file=sys.stderr)
        raise
    print('EXCLUDED\t'+a.exclude)
    print('COHORT\t42 comparative; 4 focal libraries')
    print('VALIDATION\t'+','.join(map(str,validation)))
    print('PRESERVED\tAll retained BAMs and the existing calling arrays')
    print('BACKUP\t'+str(backup))

if __name__=='__main__':main()
