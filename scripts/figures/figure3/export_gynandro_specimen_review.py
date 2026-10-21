#!/usr/bin/env python3
"""Sheet-ready specimen ledger from explicit identity and phenotype sources.
No phenotype inference from genotypes; no modifications to analysis inputs.
"""
import argparse,csv,gzip,hashlib,json,re
from collections import defaultdict
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
RUN=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
DEN=RUN/'Optix_Arnaud_review_nonPK_v2/merged'
V7=ROOT/'analysis/ybv_short_ALT_geography_v7'
META=ROOT/'data/derived/region_panel_metadata.tsv'
TRAITS={'B':'B','D':'D','R':'R','Yb_dorsal':'Yb-dorsal','Yb_ventral':'Yb-ventral','N':'N','Sb':'Sb','Br':'Br'}
def clean(v):return ' '.join(str(v if v is not None else '').split())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def binary(v):
    try:return str(int(float(v))) if float(v) in (0,1) else ''
    except (ValueError,TypeError):return ''
def read(p):
    with Path(p).open(newline='') as f:return [{clean(k):clean(v) for k,v in r.items()} for r in csv.DictReader(f,delimiter='\t')]
def write(p,rows,fields=None):
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def provided_pk_rows():
    # User-supplied specimen table in this analysis thread, retained as recorded
    # phenotype evidence rather than substituting it for effective analysis codes.
    runs=[376,375,364,353,344,343,342,341,340,339,374,373,372,371,370,369,368,367,366,365,363,362,361,360,359,358,357,356,355,354,352,351,350,349,348,347,346,345]
    names=['P1001','P1002','P1006','P1007','P1008','P1009','P1010','P1011','P1012','P1013','P1014','P1015','P2001','P2002','P2003','P2004','P2005','P2006','P2007','P2008','P2009','P2010','P2011','P3005','P3006','P3007','P3016','P3017','P3019','P4001','P4002','P4003','P4004','P4005','P4006','P4007','P4008','P4009']
    dark={'P1001','P1002','P1006','P1007','P1008','P3005','P3006','P3007','P3016','P3017','P3019','P4001'}
    pale={'P2010','P4002','P4003','P4004','P4005','P4006','P4007','P4008','P4009'}
    rows=[]
    for i,(run,name) in enumerate(zip(runs,names)):
        pk=name in dark|pale;ivory=name.startswith('P200') or name=='P2011'
        row=dict(sample=f'SRR11243{run:03d}',BioSample=f'SAMN{14276602+i}',specimen_name=name,
                 phenotype_label='PK dark' if name in dark else 'PK pale' if name in pale else 'ivory' if ivory else 'cydno / cross',
                 Ivory_state_recorded='mutant homozygote' if ivory else 'WT/-' if name in pale or name=='P4001' else 'WT/WT' if pk else '')
        if pk:row.update(B='1',D='1' if name in {'P1007','P1008','P3006','P3007','P4002','P4006'} else '0',R='1' if name in {'P3006','P3007','P4002','P4006'} else '0')
        if name=='P1008':row['recorded_genotype_note']='D/d in user specimen table'
        rows.append(row)
    return rows

class Identities:
    def __init__(self):self.parent={}
    def find(self,x):
        self.parent.setdefault(x,x)
        if self.parent[x]!=x:self.parent[x]=self.find(self.parent[x])
        return self.parent[x]
    def join(self,ids):
        ids=list(dict.fromkeys(x for x in ids if x))
        if not ids:return
        root=self.find(ids[0])
        for x in ids[1:]:self.parent[self.find(x)]=root

def export(out):
    import pandas as pd,numpy as np
    out.mkdir(parents=True,exist_ok=True);target=out/'specimens_for_review.tsv';manifest=out/'review_manifest.json'
    if manifest.exists():
        previous=json.loads(manifest.read_text());csvfile=out/'specimens_for_review.csv'
        changed=(target.exists() and sha(target)!=previous.get('review_file_sha256')) or (csvfile.exists() and previous.get('review_csv_sha256') and sha(csvfile)!=previous['review_csv_sha256'])
        if changed:print('REVIEW_PRESERVED\tExisting review table contains edits; not overwritten',flush=True);return
    identity=Identities();records=[];sources=[];presence=defaultdict(set)
    def add(kind,path,rows,id_fields):
        sources.append(dict(kind=kind,path=str(path),sha256=sha(path) if isinstance(path,Path) and path.is_file() else '',rows=len(rows)))
        for i,row in enumerate(rows):
            ids=[]
            for field in id_fields:
                for token in re.split(r'[;,|]+',row.get(field,'')):
                    token=clean(token)
                    if token and token.lower() not in {'na','nan','none','unknown','?'}:
                        accessions=re.findall(r'\b(?:ERR|SRR|DRR|SAMN|SAMEA|ERS)\d+\b',token)
                        ids.extend(accessions or [token])
            if not ids:ids=[f'{kind}_unidentified_{i}']
            identity.join(ids);records.append(dict(kind=kind,path=str(path),ids=ids,row=row))
    canonical=read(META);add('sheet',META,canonical,['Accession','SequenceID','EarthCapeID'])
    dm=read(DEN/'phenotype_manifest.tsv');add('optix',DEN/'phenotype_manifest.tsv',dm,['sample'])
    ym=read(V7/'membership_geography.tsv');add('yb',V7/'membership_geography.tsv',ym,['sample','individual_id','aliases'])
    excluded=V7/'excluded_samples.tsv'
    if excluded.is_file():add('yb_exclusion',excluded,read(excluded),['sample','individual_id','aliases'])
    manual=provided_pk_rows();add('provided_specimens','user-provided specimen table in analysis conversation',manual,['sample','BioSample','specimen_name'])
    ledger=ROOT/'analysis/yb_completed_reads_recovery_v1/resolved_phenotypes.tsv'
    if ledger.is_file():add('phenotype_ledger',ledger,read(ledger),['sample','individual_id'])
    identities=ROOT/'analysis/yb_resource_identities_v1/identity_resolution.tsv'
    if identities.is_file():add('identity_map',identities,read(identities),['bam_sample','run_accessions'])
    for path in sorted((ROOT/'analysis/yb_resource_identities_v1/ena').glob('SAMEA*.tsv')):
        add('ENA_identity',path,read(path),['run_accession','sample_accession','secondary_sample_accession'])
    inventory=ROOT/'analysis/yb_existing_resources_v1/bam_candidates.tsv'
    if inventory.is_file():
        rr=read(inventory);add('BAM_inventory',inventory,rr,['sample'])
        for r in rr:
            if Path(r.get('resolved_path',r.get('path',''))).is_file():presence[r['sample']].add('existing_BAM')
    vcf=ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf'
    if vcf.is_file():
        with vcf.open() as f:
            for line in f:
                if line.startswith('#CHROM'):
                    ids=line.rstrip().split('\t')[9:];break
            else:raise RuntimeError('Original VCF sample header missing')
        add('legacy_VCF',vcf,[dict(sample=s) for s in ids],['sample'])
        for s in ids:presence[s].add('legacy_Optix_VCF')
    bam_manifest=RUN/'bam_manifest.tsv'
    if bam_manifest.is_file():
        bm=read(bam_manifest);add('shared_bam_manifest',bam_manifest,bm,['sample'])
        for r in bm:
            if r.get('sample'):presence[r['sample']].add('shared_Rosina_BAM_manifest')
    focal=['TN2L_S1','TN2R_S2','TN3L_S3','TN3R_S4','TN_Left','TN_Right','TN_Left_T2','TN_Left_T3','TN_Right_T2','TN_Right_T3','TN_Left_merged','TN_Right_merged']
    identity.join(['GYNANDROMORPH_TN']+focal)
    records.append(dict(kind='focal_identity',path='documented focal tissue libraries',ids=['GYNANDROMORPH_TN']+focal,row={}))
    groups=defaultdict(list)
    for r in records:groups[identity.find(r['ids'][0])].append(r)
    den_groups=json.loads((DEN/'Dennis_frequency_group_membership.json').read_text())
    ba_a={s:str(int(g['phenotype'])) for g in den_groups for s in g['samples']}
    ba_b={r['sample']:binary(r.get('v')) for r in ym if binary(r.get('v'))!=''}
    cache_d=pd.read_csv(ROOT/'analysis/Dennis_figure_set_v1/all_individual_associated_dosages.tsv.gz',sep='\t',index_col=0)
    cache_y=pd.read_csv(ROOT/'analysis/Ybv_figure_set_v1/individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
    def one(values):
        values=sorted({clean(v) for v in values if clean(v)})
        return values[0] if len(values)==1 else 'CONFLICT: '+'; '.join(values) if values else ''
    rows=[];conflicts=[]
    for root,rr in groups.items():
        ids=sorted({s for r in rr for s in r['ids']});runs=[x for x in ids if re.fullmatch(r'(?:ERR|SRR|DRR)\d+',x)]
        bios=[x for x in ids if re.fullmatch(r'(?:SAMN|SAMEA)\d+',x)];focal_row='GYNANDROMORPH_TN' in ids
        stable='GYNANDROMORPH_TN' if focal_row else bios[0] if len(bios)==1 else runs[0] if runs else ids[0]
        sheets=[r['row'] for r in rr if r['kind']=='sheet'];optix=[r['row'] for r in rr if r['kind']=='optix'];yb=[r['row'] for r in rr if r['kind']=='yb'];led=[r['row'] for r in rr if r['kind']=='phenotype_ledger'];provided=[r['row'] for r in rr if r['kind']=='provided_specimens'];ena=[r['row'] for r in rr if r['kind']=='ENA_identity']
        def field(records,key):return one([r.get(key,'') for r in records])
        taxon=useful(field(optix,'taxon')) or useful(field(yb,'race')) or useful(field(sheets,'Taxon Name')) or useful(field(led,'taxon')) or useful(field(ena,'scientific_name'))
        row=dict(specimen_id=stable,specimen_name=field(provided,'specimen_name') or field(sheets,'SequenceID') or field(ena,'sample_alias') or stable,
            taxon=taxon,population=field(yb,'race') or field(optix,'group') or taxon)
        flags=[]
        for code,original in TRAITS.items():
            raw=one([binary(r.get(original)) for r in sheets]);row[code+'_sheet']=raw
            if code in ['B','D','R']:row[code+'_specimen_table']=one([binary(r.get(code)) for r in provided])
            if code in ['B','D','R']:effective=one([binary(r.get(code)) for r in optix]) or raw
            elif code=='Yb_ventral':effective=one([binary(r.get('v')) for r in yb]) or one([binary(r.get('v')) for r in led]) or raw
            elif code=='Yb_dorsal':effective=one([binary(r.get('d')) for r in led]) or raw
            else:effective=raw
            row[code]=effective
            if effective and raw and effective!=raw:flags.append(code+': sheet differs from analysis')
            if 'CONFLICT' in effective:flags.append(code+': conflicting values')
        aa={ba_a[s] for s in ids if s in ba_a};bb={ba_b[s] for s in ids if s in ba_b}
        row.update(BA_Fig3a='TRUE' if aa else 'FALSE',BA_Fig3b='TRUE' if bb else 'FALSE',
            D_used_Fig3a=one(aa),Ybv_used_Fig3b=one(bb),BA_Fig3a_accessions='; '.join(s for s in ids if s in ba_a),BA_Fig3b_accessions='; '.join(s for s in ids if s in ba_b))
        row.update(review_notes='',population_used_Fig3a=field(optix,'group'),population_used_Fig3b=field(yb,'race'),
            sex=field(sheets,'Sex') or field(sheets,'sex'),
            country=field(sheets,'Country'),latitude=field(sheets,'Latitude'),longitude=field(sheets,'Longitude'),geography_used_Fig3b=field(yb,'geography'),
            run_accessions='; '.join(runs),BioSample_accessions='; '.join(bios),SequenceID=field(sheets,'SequenceID'),EarthCapeID=field(sheets,'EarthCapeID'),
            aliases='; '.join(ids),data_type=field(sheets,'Data Type'),recorded_phenotype_label=field(provided,'phenotype_label'),
            Ivory_state_recorded=field(provided,'Ivory_state_recorded'),source_note=field(sheets,'note'),recorded_genotype_note=field(provided,'recorded_genotype_note'))
        if len(bios)>1:flags.append('Multiple BioSample accessions: verify identity grouping')
        if sum(s in ba_a for s in ids)>1:flags.append('Multiple Fig3a cohort accessions grouped as one specimen: verify independence')
        if sum(s in ba_b for s in ids)>1:flags.append('Multiple Fig3b representatives grouped as one specimen: verify independence')
        if focal_row:flags.append('One focal specimen with multiple tissues; side-specific phenotypes require review')
        if not runs and not focal_row:flags.append('No run accession resolved')
        row['BA_Fig3a_exclusion_reason']='' if aa else 'focal specimen' if focal_row else 'not in frozen Dennis BA cohort'
        excluded_rows=[r['row'] for r in rr if r['kind']=='yb_exclusion']
        exclusion=one([r.get(k,'') for r in excluded_rows for k in ['reason','exclusion','exclusion_reason']])
        row['BA_Fig3b_exclusion_reason']='' if bb else 'focal specimen' if focal_row else exclusion or 'not in retained Yb(v) BA cohort'
        row['callable_fraction_Dennis_focus_24']=max((float(np.isfinite(cache_d.loc[s].to_numpy(float)).mean()) for s in ids if s in cache_d.index),default='')
        row['callable_fraction_Ybv_focus_27']=max((float(np.isfinite(cache_y.loc[s].to_numpy(float)).mean()) for s in ids if s in cache_y.index),default='')
        row['resources_recorded']='; '.join(sorted({v for s in ids for v in presence[s]}))
        for k in dict.fromkeys(k for r in sheets for k in r):
            row['sheet_'+re.sub(r'[^A-Za-z0-9]+','_',k).strip('_')]=field(sheets,k)
        row['metadata_sources']='; '.join(sorted({r['path'] for r in rr}));row['review_flags']='; '.join(flags)
        for fieldname in ['taxon','population',*TRAITS]:
            if str(row[fieldname]).startswith('CONFLICT:'):conflicts.append(dict(specimen_id=stable,field=fieldname,values=row[fieldname]))
        rows.append(row)
    rows.sort(key=lambda r:(r['taxon'],r['population'],r['specimen_name'],r['specimen_id']))
    if len({r['specimen_id'] for r in rows})!=len(rows):raise RuntimeError('Duplicate stable review ID')
    primary=['specimen_id','specimen_name','taxon','population',*TRAITS,'sex','recorded_phenotype_label','Ivory_state_recorded','BA_Fig3a','BA_Fig3b','review_notes']
    fields=primary+[k for r in rows for k in r if k not in primary];fields=list(dict.fromkeys(fields))
    write(target,rows,fields)
    with (out/'specimens_for_review.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    meanings={'specimen_id':'Stable identifier. Please retain unchanged so returned edits can be matched.',
        'specimen_name':'Editable display name.','taxon':'Editable taxon name.','population':'Editable population/race display name.',
        'BA_Fig3a':'TRUE/FALSE: member of the current Dennis BA cohort. Editable request for future inclusion/exclusion.',
        'BA_Fig3b':'TRUE/FALSE: member of the current Yb(v) BA cohort. Editable request for future inclusion/exclusion.',
        'B':'Red forewing band code: 1 present, 0 absent.','D':'Dennis red basal patch: 1 present, 0 absent.',
        'R':'Ray pattern: 1 present, 0 absent.','Yb_dorsal':'Dorsal yellow band: 1 present, 0 absent.',
        'Yb_ventral':'Ventral yellow band: 1 present, 0 absent.',
        'N':'Original N phenotype code; authors should confirm biological definition.',
        'Sb':'Original Sb phenotype code; authors should confirm biological definition.',
        'Br':'Original Br phenotype code; authors should confirm biological definition.',
        'review_notes':'Add decisions, primary phenotype evidence, exclusions or identity corrections here.'}
    write(out/'column_dictionary.tsv',[dict(column=k,meaning=meanings.get(k,'Source/audit column; retained for checking against the editable fields.')) for k in fields])
    write(out/'metadata_conflicts.tsv',conflicts,['specimen_id','field','values'])
    (out/'README.txt').write_text('Import specimens_for_review.tsv into Google Sheets using Tab as separator. CSV is also supplied.\nEdit names, taxon, population, phenotype codes, BA flags and review_notes. Keep specimen_id unchanged. Blank phenotype = unknown, not absence. CONFLICT means records disagree and requires review.\nFig3a = Dennis BA; Fig3b = Yb(v) BA. Flags mean cohort membership, not that every SNP was callable. *_sheet fields preserve source-sheet codes; *_used_Fig3* preserve the effective analysis values. *_specimen_table fields copy codes from the supplied inferred-genotype table for author checking, without treating them as additional measured phenotypes. Source conflicts and aliases are retained.\nRun aliases are joined only using explicit source identity links. The TN tissue libraries are one documented gynandromorph specimen; side-specific phenotypes remain for author review. Multiple BioSample IDs are flagged. Metadata-only specimens are retained for review; presence in this table does not imply successful sequencing.\nReturn the edited TSV or CSV. Edits are not applied to analyses automatically. This exporter will preserve a review TSV if its content has been manually changed.\n')
    result=dict(status='COMPLETE',specimens=len(rows),Fig3a_cohort_accessions=len(ba_a),Fig3b_cohort_representatives=len(ba_b),
        Fig3a_review_rows=sum(r['BA_Fig3a']=='TRUE' for r in rows),Fig3b_review_rows=sum(r['BA_Fig3b']=='TRUE' for r in rows),
        conflicts=len(conflicts),sources=sources,review_file_sha256=sha(target),review_csv_sha256=sha(out/'specimens_for_review.csv'),phenotypes_inferred_from_genotypes=False)
    review,repair=repair_export(out,rows)
    result.update(review_schema='compact_with_evidence',repair=repair,review_file_sha256=sha(target),review_csv_sha256=sha(out/'specimens_for_review.csv'))
    manifest.write_text(json.dumps(result,indent=2)+'\n');print('SPECIMEN_REVIEW\t'+str(target),flush=True)


REVIEW_FIELDS=['specimen_id','specimen_name','taxon','population','B','D','R',
    'Yb_dorsal','Yb_ventral','N','Sb','Br','recorded_phenotype_label',
    'Ivory_state_recorded','BA_Fig3a','BA_Fig3b','D_used_Fig3a','Ybv_used_Fig3b',
    'country','geography_used_Fig3b','run_accessions','BioSample_accessions',
    'data_type','review_status','review_notes']

def useful(value):
    value=clean(value)
    return '' if value.lower() in {'unresolved taxon','unknown','unresolved','na','nan','none','?'} else value

def repair_export(out, original):
    """Presentation/evidence repair only. Frozen cohorts and genotype inputs are never edited."""
    out.mkdir(parents=True,exist_ok=True)
    manual={r['specimen_name']:r for r in provided_pk_rows()}
    labels={'P1009':'cydno galanthus','P1010':'cydno x pale PK',
        'P1011':'F1 pale PK x cydno','P1012':'F1 pale PK x cydno',
        'P1013':'BC of F1 x cydno','P1014':'BC of F1 x cydno','P1015':'BC of F1 x cydno'}
    countries={'French Guiana','Colombia','Panama','Panamá','Venezuela','Peru','Ecuador','Brazil','Costa Rica'}
    review=[];evidence=[];issues=[];changes=[];members=[]
    for source in original:
        r=dict(source);sid=r['specimen_id'];name=r['specimen_name'];m=manual.get(name,{})
        flags=[x for x in clean(r.get('review_flags')).split('; ') if x]
        def flag(message):
            if message not in flags:flags.append(message)
        def change(field,value,reason):
            old=r.get(field,'')
            if old!=value:
                changes.append(dict(specimen_id=sid,field=field,previous=old,corrected=value,reason=reason))
                r[field]=value
        if m:
            change('recorded_phenotype_label',labels.get(name,m['phenotype_label']),'User-provided specimen table')
            if not useful(r.get('taxon')):
                cohort='Butterfly world, piano keys' if m['phenotype_label'].startswith('PK') else 'Butterfly world, ivory stock' if m['phenotype_label']=='ivory' else 'Heliconius cydno ssp. galanthus' if name=='P1009' else 'Butterfly world, cydno / PK cross'
                change('taxon',cohort,'Recorded stock/phenotype identification; not independent taxonomic verification')
                if name!='P1009':flag('Taxon uses laboratory stock label; formal species unverified')
        pop=next((useful(r.get(k)) for k in ['population','population_used_Fig3b','population_used_Fig3a','taxon'] if useful(r.get(k))), '?')
        change('population',pop,'Known name takes precedence over placeholder')
        for trait in TRAITS:
            value=useful(r.get(trait));supplied=useful(r.get(trait+'_specimen_table'))
            provenance='Existing ledger / analysis annotation' if value else ''
            if not value and supplied:
                value=supplied;provenance='Author-supplied specimen table (recorded B/D/R code)'
                flag(trait+': restored author-supplied code; existing analysis unchanged')
            if value and supplied and value!=supplied:flag(trait+': analysis and author-supplied table disagree')
            if not value:value='?';provenance='Unresolved in exported evidence; requires historical source reconciliation'
            change(trait,value,'Restore explicit supplied annotation or label unknown')
            r[trait+'_review_source']=provenance
        if m and m['phenotype_label'].startswith('PK') and any(r[k]=='?' for k in ('Yb_dorsal','Yb_ventral')):
            flag('PK Yb surface annotations missing from exported sources; do not infer from optix/ivory genotype')
        for field in ['country','sheet_Country']:
            if r.get(field):change(field,r[field].replace('PanamÃ¡','Panamá'),'Correct UTF-8 mojibake')
        bad=r.get('sheet_Called_sites','')
        if bad:
            try:float(bad)
            except (ValueError,TypeError):
                flag('Malformed historical source statistics; original fields retained, no numeric values invented')
                pubs=[v for k,v in r.items() if k.startswith('sheet_Publication')]
                if not useful(r.get('country')) and len(pubs)==1 and pubs[0] in countries:
                    change('country',pubs[0],'Country explicitly recorded in malformed source row')
        if not r.get('resources_recorded'):flag('No resource recorded in ledger; data absence not established')
        for fig in ['3a','3b']:
            accessions=[s for s in str(r.get('BA_Fig'+fig+'_accessions','')).split('; ') if s]
            if len(accessions)>1:flag('Fig'+fig+': multiple cohort accessions share this specimen identity; verify independence/duplicate contribution')
            for accession in accessions:
                members.append(dict(figure='Fig'+fig,cohort_accession=accession,specimen_id=sid,
                    specimen_name=name,phenotype_used=r.get('D_used_Fig3a' if fig=='3a' else 'Ybv_used_Fig3b',''),
                    identity_review='Verify duplicate contribution' if len(accessions)>1 else ''))
        missing=[k for k in TRAITS if r[k]=='?']
        if missing:flag('Unresolved phenotype scores: '+', '.join(missing))
        if not useful(r.get('taxon')):change('taxon','?','No taxon recovered');flag('Taxon unresolved')
        r['review_status']='; '.join(flags) or 'No issue identified by ledger checks'
        issues.extend(dict(specimen_id=sid,specimen_name=name,issue=f) for f in flags)
        review.append({k:r.get(k,'') for k in REVIEW_FIELDS})
        audit=dict(source)
        for k,v in r.items():
            if k not in source or source.get(k)!=v:audit['corrected_'+k]=v
        evidence.append(audit)
    if len({r['specimen_id'] for r in review})!=len(review):raise RuntimeError('Duplicate stable specimen ID')
    write(out/'specimens_for_review.tsv',review,REVIEW_FIELDS)
    with (out/'specimens_for_review.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=REVIEW_FIELDS);w.writeheader();w.writerows(review)
    write(out/'specimen_source_evidence.tsv',evidence)
    write(out/'specimen_review_issues.tsv',issues,['specimen_id','specimen_name','issue'])
    write(out/'specimen_review_changes.tsv',changes,['specimen_id','field','previous','corrected','reason'])
    write(out/'specimen_BA_membership.tsv',members,['figure','cohort_accession','specimen_id','specimen_name','phenotype_used','identity_review'])
    meanings={'specimen_id':'Stable identifier: retain unchanged.',
        'B':'Red band: 1 present, 0 absent, ? unresolved.','D':'Dennis: 1 present, 0 absent, ? unresolved.',
        'R':'Ray: 1 present, 0 absent, ? unresolved.','Yb_dorsal':'Dorsal yellow band: 1 present, 0 absent, ? unresolved.',
        'Yb_ventral':'Ventral yellow band: 1 present, 0 absent, ? unresolved.',
        'BA_Fig3a':'Current Dennis cohort membership, not callability at every SNP.',
        'BA_Fig3b':'Current Yb(v) cohort membership, not callability at every SNP.',
        'review_notes':'Author edits/decisions.','review_status':'Specific unresolved metadata/identity issues.',
        'D_used_Fig3a':'Frozen Dennis analysis phenotype.','Ybv_used_Fig3b':'Frozen Yb(v) analysis phenotype.'}
    write(out/'column_dictionary.tsv',[dict(column=k,meaning=meanings.get(k,'Recorded review field; source evidence retained separately.')) for k in REVIEW_FIELDS])
    (out/'README.txt').write_text('1=present, 0=absent, ?=unresolved; unknown never means absence.\nFig3a=Dennis; Fig3b=Yb(v). BA flags and *_used_Fig3* record existing frozen analyses, not per-SNP callability.\nMissing B/D/R codes are restored only from explicit author-supplied codes. Missing Yb surface scores are not inferred from optix/ivory genotypes. Laboratory stock labels do not establish formal species.\nEdit review fields/notes; retain specimen_id. Nothing is applied to analyses automatically.\nSource evidence preserves original fields, including malformed historical statistics. BA membership has one row per cohort accession; shared specimen identities remain flagged rather than silently split.\n')
    summary=dict(status='REPAIRED_WITH_EXPLICIT_UNRESOLVED',specimens=len(review),columns=len(REVIEW_FIELDS),
        changes=len(changes),unknown_phenotypes={k:sum(r[k]=='?' for r in review) for k in TRAITS},
        Fig3a_cohort_accessions=sum(r['figure']=='Fig3a' for r in members),
        Fig3b_cohort_representatives=sum(r['figure']=='Fig3b' for r in members),analyses_changed=False)
    (out/'repair_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return review,summary


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();export(a.output)
if __name__=='__main__':main()
