import argparse, csv, gzip, importlib.util, itertools, json, os, time
from pathlib import Path
import numpy as np

ROOT=Path('analysis/genomewide_mixing_corrected_v1')
OUT=ROOT/'expanded_calibration'
CONTROLS='''ERR3790122 ERR3790123 ERR2298257 ERR2298260 ERR3790118 ERR2298201 ERR2298265 ERR2298235 ERR2298236 ERR2298216 ERR1307113 ERR1307112 ERR1307111 ERR3790124 ERR2298208 ERR2298263 ERR2298221 SRR11243339'''.split()
GRID=np.array([0,.05,.10,.15,.20,.25,.30])
FULL_PER_CHROM=50000
CAL_PER_CHROM=1000
REPLICATES_PER_TASK=20
TOTAL_TASKS=100

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

pilot=module('scripts/185_restored_model_pilot.py','pilot185')
qc=pilot.qc
pilot.OUT=OUT

def save_npz(path,**values):
    temporary=str(path)+'.new'
    with open(temporary,'wb') as handle:
        np.savez_compressed(handle,**values)
    os.replace(temporary,path)

def prepare(chrom):
    source=ROOT/'inputs'/f'chr{chrom}_NC_{140388+chrom}.1.qc.tsv.gz'
    rng=np.random.default_rng(186000+chrom)
    reservoir=[]
    eligible=0
    reasons={}

    with gzip.open(source,'rt') as handle:
        header=next(handle).rstrip('\n').split('\t')
        ci=[header.index(s) for s in CONTROLS]
        li=header.index('TN_Left_merged')
        ri=header.index('TN_Right_merged')

        for input_rows,line in enumerate(handle,1):
            row=line.rstrip('\n').split('\t')

            try:
                la,ln,ld=pilot.parse_focal_counts(row[li])
                ra,rn,rd=pilot.parse_focal_counts(row[ri])
            except ValueError:
                reasons['invalid_focal_counts']=reasons.get(
                    'invalid_focal_counts',0
                )+1
                continue

            if min(ln,rn,ld,rd)<10:
                reasons['focal_depth_below_10']=reasons.get(
                    'focal_depth_below_10',0
                )+1
                continue

            values=[qc.classify(row[i])[1] for i in ci]
            dosage=np.array(
                [v[0] if v is not None else -1 for v in values],
                dtype=np.int8
            )
            called=dosage>=0

            if called.sum()<13:
                reasons['fewer_than_13_confident_controls']=reasons.get(
                    'fewer_than_13_confident_controls',0
                )+1
                continue

            alternate=int(dosage[called].sum())

            if alternate==0 or alternate==2*int(called.sum()):
                reasons['controls_monomorphic_after_QC']=reasons.get(
                    'controls_monomorphic_after_QC',0
                )+1
                continue

            record=(
                int(row[1]),row[2],row[3],
                la,ln,ra,rn,dosage
            )
            eligible+=1

            if len(reservoir)<FULL_PER_CHROM:
                reservoir.append(record)
            else:
                replace=int(rng.integers(0,eligible))
                if replace<FULL_PER_CHROM:
                    reservoir[replace]=record

            if input_rows%500000==0:
                print(
                    f'SCAN chr={chrom} rows={input_rows} '
                    f'eligible={eligible}',
                    flush=True
                )

    if len(reservoir)<CAL_PER_CHROM:
        raise RuntimeError(
            f'Insufficient sites chr{chrom}: {len(reservoir)}'
        )

    reservoir.sort(key=lambda r:r[0])

    data=dict(
        chrom=np.full(len(reservoir),chrom,dtype=np.int8),
        pos=np.array([r[0] for r in reservoir],dtype=np.int32),
        ref=np.array([r[1] for r in reservoir]),
        alt=np.array([r[2] for r in reservoir]),
        left_alt=np.array([r[3] for r in reservoir],dtype=np.int32),
        left_total=np.array([r[4] for r in reservoir],dtype=np.int32),
        right_alt=np.array([r[5] for r in reservoir],dtype=np.int32),
        right_total=np.array([r[6] for r in reservoir],dtype=np.int32),
        controls=np.array([r[7] for r in reservoir],dtype=np.int8)
    )

    save_npz(OUT/f'chr{chrom}.npz',**data)

    cal_rng=np.random.default_rng(186500+chrom)
    chosen=np.sort(
        cal_rng.choice(
            len(reservoir),
            CAL_PER_CHROM,
            replace=False
        )
    )

    save_npz(
        OUT/f'chr{chrom}.calibration.npz',
        **{key:value[chosen] for key,value in data.items()}
    )

    audit=dict(
        chromosome=chrom,
        input_rows=input_rows,
        eligible=eligible,
        retained_full=len(reservoir),
        retained_calibration=CAL_PER_CHROM,
        rejections=reasons,
        selection='deterministic reservoir sample',
        control_rule='PL-consistent gap>=20; >=13 called; polymorphic',
        focal_depth='AD sum >=10 and FORMAT/DP >=10',
        depth_cap=None
    )

    (OUT/f'chr{chrom}.audit.json').write_text(
        json.dumps(audit,indent=2)+'\n'
    )
    print(json.dumps(audit),flush=True)

def dataset(calibration=False):
    suffix='.calibration.npz' if calibration else '.npz'
    parts=[]

    for chrom in range(1,21):
        with np.load(
            OUT/f'chr{chrom}{suffix}',
            allow_pickle=False
        ) as source:
            parts.append({
                key:source[key] for key in source.files
            })

    return {
        key:np.concatenate([part[key] for part in parts])
        for key in parts[0]
    }

def observed():
    emission=pilot.load_emission()
    results={}

    for name,is_cal in [('full',False),('calibration',True)]:
        data=dataset(is_cal)
        results[name]=pilot.fit_pair(
            data,emission,f'observed_{name}'
        )

    results['scope']='expanded PL-QC control-polymorphic autosomal sites'
    results['calibration_panel_sites']=20*CAL_PER_CHROM
    results['full_sites']=20*FULL_PER_CHROM

    (OUT/'observed_expanded.json').write_text(
        json.dumps(results,indent=2)+'\n'
    )
    print(json.dumps(results),flush=True)

def imputed(matrix,index,q,rng):
    result=matrix[:,index].copy()
    missing=result<0
    result[missing]=rng.binomial(2,q[missing])
    return result

def gamete(genotype,rng):
    return rng.binomial(1,genotype/2)

def simulations(task):
    data=dataset(True)
    matrix=data['controls']
    called=matrix>=0
    q=(
        np.where(called,matrix,0).sum(axis=1)+.5
    )/(
        2*called.sum(axis=1)+1
    )
    pairs=list(itertools.combinations(range(len(CONTROLS)),2))
    grid=list(itertools.product(GRID,GRID))
    emission=pilot.load_emission()
    rows=[]

    for local in range(REPLICATES_PER_TASK):
        replicate=task*REPLICATES_PER_TASK+local
        a,b=pairs[replicate%len(pairs)]
        cl,cr=grid[replicate%len(grid)]

        for founder_index,founder in enumerate(['three','four']):
            rng=np.random.default_rng(
                186000000+replicate*2+founder_index
            )
            ga=imputed(matrix,a,q,rng)
            gb=imputed(matrix,b,q,rng)

            if founder=='three':
                constructed=gamete(ga,rng)+gamete(gb,rng)
                if replicate%2==0:
                    gl,gr=ga,constructed
                else:
                    gl,gr=constructed,ga

                if np.any(np.abs(gl-gr)==2):
                    raise AssertionError(
                        'Three-founder construction violated'
                    )
            else:
                gl,gr=ga,gb

            synthetic=pilot.simulate_counts(
                gl,gr,float(cl),float(cr),
                data,emission,rng
            )
            fit=pilot.fit_pair(
                synthetic,
                emission,
                f'sim_t{task}_r{replicate}_{founder}'
            )

            rows.append(dict(
                task=task,
                replicate=replicate,
                generating_model=founder,
                donor_a=CONTROLS[a],
                donor_b=CONTROLS[b],
                true_c_left=cl,
                true_c_right=cr,
                genotype_contrast=float(
                    np.mean(np.abs(gl-gr)/2)
                ),
                estimated_c_left_null=fit['null']['c_left'],
                estimated_c_right_null=fit['null']['c_right'],
                estimated_c_left_alt=fit['alternative']['c_left'],
                estimated_c_right_alt=fit['alternative']['c_right'],
                epsilon=fit['alternative']['epsilon'],
                LR=fit['LR'],
                sites=len(gl)
            ))

        if (local+1)%5==0:
            print(
                f'PROGRESS task={task} replicates={local+1}',
                flush=True
            )

    path=OUT/'simulation_tasks'/f'task_{task:03d}.tsv'
    path.parent.mkdir(parents=True,exist_ok=True)

    with path.open('w') as handle:
        writer=csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter='\t'
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f'COMPLETE task={task} rows={len(rows)}',
        flush=True
    )

def collect():
    rows=[]

    for task in range(TOTAL_TASKS):
        path=OUT/'simulation_tasks'/f'task_{task:03d}.tsv'
        with path.open() as handle:
            rows.extend(
                csv.DictReader(handle,delimiter='\t')
            )

    numeric=[
        'task','replicate','true_c_left','true_c_right',
        'genotype_contrast','estimated_c_left_null',
        'estimated_c_right_null','estimated_c_left_alt',
        'estimated_c_right_alt','epsilon','LR','sites'
    ]

    for row in rows:
        for key in numeric:
            row[key]=float(row[key])

    expected=TOTAL_TASKS*REPLICATES_PER_TASK*2

    if len(rows)!=expected:
        raise RuntimeError(
            f'Unexpected simulation rows: {len(rows)}; '
            f'expected {expected}'
        )

    output=OUT/'all_simulations.tsv.gz'

    with gzip.open(str(output)+'.new','wt') as handle:
        writer=csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter='\t'
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(str(output)+'.new',output)

    observed_result=json.loads(
        (OUT/'observed_expanded.json').read_text()
    )['calibration']
    observed_lr=observed_result['LR']

    three=[
        row for row in rows
        if row['generating_model']=='three'
    ]
    four=[
        row for row in rows
        if row['generating_model']=='four'
    ]

    def near_observed(row):
        return (
            abs(
                row['true_c_left']
                - observed_result['alternative']['c_left']
            )<=.05
            and abs(
                row['true_c_right']
                - observed_result['alternative']['c_right']
            )<=.05
        )

    three_near=[
        row for row in three
        if near_observed(row)
    ]
    four_near=[
        row for row in four
        if near_observed(row)
    ]

    report=dict(
        observed_calibration_LR=observed_lr,
        simulation_sites=int(three[0]['sites']),
        three_replicates=len(three),
        four_replicates=len(four),
        empirical_three_null_p_all=(
            1+sum(row['LR']>=observed_lr for row in three)
        )/(1+len(three)),
        empirical_three_null_p_near_observed_mixing=(
            1+sum(
                row['LR']>=observed_lr
                for row in three_near
            )
        )/(1+len(three_near)),
        nearby_three_replicates=len(three_near),
        nearby_four_replicates=len(four_near),
        four_fraction_LR_above_observed_all=(
            sum(row['LR']>=observed_lr for row in four)
            /len(four)
        ),
        note=(
            'Empirical tail and power summaries only; '
            'no posterior model probability assigned'
        )
    )

    (
        OUT/'simulation_calibration_summary.json'
    ).write_text(json.dumps(report,indent=2)+'\n')

    print(json.dumps(report,indent=2),flush=True)

def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    observed_result=json.loads(
        (OUT/'observed_expanded.json').read_text()
    )
    full=dataset(False)
    cal=dataset(True)

    with gzip.open(
        OUT/'all_simulations.tsv.gz','rt'
    ) as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))

    candidates={}

    for founder in ['three','four']:
        subset=[
            row for row in rows
            if row['generating_model']==founder
            and float(row['true_c_left'])==.2
            and float(row['true_c_right'])==.2
        ]
        candidates[founder]=min(
            subset,
            key=lambda row:
                abs(float(row['estimated_c_left_alt'])-.2)
                +abs(float(row['estimated_c_right_alt'])-.2)
        )

    matrix=cal['controls']
    called=matrix>=0
    q=(
        np.where(called,matrix,0).sum(axis=1)+.5
    )/(
        2*called.sum(axis=1)+1
    )

    panels=[
        (
            'Observed: expanded full set',
            full['left_alt']/full['left_total'],
            full['right_alt']/full['right_total'],
            'observed'
        )
    ]

    for founder_index,founder in enumerate(['three','four']):
        row=candidates[founder]
        replicate=int(float(row['replicate']))
        a=CONTROLS.index(row['donor_a'])
        b=CONTROLS.index(row['donor_b'])
        rng=np.random.default_rng(
            186000000+replicate*2+founder_index
        )
        ga=imputed(matrix,a,q,rng)
        gb=imputed(matrix,b,q,rng)

        if founder=='three':
            constructed=gamete(ga,rng)+gamete(gb,rng)
            if replicate%2==0:
                gl,gr=ga,constructed
            else:
                gl,gr=constructed,ga
        else:
            gl,gr=ga,gb

        synthetic=pilot.simulate_counts(
            gl,gr,.2,.2,cal,pilot.load_emission(),rng
        )

        panels.append(
            (
                f'Simulated {founder}-founder: 20%',
                synthetic['left_alt']/synthetic['left_total'],
                synthetic['right_alt']/synthetic['right_total'],
                founder
            )
        )

    pdf=OUT/'expanded_AF_centroid_comparison.pdf'

    with PdfPages(pdf) as pages:
        figure,axes=plt.subplots(
            1,3,
            figsize=(11.7,4.9),
            sharex=True,
            sharey=True
        )
        figure.subplots_adjust(
            left=.07,right=.90,
            bottom=.18,top=.78,wspace=.15
        )

        for axis,(title,x,y,kind) in zip(axes,panels):
            axis.hexbin(
                x,y,
                gridsize=55,
                extent=(0,1,0,1),
                mincnt=1,
                cmap='Greys',
                bins='log',
                linewidths=0
            )
            axis.set_xlim(-.02,1.02)
            axis.set_ylim(-.02,1.02)
            axis.set_aspect('equal')
            axis.set_xlabel('Left ALT-read fraction')
            axis.set_title(
                title,loc='left',
                weight='bold',fontsize=10
            )
            axis.spines['top'].set_visible(False)
            axis.spines['right'].set_visible(False)

            if kind=='observed':
                mixing=observed_result['full']['alternative']
            else:
                mixing={'c_left':.2,'c_right':.2}

            emission=pilot.load_emission()

            expected_left=(
                (1-mixing['c_left'])*pilot.model.GL
                +mixing['c_left']*pilot.model.GR
            )
            expected_right=(
                mixing['c_right']*pilot.model.GL
                +(1-mixing['c_right'])*pilot.model.GR
            )

            x0=pilot.model.map_fraction(
                expected_left,
                emission['error_rate'],
                emission['bias_log_odds']
            )
            y0=pilot.model.map_fraction(
                expected_right,
                emission['error_rate'],
                emission['bias_log_odds']
            )

            axis.scatter(
                x0[pilot.COMPATIBLE],
                y0[pilot.COMPATIBLE],
                marker='s',
                s=34,
                facecolor='white',
                edgecolor='#0072B2',
                linewidth=1.1,
                zorder=4
            )

            if kind!='three':
                axis.scatter(
                    x0[pilot.OPPOSING],
                    y0[pilot.OPPOSING],
                    marker='D',
                    s=36,
                    facecolor='#D55E00',
                    edgecolor='white',
                    linewidth=.7,
                    zorder=5
                )

        axes[0].set_ylabel('Right ALT-read fraction')
        figure.suptitle(
            'Genome-wide observed and matched-depth founder benchmarks',
            fontsize=14,
            weight='bold',
            y=.96
        )
        figure.text(
            .07,.85,
            'Observed: up to 1,000,000 autosomal sites. '
            'Simulations: 20,000 chromosome-balanced sites. '
            'Log density shown independently per panel.',
            fontsize=9,
            color='#4A4A4A'
        )
        pages.savefig(figure,bbox_inches='tight')
        plt.close(figure)

    print(f'PDF={pdf}',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument(
        'mode',
        choices=[
            'prepare','observed',
            'simulate','collect','plot'
        ]
    )
    parser.add_argument('--task',type=int,default=0)
    args=parser.parse_args()

    OUT.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()

    print(
        f'START mode={args.mode} task={args.task}',
        flush=True
    )

    if args.mode=='prepare':
        prepare(args.task+1)
    elif args.mode=='observed':
        observed()
    elif args.mode=='simulate':
        simulations(args.task)
    elif args.mode=='collect':
        collect()
    else:
        plot()

    print(
        f'END elapsed_seconds={time.monotonic()-start:.1f}',
        flush=True
    )
