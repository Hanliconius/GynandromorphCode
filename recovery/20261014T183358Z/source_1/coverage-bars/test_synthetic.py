"""Run the entire pipeline on synthetic data; no Slurm required."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
import pandas as pd
import pysam

script = Path(__file__).with_name('coverage_bars.py').resolve()
with tempfile.TemporaryDirectory(prefix='coverage_bars_test_') as folder:
    root = Path(folder)
    ref = root / 'reference.fa'
    contigs = ['auto1', 'auto2', 'sexZ', 'sexW']
    ref.write_text(''.join(f'>{c}\n' + 'N'*50 + 'A'*1000 + '\n' for c in contigs))
    pysam.faidx(str(ref))
    report = root / 'assembly_report.txt'
    report.write_text(''.join(f'{c}\tassembled-molecule\t{label}\tChromosome\tna\t=\t{c}\tPrimary Assembly\t1050\n' for c,label in zip(contigs,['1','2','Z','W'])))
    for sample, depths in [('TN2L_S1',[20,20,10,10]), ('TN3L_S3',[40,40,30,10]), ('TN2R_S2',[30,30,30,0]), ('TN3R_S4',[40,40,20,10])]:
        header = {'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':c,'LN':1050} for c in contigs], 'RG':[{'ID':sample,'SM':sample}]}
        bam = root / (sample + '.bam')
        with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
            for tid, depth in enumerate(depths):
                for start in range(50,1050,50):
                    for copy in range(depth+2):
                        r = pysam.AlignedSegment()
                        r.query_name = f'r{tid}_{start}_{copy}'
                        r.query_sequence = 'A'*50
                        r.flag = 0 if copy < depth else 1024
                        r.reference_id = tid
                        r.reference_start = start
                        r.mapping_quality = 60
                        r.cigar = [(0,50)]
                        r.query_qualities = pysam.qualitystring_to_array('I'*50)
                        r.set_tag('RG',sample)
                        out.write(r)
        pysam.index(str(bam))
    wrapper = root / 'samtools'
    wrapper.write_text('#!' + sys.executable + '\nimport sys,pysam\nif sys.argv[1] == "--version":\n print("samtools " + pysam.__samtools_version__)\nelse:\n sys.stdout.write(pysam.depth(*sys.argv[2:]))\n')
    wrapper.chmod(0o755)
    run = root / 'run'
    def call(*args):
        subprocess.run([sys.executable,str(script),*map(str,args)], check=True)
    call('prepare','--project',root,'--reference',ref,'--assembly-report',report,
         '--left-t2-bam',root/'TN2L_S1.bam','--left-t3-bam',root/'TN3L_S3.bam',
         '--right-t2-bam',root/'TN2R_S2.bam','--right-t3-bam',root/'TN3R_S4.bam',
         '--samtools',wrapper,'--window','100','--run-dir',run)
    for i in range(4):
        call('scan','--run',run,'--index',i)
    call('scan','--run',run,'--index',0)
    call('plot','--run',run)
    result = pd.read_csv(run/'chromosome_coverage.tsv',sep='\t')
    assert np.allclose(result.left_t2_relative_depth,[1,1,.5,.5]), result
    assert np.allclose(result.right_t2_relative_depth,[1,1,1,0]), result
    assert np.allclose(result.left_t3_relative_depth,[1,1,.75,.25]), result
    assert np.allclose(result.right_t3_relative_depth,[1,1,.5,.25]), result
    assert (result.acgt_bases == 1000).all()
    assert (run/'chromosome_coverage.pdf').stat().st_size > 1000
    assert (run/'chromosome_coverage.png').stat().st_size > 1000
    complete = json.loads((run/'COMPLETE.json').read_text())
    assert complete['baseline'] == {'left_t2':20.,'left_t3':40.,'right_t2':30.,'right_t3':40.}
    import shutil
    shutil.copyfile(run/'chromosome_coverage.png', '/tmp/coverage_four_sample_test.png')
    print('ALL_SYNTHETIC_TESTS_PASS')
