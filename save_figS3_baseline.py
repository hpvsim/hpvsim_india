"""
Freeze figS3 plot-ready values (ASR cancer incidence, HPVsim vs Globocan).

Usage: python save_figS3_baseline.py [--outdir results/v2.2.6_baseline]
"""

import argparse
import json
from datetime import date
from pathlib import Path

import hpvsim as hpv
import pandas as pd
import sciris as sc


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--resfolder', default='results')
    parser.add_argument('--outdir', default=f'results/v{hpv.__version__}_baseline')
    parser.add_argument('--globocan', default='data/globocan_asr_cancer_incidence.csv')
    parser.add_argument('--start-year', type=int, default=2005)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    res = sc.loadobj(f'{args.resfolder}/india_msim.obj')
    ind = sc.findinds(res['year'], args.start_year)[0]
    years = res['year'][ind:]

    hpvsim_df = pd.DataFrame({
        'year': years,
        'asr_cancer_incidence': res['asr_cancer_incidence'].values[ind:],
        'asr_cancer_incidence_low': res['asr_cancer_incidence'].low[ind:],
        'asr_cancer_incidence_high': res['asr_cancer_incidence'].high[ind:],
    })
    hpvsim_df.to_csv(outdir / 'figS3_hpvsim.csv', index=False)

    globocan = pd.read_csv(args.globocan)
    globocan[globocan['year'] >= args.start_year].to_csv(outdir / 'figS3_globocan.csv', index=False)

    manifest = {
        'figure': 'figS3',
        'source_objects': [f'{args.resfolder}/india_msim.obj'],
        'hpvsim_version': hpv.__version__,
        'date': date.today().isoformat(),
    }
    (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2))

    print(f'Saved figS3 baseline to {outdir}')
