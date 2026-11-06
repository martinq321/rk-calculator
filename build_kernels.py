"""Rebuild the browser finite-cut quadrature for k=1..30; not a certificate."""
from pathlib import Path
import time, json, hashlib
import numpy as np
from rk_calculator import RealCutRk
out=Path(__file__).resolve().parent
meta=json.loads((out/'kernels.json').read_text())
for name,h in [('fine',0.05),('mesh',0.1)]:
    start=time.monotonic()
    print(f'Building {name} weights for ranks 1..30', flush=True)
    engine=RealCutRk(kmax=30,Q=128,h=h)
    a=np.vstack([engine.sigma[None,:],engine.weights]).astype('<f8')
    path=out/f'kernels_{name}.bin'
    a.tofile(path)
    info=engine.metadata()
    info['nodes_per_cell']=info['nodes']//info['Q']
    info['cutoff_start_index']=int(np.flatnonzero(engine.cells>64)[0])
    info.update(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    meta[name]=info
    (out/'kernels.json').write_text(json.dumps(meta,indent=2)+'\n')
    xs=[2,10000,1e6,1e10,1e30,1e100,1e300]
    values=engine.evaluate(xs)
    print(f'{name}: {a.shape}, {path.stat().st_size} bytes, {time.monotonic()-start:.1f}s',flush=True)
    for x,v in zip(xs,values):
        print(x, ' '.join(f'k{k}={v[k-1]:.12g}' for k in [1,2,6,10,15,20,25,30]),flush=True)
meta['bounds']['k_max']=30
meta['browser_revision']='ranks30-ungrouped-20261106-v1'
meta['large_input_note']='Scaled summation supports 1 < x <= 1e300 and ranks 1..30; not arbitrary precision or a correctly-rounded Rk certificate. High-rank cancellation and endpoint effects can be severe.'
(out/'kernels.json').write_text(json.dumps(meta,indent=2)+'\n')
header = '#pragma once\n#include <cstddef>\n'
header += 'static constexpr std::size_t KMAX=30;\n'
header += 'static constexpr std::size_t NFINE=%d;\n' % meta['fine']['nodes']
header += 'static constexpr std::size_t NMESH=%d;\n' % meta['mesh']['nodes']
header += 'static constexpr double GRAM_DEN[200]={' + ','.join(format(v,'.17g') for v in meta['gram_denominators']) + '};\n'
(out/'kernels_config.h').write_text(header)
print('Finished both kernels.',flush=True)
