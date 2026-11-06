/* Rk browser evaluator: large inputs + genuine work progress.
   Uses exactly the supplied Q=128 float64 kernels, not a certified hybrid.
   No input is sent to a server; fetch only loads the fixed kernel files.
   Values are scaled by x during summation, then converted back at the end.
   Cancel is implemented by terminating this Worker from the page. */
'use strict';
const RK_VERSION = 'large-inputs-progress-20261106-v1';
const X_MAX = 1e300;
let state = null;
let busy = false;

function makeProgress(id) {
  let last = -1, lastTime = 0;
  return (percent, stage, force = false) => {
    const p = Math.min(99, Math.max(last, Math.max(0, percent)));
    const now = performance.now();
    if (force || now - lastTime > 65 || p - last >= 3) {
      self.postMessage({type:'progress', version:RK_VERSION, id, percent:p, stage});
      last = p; lastTime = now;
    }
  };
}

async function initialise(progress) {
  progress(1, 'Loading quadrature settings', true);
  const response = await fetch('kernels.json');
  if (!response.ok) throw new Error('Cannot load kernels.json. Upload all calculator data files.');
  const meta = await response.json();
  const settings = [meta.fine, meta.mesh];
  for (const m of settings) {
    if (!m || !Number.isInteger(m.nodes) || m.nodes < 1 || m.nodes > 2000000 || m.kmax !== 6) {
      throw new Error('Unsupported quadrature metadata. Use the supplied kernels.json.');
    }
  }
  const totalBytes = settings.reduce((s,m) => s + 8 * (m.kmax+1) * m.nodes, 0);
  let downloaded = 0;
  const load = async (m, label) => {
    const r = await fetch(m.file);
    if (!r.ok) throw new Error('Cannot load ' + m.file + '. Upload it beside index.html.');
    const expected = 8 * (m.kmax+1) * m.nodes;
    let buffer;
    if (r.body && typeof r.body.getReader === 'function') {
      const reader = r.body.getReader();
      const bytes = new Uint8Array(expected);
      let offset = 0;
      try {
        for (;;) {
          const {done, value} = await reader.read();
          if (done) break;
          if (offset + value.length > expected) throw new Error('Unexpected size of ' + m.file);
          bytes.set(value, offset); offset += value.length; downloaded += value.length;
          progress(4 + 31 * downloaded / totalBytes,
            'Loading ' + label + ' (' + Math.round(100*downloaded/totalBytes) + '% of data)');
        }
      } finally { reader.releaseLock(); }
      if (offset !== expected) throw new Error('Incomplete download of ' + m.file + '. Retry.');
      buffer = bytes.buffer;
    } else {
      buffer = await r.arrayBuffer(); downloaded += buffer.byteLength;
      progress(4 + 31 * downloaded / totalBytes, 'Loaded ' + label, true);
    }
    if (buffer.byteLength !== expected) throw new Error('Invalid kernel size: ' + m.file);
    // Read the little-endian file explicitly on unusual big-endian hosts.
    const littleEndian = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
    let a;
    if (littleEndian) a = new Float64Array(buffer);
    else {
      a = new Float64Array(expected/8);
      const view = new DataView(buffer);
      for (let i=0; i<a.length; ++i) a[i] = view.getFloat64(8*i, true);
    }
    for (let i=0; i<m.nodes; ++i) {
      if (!(a[i]>0 && a[i]<1)) throw new Error('Invalid quadrature node in ' + m.file);
    }
    return a;
  };
  const fine = await load(meta.fine, 'fine quadrature data');
  const mesh = await load(meta.mesh, 'comparison quadrature data');
  if (!Array.isArray(meta.gram_denominators) || meta.gram_denominators.length < 53) {
    throw new Error('The Gram-series coefficient table is incomplete.');
  }
  return {meta, fine, mesh};
}

// Kahan summation in a fixed scale: terms use exp(log(x)*(sigma-1)).
// The cutoff diagnostic sums the additional cells directly, rather than
// subtracting two nearly equal huge totals (which would round to zero).
function integrateScaled(a, n, k, logx, progress, start, end, label, withCutoff) {
  let sum=0, correction=0, low=0, lowCorrection=0, absoluteSum=0;
  const base=k*n;
  for (let i=0; i<n; ++i) {
    const term = Math.exp((a[i]-1)*logx) * a[base+i];
    if (!Number.isFinite(term)) throw new Error('Non-finite quadrature term.');
    let y=term-correction, t=sum+y;
    correction=(t-sum)-y; sum=t; absoluteSum+=Math.abs(term);
    if (withCutoff && a[i] <= 1/65) {
      y=term-lowCorrection; t=low+y;
      lowCorrection=(t-low)-y; low=t;
    }
    if ((i & 255) === 0) progress(start+(end-start)*(i+1)/n, label);
  }
  progress(end, label + ' complete', true);
  return {sum, low, absoluteSum};
}

function gramScaled(x, den, progress) {
  const logx = Math.log(x);
  let term=1/x, sum=term, correction=0, used=0;
  const cap=Math.ceil(logx+20*Math.sqrt(logx+1)+250);
  for (let n=1; n<=cap; ++n) {
    term *= logx/n;
    // Beyond the saved table, zeta(n+1)-1 is far below float64 resolution.
    const add = term/(n <= den.length ? den[n-1] : n);
    const y=add-correction, t=sum+y; correction=(t-sum)-y; sum=t;
    used=n;
    if ((n & 15) === 0) progress(36+52*n/cap, 'Summing the convergent Gram series');
    if (n > logx+10 && Math.abs(add) < 1e-16*Math.abs(sum)) {
      progress(88, 'Gram-series stopping test reached', true);
      return {sum, terms:used};
    }
  }
  throw new Error('Gram-series stopping test not reached; no result returned.');
}

function exactNk(k, n, progress) {
  if (n > 1000000) return null;
  progress(89, 'Computing the exact count for comparison', true);
  const omega=new Uint8Array(n+1);
  for (let p=2; p<=n; ++p) {
    if (omega[p] === 0) {
      for (let power=p; power<=n;) {
        for (let m=power; m<=n; m+=power) omega[m]++;
        if (power > Math.floor(n/p)) break;
        power *= p;
      }
    }
    if ((p & 4095) === 0) progress(89+7*p/n, 'Sieving the exact almost-prime count');
  }
  let count=0;
  for (let m=2; m<=n; ++m) {
    if (omega[m]===k) count++;
    if ((m & 16383) === 0) progress(96+2*m/n, 'Counting the sieve entries');
  }
  progress(98, 'Exact-count comparison complete', true);
  return count;
}

self.onmessage = async event => {
  const {id,k,x} = event.data || {};
  if (busy) {
    self.postMessage({type:'error',version:RK_VERSION,id,error:'A calculation is already running.'});
    return;
  }
  busy=true;
  const started=performance.now();
  const progress=makeProgress(id);
  try {
    if (!Number.isInteger(k) || k<1 || k>6 || !Number.isFinite(x) || x<=1 || x>X_MAX) {
      throw new Error('Use integer k = 1…6 and 1 < x ≤ 10^300.');
    }
    if (!state) state=await initialise(progress);
    else progress(35, 'Using quadrature data already loaded', true);
    const logx=Math.log(x);
    let value, cutoff=null, mesh=null, cutoffRelative=null, meshRelative=null, terms=null;
    if (k===1) {
      const g=gramScaled(x,state.meta.gram_denominators,progress);
      value=x*g.sum; terms=g.terms;
    } else {
      const fine=integrateScaled(state.fine,state.meta.fine.nodes,k,logx,progress,
        36,70,'Evaluating the fine real-cut quadrature',true);
      const coarse=integrateScaled(state.mesh,state.meta.mesh.nodes,k,logx,progress,
        70,88,'Checking the integration mesh',false);
      value=x*fine.sum;
      cutoff=x*Math.abs(fine.low);
      mesh=x*Math.abs(fine.sum-coarse.sum);
      if (fine.sum !== 0) {
        cutoffRelative=Math.abs(fine.low/fine.sum);
        meshRelative=Math.abs((fine.sum-coarse.sum)/fine.sum);
      }
    }
    if (!Number.isFinite(value)) throw new Error('The result exceeded the numerical range.');
    const count=exactNk(k,Math.floor(x),progress);
    progress(99, 'Formatting the result',true);
    self.postMessage({type:'result',version:RK_VERSION,id,k,x,value,cutoff,mesh,
      cutoffRelative,meshRelative,count,terms,certified:false,
      seconds:(performance.now()-started)/1000});
  } catch (e) {
    self.postMessage({type:'error',version:RK_VERSION,id,error:String(e.message||e)});
  } finally { busy=false; }
};
