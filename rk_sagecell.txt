# Self-contained Sage cell. Core executed as Python to avoid preparser changes.
exec(r'''
"""Exploratory real-cut R_k evaluator, usable in SageMath or Python.

Mathematical target: R_k(x)=integral_0^1 Im(F_k(s-i0))/pi*x**s/s ds.
This implementation returns a finite-cut numerical approximation for k>=2:
  R_k,Q(x) = integral_(1/(Q+1))^1 Delta_k(s)*x**s/s ds.
It does NOT enclose the omitted lower cut and is NOT the certified hybrid.
R_1 is instead evaluated from the convergent Gram series. Ordinary Selberg
asymptotic series are not used. All numerical integration uses float64.

Dependencies: numpy, scipy; matplotlib only for plots. No C++ compiler needed.
In Sage: load('rk_calculator.py'); print(rk_report(2, 10**6))
"""
from __future__ import annotations
import math
import csv
import json
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy.special import zetac


def _integer(value, name, lo, hi):
    n = int(value)
    if n != value or not lo <= n <= hi:
        raise ValueError('%s must be an integer in [%d,%d]' % (name, lo, hi))
    return n


def mobius_table(n):
    n = int(n)
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    prime = np.ones(n + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, n + 1):
        if prime[p]:
            mu[p::p] *= -1
            if p*p <= n:
                mu[p*p::p*p] = 0
                prime[p*p::p] = False
    return mu


def riemann_R(x):
    """Convergent Gram series, evaluated in double precision; x>=1."""
    a = np.asarray(x, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(a < 1) or np.any(a > 1e300):
        raise ValueError('This demonstration supports 1 <= x <= 1e300.')
    y = np.log(a)
    term = 1.0/a
    total = term.copy()
    correction = np.zeros_like(y)
    for n in range(1, int(np.max(y) + 20*np.sqrt(np.max(y)+1) + 251)):
        term = term*y/n
        add = term/(n*(1.0 + float(zetac(n+1.0))))
        residual = add-correction
        updated = total+residual
        correction = (updated-total)-residual
        total = updated
        if n > float(np.max(y)) + 10 and np.all(np.abs(add) < 1e-16*np.maximum(1e-308, np.abs(total))):
            break
    else:
        raise ArithmeticError('Gram series did not meet its numerical stopping rule.')
    total = total*a
    return float(total) if total.ndim == 0 else total


class RealCutRk:
    """Cache a finite-cut quadrature once, then reuse it for many x values.

    Q is the number of reciprocal cells. h is the tanh-sinh mesh spacing.
    cutoff ensures every retained Mobius sum reaches n*m*s approximately
    cutoff or larger (45 by default). These are convergence controls, NOT
    significant-digit guarantees. Results can be poor at large k or x.
    """
    def __init__(self, kmax=3, Q=64, h=0.10, cutoff=45.0, umax=3.15):
        self.kmax = _integer(kmax, 'kmax', 1, 12)
        self.Q = _integer(Q, 'Q', 8, 256)
        self.h, self.cutoff, self.umax = float(h), float(cutoff), float(umax)
        if not 0.025 <= self.h <= 0.20:
            raise ValueError('Use 0.025 <= h <= 0.20.')
        if not 30 <= self.cutoff <= 60 or not 2.8 <= self.umax <= 3.2:
            raise ValueError('Use 30 <= cutoff <= 60 and 2.8 <= umax <= 3.2.')
        M = int(math.ceil(self.cutoff*(self.Q+1)))+1
        mu = mobius_table(M)
        B = np.r_[0.0, np.cumsum(mu[1:]/np.arange(1,M+1))]
        u = np.arange(-int(self.umax/self.h), int(self.umax/self.h)+1)*self.h
        v = (math.pi/2)*np.sinh(u)
        nodes = np.tanh(v)
        jac = (math.pi/2)*np.cosh(u)/np.cosh(v)**2
        sigmas, weights, cells = [], [], []
        for q in range(1, self.Q+1):
            lo, hi = 1.0/(q+1), 1.0/q
            half = (hi-lo)/2
            sig = (hi+lo)/2 + half*nodes
            keep = (sig>lo) & (sig<hi)
            sig, ds = sig[keep], (self.h*half*jac)[keep]
            prime_zeta = []
            for m in range(1, self.kmax+1):
                nmax = int(math.ceil(self.cutoff/(m*lo)))
                ns = np.nonzero(mu[1:nmax+1])[0]+1
                arguments = (m*ns[:,None])*sig[None,:]
                # Never evaluate a rounded node exactly at a reciprocal pole.
                keep = ~np.any(arguments == 1.0, axis=0)
                if not np.all(keep):
                    sig, ds = sig[keep], ds[keep]
                    prime_zeta = [p[keep] for p in prime_zeta]
                    arguments = arguments[:,keep]
                z = zetac(arguments)
                logs = np.empty_like(z)
                above = arguments>1.0
                logs[above] = np.log1p(z[above])
                logs[~above] = np.log(np.abs(1.0+z[~above]))
                real_part = np.sum((mu[ns]/ns)[:,None]*logs, axis=0)
                bank_jump = math.pi*B[np.floor(1.0/(m*sig)).astype(int)]
                prime_zeta.append(real_part + 1j*bank_jump)
            # If exp(sum_m P(m*s)*z**m/m)=sum_n f_n*z**n,
            # n*f_n=sum_(m=1)^n P(m*s)*f_(n-m).
            f = [np.ones_like(sig,dtype=complex)]
            for n in range(1,self.kmax+1):
                f.append(sum(prime_zeta[m-1]*f[n-m] for m in range(1,n+1))/n)
            ws = np.array([z.imag*ds/(math.pi*sig) for z in f[1:]])
            if not np.all(np.isfinite(ws)):
                raise ArithmeticError('Nonfinite quadrature weights; try a coarser h.')
            sigmas.append(sig); weights.append(ws)
            cells.extend([q]*len(sig))
        self.sigma = np.concatenate(sigmas)
        self.weights = np.concatenate(weights, axis=1)
        self.cells = np.asarray(cells)

    def evaluate(self, x, full_r1=True):
        """Rows correspond to x, columns to ranks 1..kmax.

        For scalar x return a 1-D vector. Only rank 1 can use the complete
        convergent Gram series; k>=2 still omit 0<s<1/(Q+1).
        """
        original = np.asarray(x, dtype=float)
        xs = np.atleast_1d(original)
        if xs.ndim != 1 or np.any(~np.isfinite(xs)) or np.any(xs<=1) or np.any(xs>1e300):
            raise ValueError('Supply a scalar or 1-D grid with 1 < x <= 1e300.')
        values = np.empty((len(xs),self.kmax))
        for start in range(0,len(xs),64):
            stop = min(start+64,len(xs))
            kernel = np.exp(np.log(xs[start:stop,None])*(self.sigma[None,:]-1.0))
            # einsum avoids forcing a multithreaded BLAS session in SageCell.
            values[start:stop] = xs[start:stop,None]*np.einsum('ij,kj->ik',kernel,self.weights,optimize=False)
        if full_r1:
            values[:,0] = riemann_R(xs)
        return values[0] if original.ndim == 0 else values

    def metadata(self):
        return {'method':'finite positive real-cut quadrature; Gram series for rank 1',
                'Q':self.Q,'sigma_min':1.0/(self.Q+1),'h':self.h,
                'mobius_cutoff_product':self.cutoff,'umax':self.umax,
                'kmax':self.kmax,'nodes':int(len(self.sigma)),
                'arithmetic':'IEEE float64','certified_full_Rk':False,
                'warning':'For k>=2 the lower cut is omitted. Sensitivity is not an error bound.'}


@lru_cache(maxsize=3)
def _cached_engine(kmax, Q, h, cutoff):
    return RealCutRk(kmax=kmax,Q=Q,h=h,cutoff=cutoff)


def rk_report(k, x, Q=64):
    """Compare Q and 2Q, then halve the integration mesh. NOT a certificate."""
    k = _integer(k,'k',1,12)
    Q = _integer(Q,'Q',8,128)
    x = float(x)
    if not 1<x<=1e300:
        raise ValueError('This demonstration supports 1 < x <= 1e300.')
    if k == 1:
        return {'k':1,'x':x,'value':riemann_R(x),
                'method':'convergent Gram series, float64; not an interval certificate'}
    coarse = _cached_engine(k,Q,0.10,45.0)
    deep = _cached_engine(k,2*Q,0.10,45.0)
    fine = _cached_engine(k,2*Q,0.05,45.0)
    v0,v1,v2 = [float(e.evaluate(x)[k-1]) for e in [coarse,deep,fine]]
    result = {'k':k,'x':x,'value':v2,
              'cutoff_change_abs':float(abs(x*np.sum(np.exp(math.log(x)*(fine.sigma[fine.cells>Q]-1.0))*fine.weights[k-1,fine.cells>Q]))),
              'mesh_change_abs':abs(v2-v1),
              'certificate':False,**fine.metadata()}
    return result


def exact_counts(xmax, kmax=3):
    """Exact N_k(n), factors counted WITH multiplicity, for all n<=xmax."""
    xmax=_integer(xmax,'xmax',2,2000000)
    kmax=_integer(kmax,'kmax',1,20)
    prime=np.ones(xmax+1,dtype=bool); prime[:2]=False
    for p in range(2,math.isqrt(xmax)+1):
        if prime[p]: prime[p*p::p]=False
    omega=np.zeros(xmax+1,dtype=np.uint8)
    for p0 in np.flatnonzero(prime):
        p=int(p0); power=p
        while power<=xmax:
            omega[power::power] += 1
            if power > xmax//p: break
            power *= p
    return np.column_stack([np.cumsum(omega==k,dtype=np.int64)
                            for k in range(1,kmax+1)])


def first_zero_scale(k,x):
    """s_k for k=1 or 2, from the user's stated first-zero constants."""
    k = _integer(k,'k',1,2)
    t = (-1.0986500372612+0.24119591642563j)-np.log(np.log(np.asarray(x,dtype=float)))
    if k==1: return -t.real
    c2=-0.95484578810794-0.23898284814566j
    return -np.real((t*t+c2)*np.conj(t)/2)/np.abs(t)**2


def adjacent_data(xmax=1000000, step=100, Q=128, h=0.05):
    """One x-grid for all three numerical R ranks and all exact counts."""
    xmax=_integer(xmax,'xmax',10,2000000)
    step=_integer(step,'step',1,xmax-2)
    xs=np.unique(np.r_[np.arange(10,xmax+1,step,dtype=np.int64),xmax])
    engine=RealCutRk(kmax=3,Q=Q,h=h)
    approx=engine.evaluate(xs)
    counts=exact_counts(xmax,3)[xs]
    return xs,approx,counts,engine.metadata()


def save_adjacent_plots(xs,approx,counts,outdir='.',show=False):
    """Two separate unsmoothed plots. No data are reconstructed from images."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    errors=np.asarray(approx)-np.asarray(counts)
    saved=[]
    for k in (1,2):
        predicted=first_zero_scale(k,xs)*errors[:,k-1]
        target=-errors[:,k]
        fig,ax=plt.subplots(figsize=(9.0,4.0))
        ax.plot(xs,predicted,color='#2455a4',lw=0.58,label=r'$s_%d(x)D_%d(x)$ (scaled)'%(k,k))
        ax.plot(xs,target,color='#c33f3f',lw=0.58,alpha=0.8,label=r'$-D_%d(x)$ (unscaled)'%(k+1))
        ax.axhline(0,color='0.55',lw=0.45,zorder=0)
        ax.set_xlabel(r'$x$',fontsize=11)
        ax.set_ylabel('Counting-error units',fontsize=10)
        ax.set_title('Adjacent almost-prime errors: k = %d and %d'%(k,k+1),fontsize=12,pad=11)
        ax.set_xlim(float(xs[0]),float(xs[-1]))
        ax.tick_params(labelsize=9,width=0.6,length=3)
        for spine in ax.spines.values(): spine.set_linewidth(0.6); spine.set_color('0.35')
        ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0),useMathText=True)
        ax.legend(loc='upper center',bbox_to_anchor=(0.5,1.01),ncol=2,
                  frameon=False,fontsize=10,handlelength=2.7)
        fig.subplots_adjust(left=.095,right=.985,bottom=.19,top=.83)
        fig.text(.095,.035,'Finite-cut numerical illustration; lower-cut error is not certified. No smoothing.',fontsize=8,color='0.35')
        stem=out/('adjacent_k%d_k%d'%(k,k+1))
        for ext in ['png','pdf','svg']:
            path=str(stem)+'.'+ext
            fig.savefig(path,dpi=300,facecolor='white')
            saved.append(path)
        if show: plt.show()
        plt.close(fig)
    return saved


def export_demo(outdir='rk_demo',xmax=1000000,step=100,Q=128,h=0.05,show=False):
    xs,approx,counts,meta=adjacent_data(xmax,step,Q,h)
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    with open(out/'adjacent_data.csv','w',newline='') as f:
        writer=csv.writer(f); writer.writerow(['x','R1','R2_approx','R3_approx','N1','N2','N3'])
        for x,r,n in zip(xs,approx,counts): writer.writerow([int(x),*map(float,r),*map(int,n)])
    with open(out/'metadata.json','w') as f: json.dump(meta,f,indent=2)
    return save_adjacent_plots(xs,approx,counts,out,show)


''' )

report = rk_report(k=2, x=10**100, Q=64)
print('Rounded numerical approximation:', round(report['value']))
print('Approximate scientific form: {:.9e}'.format(report['value']))
print('Not certified. Large trailing integer digits are not reliable.')
