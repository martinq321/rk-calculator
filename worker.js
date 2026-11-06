/* All calculator arithmetic runs locally in this worker. No user input is
   sent to a compute service. The optional Sage page is separate. */
let state;
async function initialise(){
 const meta=await (await fetch('kernels.json')).json();
 const load=async m=>{const r=await fetch(m.file);if(!r.ok)throw Error('Kernel download failed'); const b=await r.arrayBuffer();const a=new Float64Array(b);if(a.length!==(m.kmax+1)*m.nodes)throw Error('Invalid kernel size');return a;};
 return {meta,fine:await load(meta.fine),mesh:await load(meta.mesh)};
}
function R1(x,den){let y=Math.log(x),t=1,s=1;for(let n=1;n<=den.length;n++){t*=y/n;let a=t/den[n-1];s+=a;if(n>y+10&&Math.abs(a)<2e-16*Math.abs(s))return s;}throw Error('Gram series stopping rule not reached');}
function dot(a,n,k,x,minSigma=0){let s=0,c=0,L=Math.log(x);for(let i=0;i<n;i++){if(a[i]<=minSigma)continue;let y=Math.exp(L*a[i])*a[k*n+i]-c,t=s+y;c=(t-s)-y;s=t;}return s;}
function exactNk(k,n){if(n>1000000)return null;const omega=new Uint8Array(n+1);for(let p=2;p<=n;p++){if(omega[p]!==0)continue;for(let power=p;power<=n;){for(let m=power;m<=n;m+=power)omega[m]++;if(power>Math.floor(n/p))break;power*=p;}}let count=0;for(let m=2;m<=n;m++)if(omega[m]===k)count++;return count;}
self.onmessage=async ev=>{const {id,k,x}=ev.data;try{if(!Number.isInteger(k)||k<1||k>6||!Number.isFinite(x)||x<=1||x>1e12)throw Error('Use k=1..6 and 1 < x <= 10^12.');if(!state)state=initialise();const a=await state;let value,cutoff=null,mesh=null;if(k===1)value=R1(x,a.meta.gram_denominators);else{value=dot(a.fine,a.meta.fine.nodes,k,x);cutoff=Math.abs(value-dot(a.fine,a.meta.fine.nodes,k,x,1/65));mesh=Math.abs(value-dot(a.mesh,a.meta.mesh.nodes,k,x));}self.postMessage({id,value,cutoff,mesh,count:exactNk(k,Math.floor(x)),k,x});}catch(e){state=undefined;self.postMessage({id,error:String(e.message||e)});}};
