/* HRTF Explorer — EarSim vs BTE vs P0006 (CherISH) application code.
 *
 * Data (DATA_ALL, IRD) is loaded from JSON via fetch() at startup, so this
 * page must be served over HTTP(S) — e.g. GitHub Pages or
 * `python3 -m http.server` in this folder. Opening index.html via file://
 * will fail (browsers block fetch on file://).
 */
(async () => {
const DATA_ALL = await fetch('data/data_all.json').then(r => {
  if (!r.ok) throw new Error('data_all.json: HTTP ' + r.status);
  return r.json();
});
const IRD = await fetch('data/ird.json').then(r => {
  if (!r.ok) throw new Error('ird.json: HTTP ' + r.status);
  return r.json();
});


const F = DATA_ALL.freqs, NF = F.length;
const YMAG = DATA_ALL.ranges.mag, YDIF = DATA_ALL.ranges.dif,
      YNORM = DATA_ALL.ranges.norm, LVL = DATA_ALL.ranges.lvl;
const DF = DATA_ALL.diffuse;  // global diffuse-field responses (dB, log-f grid)
const $ = id => document.getElementById(id);
const COLM = '#4da3ff', COLS = '#ff7a59', COLH = '#22d3ee';
const INK = '#d7e0ea', DIM = '#8fa3b8', GRID = 'rgba(130,160,190,0.14)';
const FONT = '-apple-system,"Segoe UI",Roboto,sans-serif';
const GAPRGB = [28, 37, 48];
function hx(h){return [parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)];}
const lerp=(a,b,t)=>a+(b-a)*t;
function cmap(stops,t){t=Math.min(1,Math.max(0,t));
  for(let i=1;i<stops.length;i++){if(t<=stops[i][0]){
    const [t0,c0]=stops[i-1],[t1,c1]=stops[i],u=(t-t0)/(t1-t0),A=hx(c0),B=hx(c1);
    return [lerp(A[0],B[0],u)|0,lerp(A[1],B[1],u)|0,lerp(A[2],B[2],u)|0];}}
  return hx(stops[stops.length-1][1]);}
function lut(stops){const L=[];for(let i=0;i<256;i++)L.push(cmap(stops,i/255));return L;}
const INFERNO=[[0,'#000004'],[.25,'#3b0f4f'],[.5,'#832681'],[.75,'#f1605d'],[1,'#fcffa4']];
const RDBU=[[0,'#2166ac'],[.5,'#f7f7f7'],[1,'#b2182b']];

/* ---------- Gamper-2013a interpolation: Delaunay simplices + barycentric weights ----------
   Mirrors the reference implementation (bela-hrir-convolver, raw-hrir-interpolation):
   precomputed (V^-1)^T per simplex (vertices as rows), first all-positive weight set wins,
   weights normalized. Blends dB magnitude spectra + scalar ITD/ILD (reference blends HRIRs).
   A hit is rejected when the query elevation is below the simplex's lowest vertex
   elevation (-1 deg): this neutralizes hull cap facets so below-rim queries are omitted
   instead of extrapolated. */
const D2R = Math.PI/180;
// Azimuth sign convention: the datasets store +azimuth = listener's LEFT
// (counterclockwise seen from above), but the UI presents +azimuth = RIGHT so
// the slider moves the way it looks. Negate once, at the UI->data boundary;
// everything data-side (unit(), buildInterp(), query(), SH basis) is untouched.
const toDataAz = a => -a;
function unit(az,el){const a=az*D2R,e=el*D2R,ce=Math.cos(e);
  return [ce*Math.cos(a),ce*Math.sin(a),Math.sin(e)];}
function inv3T(ax,ay,az,bx,by,bz,cx,cy,cz){
  // (V^-1)^T flattened row-major, V rows = a,b,c (same convention as reference)
  const c11=(by*cz-bz*cy), c12=-(bx*cz-bz*cx), c13=bx*cy-by*cx;
  const det=ax*c11+ay*c12+az*c13, id=1/det;
  return [c11*id, c12*id, c13*id,
    -(ay*cz-az*cy)*id, (ax*cz-az*cx)*id, -(ax*cy-ay*cx)*id,
    (ay*bz-az*by)*id, -(ax*bz-az*bx)*id, (ax*by-ay*bx)*id];
}
function buildInterp(ds, triKey){
  const arr=DATA_ALL[ds], n=arr.length, P=new Float64Array(n*3);
  for(let i=0;i<n;i++){const v=unit(arr[i].az,arr[i].el);
    P[i*3]=v[0];P[i*3+1]=v[1];P[i*3+2]=v[2];}
  const tris=DATA_ALL[triKey], m=tris.length;
  const Minv=new Float64Array(m*9), minEl=new Float64Array(m);
  for(let t=0;t<m;t++){const a=tris[t][0],b=tris[t][1],c=tris[t][2];
    const M=inv3T(P[a*3],P[a*3+1],P[a*3+2],P[b*3],P[b*3+1],P[b*3+2],P[c*3],P[c*3+1],P[c*3+2]);
    for(let k=0;k<9;k++)Minv[t*9+k]=M[k];
    minEl[t]=Math.min(arr[a].el,arr[b].el,arr[c].el);}
  return {arr:arr,tris:tris,Minv:Minv,minEl:minEl,hint:0};
}
const IM=buildInterp('meta','tri_m'), IS=buildInterp('soni','tri_s');
const HAS_HUMAN = Array.isArray(DATA_ALL.human) && !!IRD.human;
const IH = HAS_HUMAN ? buildInterp('human','tri_h') : null;
const LBL = {meta:'EarSim', soni:'BTE', human:'P0006'};

function queryW(I,az,el){
  const v=unit(az,el),x=v[0],y=v[1],z=v[2],m=I.tris.length;
  for(let n=0;n<m;n++){
    const t=(I.hint+n)%m,o=t*9;
    const g0=I.Minv[o]*x+I.Minv[o+1]*y+I.Minv[o+2]*z;
    if(g0<-1e-9)continue;
    const g1=I.Minv[o+3]*x+I.Minv[o+4]*y+I.Minv[o+5]*z;
    if(g1<-1e-9)continue;
    const g2=I.Minv[o+6]*x+I.Minv[o+7]*y+I.Minv[o+8]*z;
    if(g2<-1e-9)continue;
    if(el<I.minEl[t]-1.0)continue;  // below measured rim: omit, don't extrapolate
    const s=g0+g1+g2;
    I.hint=t;
    return {t:t,w:[g0/s,g1/s,g2/s]};
  }
  return null;
}
function query(I,az,el){
  const q=queryW(I,az,el);
  if(!q)return null;
  const w=q.w,tr=I.tris[q.t];
  const A=I.arr[tr[0]],B=I.arr[tr[1]],C=I.arr[tr[2]];
  const L=new Array(NF),R=new Array(NF);
  for(let i=0;i<NF;i++){L[i]=w[0]*A.L[i]+w[1]*B.L[i]+w[2]*C.L[i];
    R[i]=w[0]*A.R[i]+w[1]*B.R[i]+w[2]*C.R[i];}
  return {L:L,R:R,itd:w[0]*A.itd+w[1]*B.itd+w[2]*C.itd,ild:w[0]*A.ild+w[1]*B.ild+w[2]*C.ild};
}

/* ---------- time-domain HRIRs for the binaural preview player ---------- */

const IRF={};
for(const ds of Object.keys(IRD)){
  if(IRD[ds] && IRD[ds].ncoeff)continue;  // SH coefficient bank (unused in this build)
  const e=IRD[ds],n=e.npos*e.nsamp;
  const raw=atob(e.b64),bytes=new Uint8Array(raw.length);
  for(let i=0;i<raw.length;i++)bytes[i]=raw.charCodeAt(i);
  const q=new Int16Array(bytes.buffer),L=new Float32Array(n),R=new Float32Array(n);
  for(let p=0;p<e.npos;p++)for(let s=0;s<e.nsamp;s++){
    L[p*e.nsamp+s]=q[(p*2)*e.nsamp+s]*e.scale/32767;
    R[p*e.nsamp+s]=q[(p*2+1)*e.nsamp+s]*e.scale/32767;
  }
  IRF[ds]={L:L,R:R,nsamp:e.nsamp,npos:e.npos,itd:Float32Array.from(e.itd)};
}
// Fractional delay on one ear: integer part shifts taps, fractional part
// linearly interpolates between adjacent taps.
// Follows yoyolicoris/bela-hrir-convolver branch minimum-phase-itd (render.cpp).
function delayEar(h,d){
  if(d<1e-12)return h.slice();
  const delay=Math.floor(d),p=d-delay,ns=h.length,g=new Float32Array(ns+delay+1);
  g[delay]=h[0]*(1-p);
  for(let k=1;k<ns;k++)g[delay+k]=h[k]*(1-p)+h[k-1]*p;
  g[delay+ns]=h[ns-1]*p;
  return g;
}
function padTo(h,n){ // zero-pad to length n (no-op when already long enough)
  if(h.length>=n)return h;
  const g=new Float32Array(n);g.set(h);return g;
}
function queryIR(ds,az,el){
  // Barycentric interpolation on HRIRs and on the fractional
  // ITD separately; the ITD (+ve = right lags) is applied to the lagging ear.
  const I=ds==='meta'?IM:ds==='human'?IH:IS,F=IRF[ds];
  if(!I||!F)return null;  // dataset's time-domain IRs not shipped in this build
  const q=queryW(I,az,el);
  if(!q)return null;
  const w=q.w,tr=I.tris[q.t],ns=F.nsamp;
  const L=new Float32Array(ns),R=new Float32Array(ns);
  const i0=tr[0]*ns,i1=tr[1]*ns,i2=tr[2]*ns;
  for(let s=0;s<ns;s++){
    L[s]=w[0]*F.L[i0+s]+w[1]*F.L[i1+s]+w[2]*F.L[i2+s];
    R[s]=w[0]*F.R[i0+s]+w[1]*F.R[i1+s]+w[2]*F.R[i2+s];
  }
  let itd=w[0]*F.itd[tr[0]]+w[1]*F.itd[tr[1]]+w[2]*F.itd[tr[2]];
  /* no ITD bias compensation for EarSim/BTE KEMAR data */
  let oL=L,oR=R,n=ns;
  if(itd>1e-12){oR=delayEar(R,itd);n=oR.length;oL=padTo(L,n);}
  else if(itd<-1e-12){oL=delayEar(L,-itd);n=oL.length;oR=padTo(R,n);}
  return {L:oL,R:oR,ns:n}; // causal: delayed tail kept, other ear zero-padded to match
}

/* ---------- direct SH resynthesis (— dataset; no interpolation) ---------- */
// Real spherical harmonics, orthonormal, Condon-Shortley phase, ACN ordering:
//   Y_n^m = N * P_n^|m|(cos th) * sqrt2*cos(m*phi) (m>0), *1 (m=0), sqrt2*sin(|m|*phi) (m<0)
//   th = colatitude from +z = 90deg - elevation, phi = azimuth from +x toward +y
// Coefficients are the order-15 complex least-squares SH fit to the provided
// HRTF spectra (mapped to time domain for shipping); see tools/build_sh.py.
// Synthesis evaluates the basis at the exact requested direction, so there is
// no spatial interpolation.
const HAS_SH = false; // SH resynthesis disabled; third dataset is P0006 (human)
let SHC = null;
if (HAS_SH) {
  const e = IRD.meta_sh, n = e.ncoeff * e.nsamp;
  const raw = atob(e.b64), bytes = new Uint8Array(raw.length);
  for (let i = 0; i < raw.length; i++) bytes[i] = raw.charCodeAt(i);
  const q = new Int16Array(bytes.buffer);
  const L = new Float32Array(n), R = new Float32Array(n);
  for (let i = 0; i < n; i++) { L[i] = q[i] * e.scale[0] / 32767; R[i] = q[n + i] * e.scale[1] / 32767; }
  SHC = { L: L, R: R, K: e.ncoeff, N: e.nsamp, order: e.sh_order, diffuse: e.diffuse };
}
function lgamma(x) { // Lanczos approximation, x > 0
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
  x -= 1;
  let a = c[0];
  for (let i = 1; i < 9; i++) a += c[i] / (x + i);
  const t = x + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}
const SH_NORM = (() => { // orthonormalization N_n^|m|, indexed [n*(order+1)+|m|]
  if (!SHC) return null;
  const W = SHC.order + 1, T = new Float64Array(W * W);
  for (let n = 0; n <= SHC.order; n++) for (let m = 0; m <= n; m++)
    T[n * W + m] = Math.exp(0.5 * (Math.log(2 * n + 1) - Math.log(4 * Math.PI)
      + lgamma(n - m + 1) - lgamma(n + m + 1)));
  return T;
})();
function shBasis(az, el) { // -> Float64Array(K): real SH at one direction
  const order = SHC.order, K = SHC.K, W = order + 1;
  const th = (90 - el) * D2R, phi = az * D2R, ct = Math.cos(th);
  const u = Math.sqrt(Math.max(0, 1 - ct * ct));
  const P = new Float64Array(W * W);
  P[0] = 1;
  let df = 1, um = 1;
  for (let m = 1; m <= order; m++) { df *= 2 * m - 1; um *= u; P[m * W + m] = ((m % 2) ? -1 : 1) * df * um; }
  for (let m = 0; m < order; m++) P[(m + 1) * W + m] = ct * (2 * m + 1) * P[m * W + m];
  for (let m = 0; m <= order; m++) for (let n = m + 2; n <= order; n++)
    P[n * W + m] = (ct * (2 * n - 1) * P[(n - 1) * W + m] - (n + m - 1) * P[(n - 2) * W + m]) / (n - m);
  const Y = new Float64Array(K);
  for (let n = 0; n <= order; n++) for (let m = -n; m <= n; m++) {
    const am = Math.abs(m), base = SH_NORM[n * W + am] * P[n * W + am];
    const k = n * n + n + m;
    if (m > 0) Y[k] = Math.SQRT2 * base * Math.cos(m * phi);
    else if (m < 0) Y[k] = Math.SQRT2 * base * Math.sin(am * phi);
    else Y[k] = base;
  }
  return Y;
}
function synthSH(az, el) { // -> {L,R}: direct HRIR resynthesis, 384 samples each
  const Y = shBasis(az, el), K = SHC.K, N = SHC.N;
  const L = new Float32Array(N), R = new Float32Array(N);
  for (let k = 0; k < K; k++) {
    const y = Y[k]; if (y === 0) continue;
    const o = k * N;
    for (let s = 0; s < N; s++) { L[s] += y * SHC.L[o + s]; R[s] += y * SHC.R[o + s]; }
  }
  return { L: L, R: R };
}
function fft2048(re, im) { // in-place radix-2, n = 2048
  const n = 2048;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = -2 * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cwr = 1, cwi = 0;
      for (let j = 0; j < len / 2; j++) {
        const a = i + j, b = i + j + len / 2;
        const vr = re[b] * cwr - im[b] * cwi, vi = re[b] * cwi + im[b] * cwr;
        re[b] = re[a] - vr; im[b] = im[a] - vi;
        re[a] += vr; im[a] += vi;
        const nwr = cwr * wr - cwi * wi; cwi = cwr * wi + cwi * wr; cwr = nwr;
      }
    }
  }
}
function specDB(h) { // 384-sample HRIR -> dB magnitude on the shared log-f grid F
  const re = new Float64Array(2048), im = new Float64Array(2048);
  for (let s = 0; s < h.length; s++) re[s] = h[s];
  fft2048(re, im);
  const m = new Float64Array(1025);
  for (let i = 0; i < 1025; i++) m[i] = 20 * Math.log10(Math.hypot(re[i], im[i]) + 1e-12);
  const out = new Array(NF);
  let j = 1;
  for (let i = 0; i < NF; i++) {
    const lf = Math.log(F[i]);
    while (j < 1023 && Math.log(48000 * (j + 1) / 2048) < lf) j++;
    const l0 = Math.log(48000 * j / 2048), l1 = Math.log(48000 * (j + 1) / 2048);
    const t = (lf - l0) / (l1 - l0);
    out[i] = m[j] * (1 - t) + m[j + 1] * t;
  }
  return out;
}
// 4th-order Butterworth 2 kHz lowpass (matches the Python ITD estimator)
const LP_B = [0.00021313872697507864, 0.0008525549079003145, 0.0012788323618504718, 0.0008525549079003145, 0.00021313872697507864];
const LP_A = [1.0, -3.3168079106244175, 4.174245550076573, -2.3574027805622575, 0.5033753607417041];
function lfilter(b, a, x) {
  const y = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) {
    let s = b[0] * x[i];
    for (let k = 1; k < b.length; k++) s += b[k] * (i - k >= 0 ? x[i - k] : 0);
    for (let k = 1; k < a.length; k++) s -= a[k] * (i - k >= 0 ? y[i - k] : 0);
    y[i] = s;
  }
  return y;
}
function filtfilt(b, a, x) { // zero-phase, odd-extension padding like scipy
  const n = x.length, pe = 12;
  const xp = new Float64Array(n + 2 * pe);
  for (let i = 0; i < pe; i++) { xp[i] = 2 * x[0] - x[pe - i]; xp[n + pe + i] = 2 * x[n - 1] - x[n - 2 - i]; }
  for (let i = 0; i < n; i++) xp[pe + i] = x[i];
  let y = lfilter(b, a, xp);
  y = lfilter(b, a, y.reverse()).reverse();
  return y.slice(pe, pe + n);
}
function xcorrLag(rf, lf) { // integer lag maximizing <rf(t+lag), lf(t)>; +ve = rf lags
  const n = rf.length;
  let best = -1e18, bi = 0;
  for (let lag = -(n - 1); lag < n; lag++) {
    let s = 0;
    const i0 = Math.max(0, -lag), i1 = Math.min(n, n - lag);
    for (let i = i0; i < i1; i++) s += rf[i + lag] * lf[i];
    if (s > best) { best = s; bi = lag; }
  }
  return bi;
}
function cuesFromIR(L, R) { // ITD (us, +ve = right lags) + ILD (dB), same estimator as tab 1
  const lf = filtfilt(LP_B, LP_A, L), rf = filtfilt(LP_B, LP_A, R);
  const itd = xcorrLag(rf, lf) / 48000 * 1e6;
  let sl = 0, sr = 0;
  for (let s = 0; s < L.length; s++) { sl += L[s] * L[s]; sr += R[s] * R[s]; }
  const ild = 20 * Math.log10(Math.sqrt(sr / L.length) / Math.sqrt(sl / L.length) + 1e-12);
  return { itd: itd, ild: ild };
}
function shSpec(az, el) { // synthesized spectra + cues at one direction (no interpolation)
  const ir = synthSH(az, el);
  const c = cuesFromIR(ir.L, ir.R);
  return { L: specDB(ir.L), R: specDB(ir.R), itd: c.itd, ild: c.ild };
}
let shSweepCache = { key: null, irs: null };
function getSHSweep() { // 181 synthesized IR pairs along the current sweep; cached
  const key = sType + '|' + (sType === 'az' ? sFixEl : sFixAz);
  if (shSweepCache.key === key && shSweepCache.irs) return shSweepCache.irs;
  const azSw = sType === 'az', M = 181;
  const order = SHC.order, K = SHC.K, N = SHC.N, W = order + 1;
  // SH basis at all M sweep directions (M x K), then a single batched
  // matmul H = Y @ C per ear instead of a per-direction synthesis loop.
  const Ym = new Float64Array(M * K);
  if (azSw) {
    // Azimuth sweep: elevation fixed, so Legendre P_n^m(cos(el)) and the
    // normalization are identical for all directions — compute once.
    const th = (90 - sFixEl) * D2R, ct = Math.cos(th);
    const u = Math.sqrt(Math.max(0, 1 - ct * ct));
    const P = new Float64Array(W * W);
    P[0] = 1;
    let df = 1, um = 1;
    for (let m = 1; m <= order; m++) { df *= 2 * m - 1; um *= u; P[m * W + m] = ((m % 2) ? -1 : 1) * df * um; }
    for (let m = 0; m < order; m++) P[(m + 1) * W + m] = ct * (2 * m + 1) * P[m * W + m];
    for (let m = 0; m <= order; m++) for (let n = m + 2; n <= order; n++)
      P[n * W + m] = (ct * (2 * n - 1) * P[(n - 1) * W + m] - (n + m - 1) * P[(n - 2) * W + m]) / (n - m);
    const base = new Float64Array(K);
    for (let n = 0; n <= order; n++) for (let m = -n; m <= n; m++)
      base[n * n + n + m] = SH_NORM[n * W + Math.abs(m)] * P[n * W + Math.abs(m)];
    for (let i = 0; i < M; i++) {
      const phi = toDataAz(-180 + i * 2) * D2R, row = i * K;
      for (let n = 0; n <= order; n++) for (let m = -n; m <= n; m++) {
        const k = n * n + n + m, b = base[k];
        let y;
        if (m > 0) y = Math.SQRT2 * b * Math.cos(m * phi);
        else if (m < 0) y = Math.SQRT2 * b * Math.sin(-m * phi);
        else y = b;
        Ym[row + k] = y;
      }
    }
  } else {
    // Elevation sweep: azimuth fixed, so the trig factors are constant —
    // precompute them, then only the Legendre recurrence varies per direction.
    const phi = toDataAz(sFixAz) * D2R;
    const cosM = new Float64Array(W), sinM = new Float64Array(W);
    for (let m = 0; m <= order; m++) { cosM[m] = Math.cos(m * phi); sinM[m] = Math.sin(m * phi); }
    for (let i = 0; i < M; i++) {
      const el = -90 + i;
      const th = (90 - el) * D2R, ct = Math.cos(th);
      const u = Math.sqrt(Math.max(0, 1 - ct * ct));
      const P = new Float64Array(W * W);
      P[0] = 1;
      let df = 1, um = 1;
      for (let m = 1; m <= order; m++) { df *= 2 * m - 1; um *= u; P[m * W + m] = ((m % 2) ? -1 : 1) * df * um; }
      for (let m = 0; m < order; m++) P[(m + 1) * W + m] = ct * (2 * m + 1) * P[m * W + m];
      for (let m = 0; m <= order; m++) for (let n = m + 2; n <= order; n++)
        P[n * W + m] = (ct * (2 * n - 1) * P[(n - 1) * W + m] - (n + m - 1) * P[(n - 2) * W + m]) / (n - m);
      const row = i * K;
      for (let n = 0; n <= order; n++) for (let m = -n; m <= n; m++) {
        const am = Math.abs(m), b = SH_NORM[n * W + am] * P[n * W + am];
        const k = n * n + n + m;
        let y;
        if (m > 0) y = Math.SQRT2 * b * cosM[m];
        else if (m < 0) y = Math.SQRT2 * b * sinM[am];
        else y = b;
        Ym[row + k] = y;
      }
    }
  }
  // Batched synthesis H = Ym @ C for both ears (k-major: each coefficient
  // row of C is reused across all directions while hot in cache).
  const Lf = new Float32Array(M * N), Rf = new Float32Array(M * N);
  for (let k = 0; k < K; k++) {
    const cOff = k * N;
    for (let i = 0; i < M; i++) {
      const yk = Ym[i * K + k];
      if (yk === 0) continue;
      const hOff = i * N;
      for (let s = 0; s < N; s++) { Lf[hOff + s] += yk * SHC.L[cOff + s]; Rf[hOff + s] += yk * SHC.R[cOff + s]; }
    }
  }
  const L = [], R = [];
  for (let i = 0; i < M; i++) { L.push(Lf.subarray(i * N, (i + 1) * N)); R.push(Rf.subarray(i * N, (i + 1) * N)); }
  shSweepCache = { key: key, irs: { L: L, R: R } };
  return shSweepCache.irs;
}

/* ---------- shared ---------- */
function seg(sel,b){document.querySelectorAll(sel+' button').forEach(x=>x.classList.remove('on'));b.classList.add('on');}
$('tab1').onclick=()=>setTab(1); $('tab2').onclick=()=>setTab(2);
function setTab(n){
  $('tab1').classList.toggle('on',n===1); $('tab2').classList.toggle('on',n===2);
  $('panel1').classList.toggle('hidden',n!==1); $('panel2').classList.toggle('hidden',n!==2);
  $('side1').classList.toggle('hidden',n!==1); $('side2').classList.toggle('hidden',n!==2);
  if(n===1)refresh1(); else refreshSlice();
}
let raf1=false, rafS=false;
function schedule1(){if(!raf1){raf1=true;requestAnimationFrame(()=>{raf1=false;refresh1();});}}
function scheduleS(){if(!rafS){rafS=true;requestAnimationFrame(()=>{rafS=false;refreshSlice();});}}

/* ---------- panel 1: direction compare (continuous) ---------- */
let mode='overlay', fscale='log';
const plotL=$('plotL'), plotR=$('plotR');
function cur1(){const uaz=+$('az').value, el=+$('el').value, az=toDataAz(uaz);
  const o={az:uaz, el:el, m:query(IM,az,el), s:query(IS,az,el)};
  if(HAS_HUMAN&&$('c_h').checked)o.h=query(IH,az,el);
  return o;}
function aligned1(){
  const q=cur1(), a=$('align').checked, m=q.m, s=q.s, h=q.h;
  const off=e=>{if(!a||!m||!s)return 0;let sm=0,n=0;
    for(let i=0;i<NF;i++)if(F[i]>=200&&F[i]<=16000){sm+=m[e][i]-s[e][i];n++;}return sm/n;};
  const offH=e=>{if(!a||!m||!h)return 0;let sm=0,n=0;
    for(let i=0;i<NF;i++)if(F[i]>=200&&F[i]<=16000){sm+=m[e][i]-h[e][i];n++;}return sm/n;};
  const oL=off('L'),oR=off('R'), hL=offH('L'), hR=offH('R');
  return {mL:m&&m.L, mR:m&&m.R,
    sL:s&&s.L.map(v=>v+oL), sR:s&&s.R.map(v=>v+oR), m:m, s:s,
    hL:h&&h.L.map(v=>v+hL), hR:h&&h.R.map(v=>v+hR), h:h};}
function xPos(f,W,L,R){
  if(fscale==='log'){const l0=Math.log(100),l1=Math.log(20000);
    return L+(Math.log(f)-l0)/(l1-l0)*(R-L);}
  return L+(f-100)/(20000-100)*(R-L);}
function drawEar(cv,ear,C,ymin,ymax,note){
  const dpr=window.devicePixelRatio||1,W=cv.clientWidth,H=cv.clientHeight;
  if(!W||!H)return;
  cv.width=W*dpr;cv.height=H*dpr;
  const c=cv.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);
  c.fillStyle='#0b1016';c.fillRect(0,0,W,H);
  const L=56,R=W-14,T=30,B=36;
  const Y=v=>T+(1-(v-ymin)/(ymax-ymin))*(H-T-B);
  c.font='11px '+FONT;c.lineWidth=1;
  const ticks=fscale==='log'?[100,200,500,1000,2000,5000,10000,20000]:[2000,4000,6000,8000,10000,12000,14000,16000,18000,20000];
  ticks.forEach(f=>{const x=xPos(f,W,L,R);
    c.strokeStyle=GRID;c.beginPath();c.moveTo(x,T);c.lineTo(x,H-B);c.stroke();
    c.fillStyle=DIM;c.textAlign='center';c.fillText(f>=1000?(f/1000)+'k':f,x,H-B+16);});
  const step=10;
  c.textAlign='right';
  for(let v=Math.ceil(ymin/step)*step;v<=ymax;v+=step){const y=Y(v);
    c.strokeStyle=GRID;c.beginPath();c.moveTo(L,y);c.lineTo(R,y);c.stroke();
    c.fillStyle=DIM;c.fillText(v+'',L-8,y+4);}
  c.textAlign='left';c.fillStyle=DIM;
  c.fillText(mode==='overlay'?'magnitude (dB)':'EarSim \u2212 BTE (dB)',8,T-10);
  c.fillStyle=INK;c.font='600 13px '+FONT;
  c.fillText(ear==='L'?'Left ear':'Right ear',L,17);
  if(mode==='diff'&&ymin<0&&ymax>0){c.strokeStyle='rgba(215,224,234,0.35)';c.setLineDash([5,4]);
    c.beginPath();c.moveTo(L,Y(0));c.lineTo(R,Y(0));c.stroke();c.setLineDash([]);}
  const showM=$('c_m').checked, showS=$('c_s').checked, showH=HAS_HUMAN&&$('c_h').checked;
  const series=[];
  if(mode==='overlay'){
    if(showM&&C.mL)series.push(['EarSim',COLM,C['m'+ear]]);
    if(showH&&C.hL)series.push([LBL.human,COLH,C['h'+ear]]);
    if(showS&&C.sL)series.push(['BTE',COLS,C['s'+ear]]);
  }else{
    if(C.mL&&C.sL)series.push(['\u0394 EarSim\u2212BTE','#ffd166',
      F.map((_,i)=>C['m'+ear][i]-C['s'+ear][i])]);
    if(C.hL&&C.mL&&showH)series.push(['\u0394 P0006\u2212EarSim',COLH,
      F.map((_,i)=>C['h'+ear][i]-C['m'+ear][i])]);
  }
  if(!series.length){
    c.fillStyle=DIM;c.font='13px '+FONT;c.textAlign='center';
    c.fillText(note||'No data at this direction',(L+R)/2,(T+H-B)/2);
    c.textAlign='left';
  }
  series.forEach(s=>{const name=s[0],col=s[1],vals=s[2];
    c.strokeStyle=col;c.lineWidth=1.8;c.beginPath();
    vals.forEach((v,i)=>{const x=xPos(F[i],W,L,R),y=Y(v);i?c.lineTo(x,y):c.moveTo(x,y);});
    c.stroke();});
  c.font='12px '+FONT;
  let lx=R-8;
  [...series].reverse().forEach(s=>{const name=s[0],col=s[1];
    const w=c.measureText(name).width;lx-=w+26;
    c.strokeStyle=col;c.lineWidth=2.5;c.beginPath();
    c.moveTo(lx,T-12);c.lineTo(lx+16,T-12);c.stroke();
    c.fillStyle=INK;c.fillText(name,lx+20,T-8);lx-=6;});
}
function draw1(){
  const C=aligned1();
  const ymin=mode==='overlay'?YMAG[0]:-YDIF, ymax=mode==='overlay'?YMAG[1]:YDIF;
  let note='';
  if(mode==='overlay'&&!C.mL&&!C.sL)note='Outside both measurement grids';
  else if(mode==='overlay'&&!C.mL)note='EarSim: outside measurement grid';
  else if(mode==='overlay'&&!C.sL)note='BTE: outside measurement grid';
  else if(mode==='diff'&&!(C.mL&&C.sL))note='Difference needs both datasets here';
  drawEar(plotL,'L',C,ymin,ymax,note);
  drawEar(plotR,'R',C,ymin,ymax,note);
}
function stats1(){
  const q=cur1(),C=aligned1(),m=q.m,s=q.s,h=q.h;
  $('azv').textContent=q.az+'\u00B0';$('elv').textContent=q.el+'\u00B0';
  const cov=[];
  if(!m)cov.push('EarSim: outside grid');
  if(!s)cov.push('BTE: outside grid');
  $('dirinfo').innerHTML='az '+q.az+'\u00B0, el '+q.el+'\u00B0'+
    (cov.length?'<br><span class="warn">'+cov.join(' \u00B7 ')+'</span>':'<br><span class="hint">barycentric interpolation</span>');
  $('itd_m').textContent=m?(m.itd+(compITD?ITD_COMP:0)).toFixed(1):'\u2014';
  $('itd_s').textContent=s?s.itd.toFixed(1):'\u2014';
  $('itd_h').textContent=h?h.itd.toFixed(1):'\u2014';
  $('ild_m').textContent=m?m.ild.toFixed(2):'\u2014';
  $('ild_s').textContent=s?s.ild.toFixed(2):'\u2014';
  $('ild_h').textContent=h?h.ild.toFixed(2):'\u2014';
  const idx=[];for(let i=0;i<NF;i++)if(F[i]>=200&&F[i]<=16000)idx.push(i);
  if(C.mL&&C.sL){
    const dL=idx.map(i=>Math.abs(C.mL[i]-C.sL[i])),dR=idx.map(i=>Math.abs(C.mR[i]-C.sR[i]));
    const med=a=>{const t=[...a].sort((x,y)=>x-y);return t[t.length>>1];};
    let mx=0,mf=0;
    idx.forEach(i=>{const d=Math.max(Math.abs(C.mL[i]-C.sL[i]),Math.abs(C.mR[i]-C.sR[i]));
      if(d>mx){mx=d;mf=F[i];}});
    $('dl').textContent=med(dL).toFixed(2)+' dB';
    $('dr').textContent=med(dR).toFixed(2)+' dB';
    $('dmax').textContent=mx.toFixed(1)+' dB @ '+(mf>=1000?(mf/1000).toFixed(1)+'k':Math.round(mf));
  }else{
    $('dl').textContent=$('dr').textContent=$('dmax').textContent='\u2014';
  }
}
function refresh1(){draw1();stats1();}
['az','el','c_m','c_h','c_s','align'].forEach(id=>$(id).addEventListener('input',schedule1));
$('az').addEventListener('input',()=>{if(playing)pushIR();});
$('el').addEventListener('input',()=>{if(playing)pushIR();});

/* ---------- binaural preview player (Direction tab) ----------
   Sample-by-sample rendering in an AudioWorklet: the worklet convolves the
   mono source with the current stereo HRIR and, when the direction changes,
   crossfades the filter coefficients sample-by-sample (smoothstep, 512
   samples ~ 11 ms) to the new HRIR. No convolver buffer swaps, no clicks. */
const WORKLET_SRC = `
class BinauralFIR extends AudioWorkletProcessor{
constructor(){
super();
this.M=1024;
this.dl=new Float32Array(this.M);
this.pos=0;
this.hL=new Float32Array([1]);
this.hR=new Float32Array([1]);
this.xN=512;this.xLeft=0;
this.oL=null;this.oR=null;this.nL=null;this.nR=null;
this.port.onmessage=(e)=>{
const d=e.data;
if(!d||!d.L||!d.R)return;
const n=Math.max(this.hL.length,d.L.length);
const oL=new Float32Array(n);oL.set(this.hL);
const oR=new Float32Array(n);oR.set(this.hR);
const nL=new Float32Array(n);nL.set(d.L);
const nR=new Float32Array(n);nR.set(d.R);
this.oL=oL;this.oR=oR;this.nL=nL;this.nR=nR;
this.xLeft=this.xN;
};
}
process(inputs,outputs){
const ch=(inputs[0]&&inputs[0][0])||null;
const oL=outputs[0][0],oR=outputs[0][1];
const N=oL.length,M=this.M,dl=this.dl;
let pos=this.pos,hL=this.hL,hR=this.hR,xLeft=this.xLeft;
const xN=this.xN;
let oLc=this.oL,oRc=this.oR,nLc=this.nL,nRc=this.nR;
for(let n=0;n<N;n++){
const x=ch?ch[n]:0;
dl[pos]=x;
let yL=0,yR=0;
if(xLeft>0){
const t=1-xLeft/xN,a=t*t*(3-2*t),K=oLc.length;
let idx=pos;
for(let k=0;k<K;k++){
const xv=dl[idx];
yL+=(oLc[k]+a*(nLc[k]-oLc[k]))*xv;
yR+=(oRc[k]+a*(nRc[k]-oRc[k]))*xv;
if(--idx<0)idx+=M;
}
if(--xLeft===0){hL=nLc;hR=nRc;oLc=oRc=nLc=nRc=null;}
}else{
const K=hL.length;
let idx=pos;
for(let k=0;k<K;k++){
const xv=dl[idx];
yL+=hL[k]*xv;yR+=hR[k]*xv;
if(--idx<0)idx+=M;
}
}
oL[n]=yL;oR[n]=yR;
if(++pos>=M)pos=0;
}
this.pos=pos;this.hL=hL;this.hR=hR;this.xLeft=xLeft;
this.oL=oLc;this.oR=oRc;this.nL=nLc;this.nR=nRc;
return true;
}
}
registerProcessor('binaural-fir',BinauralFIR);
`;
let AC=null,binaural=null,srcGain=null,pVol=null,pSrc=null,playing=false,pDs='meta',fileBuf=null;
async function ensureAC(){
  if(AC)return;
  AC=new (window.AudioContext||window.webkitAudioContext)({sampleRate:48000});
  await AC.audioWorklet.addModule(URL.createObjectURL(new Blob([WORKLET_SRC],{type:'application/javascript'})));
  binaural=new AudioWorkletNode(AC,'binaural-fir',{numberOfInputs:1,numberOfOutputs:1,outputChannelCount:[2]});
  srcGain=AC.createGain();
  pVol=AC.createGain();
  pVol.gain.value=($('pvol').value/100);
  srcGain.connect(binaural);binaural.connect(pVol);pVol.connect(AC.destination);
}
function startSource(sbuf){
  if(pSrc){try{pSrc.stop();}catch(e){}try{pSrc.disconnect();}catch(e){}pSrc=null;}
  pSrc=AC.createBufferSource();pSrc.buffer=sbuf;pSrc.loop=true;
  pSrc.connect(srcGain);pSrc.start();
}
function swapSource(sbuf){
  // fade out, swap buffer, fade back in: no click on sound change
  const t=AC.currentTime;
  srcGain.gain.cancelScheduledValues(t);
  srcGain.gain.setValueAtTime(srcGain.gain.value,t);
  srcGain.gain.linearRampToValueAtTime(0.0001,t+0.015);
  setTimeout(()=>{
    startSource(sbuf);
    const t2=AC.currentTime;
    srcGain.gain.cancelScheduledValues(t2);
    srcGain.gain.setValueAtTime(0.0001,t2);
    srcGain.gain.linearRampToValueAtTime(1.0,t2+0.03);
  },20);
}
function normPeak(d,peak){
  let m=0;
  for(let i=0;i<d.length;i++){const a=Math.abs(d[i]);if(a>m)m=a;}
  if(m>1e-9){const g=peak/m;for(let i=0;i<d.length;i++)d[i]*=g;}
}
function synthSource(kind){
  const sr=48000,dur=kind==='sweep'?4:2,n=Math.floor(sr*dur);
  const buf=AC.createBuffer(1,n,sr),d=buf.getChannelData(0);
  if(kind==='pink'){
    let b0=0,b1=0,b2=0,b3=0,b4=0,b5=0,b6=0;
    for(let i=0;i<n;i++){const w=Math.random()*2-1;
      b0=0.99886*b0+w*0.0555179;b1=0.99332*b1+w*0.0750759;b2=0.96900*b2+w*0.1538520;
      b3=0.86650*b3+w*0.3104856;b4=0.55000*b4+w*0.5329522;b5=-0.7616*b5-w*0.0168980;
      d[i]=(b0+b1+b2+b3+b4+b5+b6+w*0.5362)*0.11;b6=w*0.115926;}
  }else if(kind==='white'){for(let i=0;i<n;i++)d[i]=(Math.random()*2-1)*0.5;}
  else if(kind==='click'){for(let i=0;i<n;i++)d[i]=0;
    for(let k=0;k<8;k++){const o=k*Math.floor(sr/4);
      for(let j=0;j<40&&o+j<n;j++)d[o+j]=(1-j/40)*(j%2?0.6:-0.6);}}
  else if(kind==='sweep'){let ph=0;
    for(let i=0;i<n;i++){const t=i/sr,f=100*Math.pow(160,t/dur);
      ph+=2*Math.PI*f/sr;
      d[i]=Math.sin(ph)*0.5*Math.min(1,i/(sr*0.05))*Math.min(1,(n-i)/(sr*0.05));}}
  normPeak(d,0.9);
  return buf;
}
function playerIR(){
  const q=queryIR(pDs,toDataAz(+$('az').value),+$('el').value);
  if(!q)return null;
  let db=null;
  if($('plmatch').checked){
    if(pDs==='soni')db=DATA_ALL.ranges.lvl;
    else if(pDs==='human'&&DATA_ALL.ranges.lvl_h!=null)db=DATA_ALL.ranges.lvl_h;
  }
  if(db!=null){
    const g=Math.pow(10,db/20);
    const L=new Float32Array(q.ns),R=new Float32Array(q.ns);
    L.set(q.L);R.set(q.R);
    for(let s=0;s<q.ns;s++){L[s]*=g;R[s]*=g;}
    return {L:L,R:R,ns:q.ns};
  }
  return q;
}
let pushScheduled=false;
function pushIR(){
  if(!AC||!playing||!binaural)return;
  if(pushScheduled)return;  // at most one update per animation frame
  pushScheduled=true;
  requestAnimationFrame(()=>{
    pushScheduled=false;
    if(!playing||!binaural)return;
    const q=playerIR();if(!q)return;  // below rim: hold last IR
    binaural.port.postMessage({L:q.L,R:q.R});
  });
}
function stopPlay(){
  if(pSrc){try{pSrc.stop();}catch(e){}try{pSrc.disconnect();}catch(e){}pSrc=null;}
  playing=false;$('pplay').innerHTML='&#9654; Play';
}
$('pplay').onclick=async()=>{
  await ensureAC();
  if(AC.state==='suspended')AC.resume();
  if(playing){stopPlay();return;}
  const kind=$('psig').value;
  let sbuf=null;
  if(kind==='file'){
    if(!fileBuf){$('pfile').click();return;}
    sbuf=fileBuf;
  }else sbuf=synthSource(kind);
  const q=playerIR();
  if(!q){$('phint').textContent='No measured data at this direction.';return;}
  binaural.port.postMessage({L:q.L,R:q.R});
  srcGain.gain.value=1.0;
  startSource(sbuf);
  playing=true;$('pplay').innerHTML='&#9632; Stop';
  $('phint').textContent=describePlaying();
};
function describePlaying(){
  return 'Playing '+(LBL[pDs]||pDs)+' IR at az '+$('az').value+'\u00B0, el '+$('el').value+'\u00B0. Drag sliders to move.';
}
$('pvol').oninput=e=>{
  $('pvolval').textContent=e.target.value+'%';
  if(pVol)pVol.gain.value=e.target.value/100;
};
$('psig').onchange=e=>{
  const kind=e.target.value;
  if(kind==='file'&&!fileBuf){$('pfile').click();return;}
  if(playing&&AC){
    swapSource(kind==='file'?fileBuf:synthSource(kind));
    $('phint').textContent=describePlaying();
  }
};
$('pfile').onchange=e=>{
  const f=e.target.files[0];if(!f)return;
  ensureAC();
  f.arrayBuffer().then(ab=>AC.decodeAudioData(ab)).then(db=>{
    const n=db.length,mono=AC.createBuffer(1,n,db.sampleRate),d=mono.getChannelData(0);
    for(let c=0;c<db.numberOfChannels;c++){const cd=db.getChannelData(c);
      for(let i=0;i<n;i++)d[i]+=cd[i]/db.numberOfChannels;}
    normPeak(d,0.9);
    fileBuf=mono;$('phint').textContent='Loaded: '+f.name+'. Press Play.';
  }).catch(()=>{$('phint').textContent='Could not decode that file.';});
};
$('psrc').querySelectorAll('button').forEach(b=>b.onclick=()=>{seg('#psrc',b);pDs=b.dataset.d;pushIR();});
$('plmatch').onchange=()=>{if(playing)pushIR();};
document.querySelectorAll('#mode button').forEach(b=>b.onclick=()=>{
  seg('#mode',b);mode=b.dataset.m;refresh1();});
document.querySelectorAll('#fscale button').forEach(b=>b.onclick=()=>{
  seg('#fscale',b);fscale=b.dataset.f;refresh1();});

/* ---------- panel 2: slice viewer (continuous) ---------- */
const sliceL=$('sliceL'), sliceR=$('sliceR');
let sDs='meta', sType='az', sNorm=false, sFscale='log';
let sFixEl=0, sFixAz=0;  // fixed angles persist across dataset / sweep switches
const LIN_N=400, LIN_F=[];
for(let li=0;li<LIN_N;li++)LIN_F.push(100+(20000-100)*li/(LIN_N-1));
function toLinear(cols){
  return cols.map(col=>{if(!col)return null;const out=new Array(LIN_N);let j=0;
    for(let i=0;i<LIN_N;i++){const f=LIN_F[i];
      while(j<NF-2&&F[j+1]<f)j++;
      const f0=F[j],f1=F[j+1],t=(f-f0)/(f1-f0);
      out[i]=col[j]*(1-t)+col[j+1]*t;}
    return out;});}
let sColsL=[], sColsR=[], sAngles=[], sVmin=0, sVmax=0, sGeoL=null, sGeoR=null, sGaps=0;
const LUT_I=lut(INFERNO), LUT_D=lut(RDBU);
function syncFixSlider(){
  const azSw=sType==='az',f=$('fix');
  $('fixlabel').textContent=azSw?'Elevation':'Azimuth';
  f.min=azSw?-90:-180;f.max=azSw?90:180;f.step=0.5;
  f.value=azSw?sFixEl:sFixAz;
  $('fixv').textContent=(azSw?sFixEl:sFixAz)+'\u00B0';
}
document.querySelectorAll('#ds button').forEach(b=>b.onclick=()=>{
  seg('#ds',b);sDs=b.dataset.d;refreshSlice();});  // angles persist: no reset
document.querySelectorAll('#stype button').forEach(b=>b.onclick=()=>{
  seg('#stype',b);sType=b.dataset.s;syncFixSlider();refreshSlice();});
document.querySelectorAll('#sfscale button').forEach(b=>b.onclick=()=>{
  seg('#sfscale',b);sFscale=b.dataset.f;refreshSlice();});
$('norm').addEventListener('input',()=>{sNorm=$('norm').checked;refreshSlice();});
$('lmatch').addEventListener('input',()=>{refreshSlice();});
['cc_m','cc_h','cc_s'].forEach(id=>$(id).addEventListener('input',()=>{refreshCues();}));
$('fix').addEventListener('input',()=>{const v=+$('fix').value;
  if(sType==='az')sFixEl=v;else sFixAz=v;
  $('fixv').textContent=v+'\u00B0';scheduleS();});
function diffuseFor(){
  // Global diffuse-field response (dB, log-f grid) for the current dataset:
  // the area-weighted power mean over ALL measurement directions.
  if(sDs==='soni')return DF.soni;
  if(sDs==='human')return DF.human;
  if(sDs==='diff')return{
    L:DF.meta.L.map((v,j)=>v-DF.soni.L[j]),
    R:DF.meta.R.map((v,j)=>v-DF.soni.R[j])};
  return DF.meta;}
function removeDiffuse(cols,df){  // subtract the diffuse-field response
  return cols.map(c=>c&&c.map((v,j)=>v-df[j]));}
function buildSlice(){
  const azSw=sType==='az', n=181;
  sAngles=new Array(n);sColsL=new Array(n);sColsR=new Array(n);sGaps=0;
  const dsLabel=sDs==='meta'?LBL.meta:sDs==='human'?LBL.human:sDs==='diff'?'\u0394 EarSim\u2212BTE (dB)':LBL.soni;
  const lm=$('lmatch').checked&&(sDs==='soni'||sDs==='human'||sDs==='diff')&&!sNorm;
  if(sDs==='diff'){
    // Per-direction EarSim-minus-BTE magnitude spectra (dB), per ear.
    // Defined only where both datasets cover the direction: below-rim
    // columns are omitted, never extrapolated. The level-match box removes
    // the global EarSim-minus-BTE level offset before differencing.
    const off=lm?LVL:0;
    for(let i=0;i<n;i++){
      const au=azSw?-180+i*2:sFixAz, e=azSw?sFixEl:-90+i*1;
      sAngles[i]=azSw?au:e;  // UI convention: + = right
      const qm=query(IM,toDataAz(au),e), qs=query(IS,toDataAz(au),e);
      if(!qm||!qs){sColsL[i]=null;sColsR[i]=null;sGaps++;}
      else{
        const cL=new Array(NF),cR=new Array(NF);
        for(let j=0;j<NF;j++){cL[j]=qm.L[j]-(qs.L[j]+off);cR[j]=qm.R[j]-(qs.R[j]+off);}
        sColsL[i]=cL;sColsR[i]=cR;
      }
    }
  }else{
    const I=sDs==='meta'?IM:sDs==='human'?IH:IS;
    const offDb=sDs==='human'?(DATA_ALL.ranges.lvl_h||0):LVL;
    for(let i=0;i<n;i++){
      const au=azSw?-180+i*2:sFixAz, e=azSw?sFixEl:-90+i*1;
      sAngles[i]=azSw?au:e;  // UI convention: + = right
      const q=query(I,toDataAz(au),e);
      if(!q){sColsL[i]=null;sColsR[i]=null;sGaps++;}
      else{sColsL[i]=q.L;sColsR[i]=q.R;}
    }
    if(lm)for(let i=0;i<n;i++){
      if(sColsL[i])sColsL[i]=sColsL[i].map(v=>v+offDb);
      if(sColsR[i])sColsR[i]=sColsR[i].map(v=>v+offDb);}
  }
  const lmchk=$('lmchk');
  const lmOk=(sDs==='soni'||sDs==='human'||sDs==='diff')&&!sNorm;
  $('lmatch').disabled=!lmOk;
  lmchk.classList.toggle('disabled',!lmOk);
  if(sNorm){const df=diffuseFor();sColsL=removeDiffuse(sColsL,df.L);sColsR=removeDiffuse(sColsR,df.R);sVmin=-YNORM;sVmax=YNORM;}
  else if(sDs==='diff'){sVmin=-YDIF;sVmax=YDIF;}
  else{sVmin=YMAG[0];sVmax=YMAG[1];}
  $('sliceinfo').innerHTML=dsLabel+' \u00B7 '+
    (azSw?'azimuth sweep @ el ':'elevation sweep @ az ')+(azSw?sFixEl:sFixAz)+'\u00B0 \u00B7 '+
    n+(sDs==='diff'?' difference':' interpolated')+' columns \u00D7 both ears<br>'+
    'fixed '+sVmin.toFixed(0)+' \u2026 '+sVmax.toFixed(0)+' dB (same every slice)'+
    (sNorm?' (diffuse-field removed)':'')+(sDs==='diff'&&!sNorm?'; red = EarSim louder, blue = BTE louder':'')+(lm?' \u00B7 level-matched':'')+
    (sGaps?'<br><span class="warn">'+sGaps+' of '+n+' columns outside measurement grid (omitted)</span>':'');
}
function drawHeat(cv,cols,freqs,earLabel){
  const dpr=window.devicePixelRatio||1,W=cv.clientWidth,H=cv.clientHeight;
  if(!W||!H||!cols.length)return null;
  cv.width=W*dpr;cv.height=H*dpr;
  const c=cv.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);
  c.fillStyle='#0b1016';c.fillRect(0,0,W,H);
  const n=cols.length,L=56,R=W-64,T=30,B=46,pw=R-L,ph=H-T-B;
  const nR=freqs.length, lin=freqs!==F;
  const LUT=(sNorm||sDs==='diff')?LUT_D:LUT_I;
  const img=c.createImageData(n,nR);
  for(let cc=0;cc<n;cc++){const col=cols[cc];
    for(let r=0;r<nR;r++){
      const o=(r*n+cc)*4;
      if(!col){img.data[o]=GAPRGB[0];img.data[o+1]=GAPRGB[1];img.data[o+2]=GAPRGB[2];img.data[o+3]=255;continue;}
      const v=col[nR-1-r],t=(v-sVmin)/(sVmax-sVmin);
      const rgb=LUT[Math.min(255,Math.max(0,Math.round(t*255)))];
      img.data[o]=rgb[0];img.data[o+1]=rgb[1];img.data[o+2]=rgb[2];img.data[o+3]=255;}}
  const off=document.createElement('canvas');off.width=n;off.height=nR;
  off.getContext('2d').putImageData(img,0,0);
  c.imageSmoothingEnabled=true;c.imageSmoothingQuality='high';
  c.drawImage(off,L,T,pw,ph);
  c.strokeStyle='#3a4d63';c.lineWidth=1;c.strokeRect(L,T,pw,ph);
  c.font='11px '+FONT;c.fillStyle=DIM;c.textAlign='left';
  [100,200,500,1000,2000,5000,10000,20000].forEach(f=>{
    const y=lin?T+(1-(f-100)/(20000-100))*ph:T+(1-(Math.log(f)-Math.log(100))/(Math.log(20000)-Math.log(100)))*ph;
    c.fillText(f>=1000?(f/1000)+'k':f,8,y+4);
    c.strokeStyle='rgba(255,255,255,0.10)';c.beginPath();c.moveTo(L,y);c.lineTo(R,y);c.stroke();});
  const stepA=sType==='az'?30:15,startA=Math.ceil(sAngles[0]/stepA)*stepA;
  c.fillStyle=INK;c.textAlign='center';
  for(let a=startA;a<=sAngles[n-1];a+=stepA){
    let bi=0,bd=1e9;for(let i=0;i<n;i++){const d=Math.abs(sAngles[i]-a);if(d<bd){bd=d;bi=i;}}
    const x=L+(bi+0.5)/n*pw;c.fillText(a+'\u00B0',x,H-B+18);}
  c.textAlign='left';c.fillStyle=DIM;
  c.fillText(sType==='az'?'azimuth (\u00B0)':'elevation (\u00B0)',L,H-8);
  c.fillStyle=INK;c.font='600 13px '+FONT;
  c.fillText(earLabel,L,17);
  const cbW=14,cbX=R+10;
  for(let y=0;y<ph;y++){const t=1-y/ph,rgb=LUT[Math.round(t*255)];
    c.fillStyle='rgb('+rgb[0]+','+rgb[1]+','+rgb[2]+')';c.fillRect(cbX,T+y,cbW,1.5);}
  c.fillStyle=DIM;c.font='11px '+FONT;c.textAlign='left';
  c.fillText(sVmax.toFixed(0)+' dB',cbX-4,T-4);
  c.fillText(((sVmin+sVmax)/2).toFixed(0),cbX+2,T+ph/2);
  c.fillText(sVmin.toFixed(0)+' dB',cbX-4,T+ph+14);
  return {L:L,T:T,pw:pw,ph:ph,n:n,freqs:freqs,cols:cols};
}
function drawSlice(){
  const lin=sFscale==='lin',fr=lin?LIN_F:F;
  sGeoL=drawHeat(sliceL,lin?toLinear(sColsL):sColsL,fr,'Left ear');
  sGeoR=drawHeat(sliceR,lin?toLinear(sColsR):sColsR,fr,'Right ear');}
function refreshSlice(){buildSlice();drawSlice();refreshCues();}
function hoverSlice(e,cv,geo,cols,earLabel){
  if(!geo||!geo.cols.length)return;
  const r=cv.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;
  const L=geo.L,T=geo.T,pw=geo.pw,ph=geo.ph,n=geo.n;
  if(x<L||x>L+pw||y<T||y>T+ph)return;
  const fr=geo.freqs,ccols=geo.cols,nR=fr.length;
  const cc=Math.min(n-1,Math.floor((x-L)/pw*n)),row=Math.min(nR-1,Math.floor((y-T)/ph*nR));
  const ang=sAngles[cc],alab=sType==='az'?'az':'el';
  const col=ccols[cc];
  if(!col){
    $('sliceread').innerHTML=earLabel+' \u00B7 '+alab+' '+ang+'\u00B0: <span class="warn">outside measurement grid</span>';
    return;
  }
  const f=fr[nR-1-row],v=col[nR-1-row];
  const fs=f>=1000?(f/1000).toFixed(1)+' kHz':Math.round(f)+' Hz';
  const tag=sNorm?' (diffuse-field removed)':'';
  $('sliceread').innerHTML=earLabel+' \u00B7 '+alab+' '+ang+'\u00B0 \u00B7 '+fs+': <b>'+v.toFixed(1)+' dB</b>'+tag;
  $('hover').innerHTML=earLabel+' \u00B7 '+alab+' '+ang+'\u00B0 \u00B7 '+fs+' \u00B7 <b>'+v.toFixed(1)+' dB</b>';
}
sliceL.addEventListener('mousemove',e=>hoverSlice(e,sliceL,sGeoL,sColsL,'Left'));
sliceR.addEventListener('mousemove',e=>hoverSlice(e,sliceR,sGeoR,sColsR,'Right'));

/* ---------- binaural cues along the slice sweep (same tab) ---------- */
const ITD_R=850, ILD_R=25;  // fixed symmetric ranges (vertex maxima: 833 us, 21.5 dB)
const ITD_COMP=0;   // +20.83 us: removes EarSim's constant 1-sample inter-channel ITD bias
let compITD=true, CUE=null, CUEG=null;
function setCompITD(v){compITD=v;
  $('itdcomp1').checked=v;$('itdcomp2').checked=v;
  refresh1();refreshCues();if(playing)pushIR();}
$('itdcomp1').onchange=e=>setCompITD(e.target.checked);
$('itdcomp2').onchange=e=>setCompITD(e.target.checked);
function drawLine(cv,title,angs,series,yr,ystep,xstep,names,cols){
  const dpr=window.devicePixelRatio||1,W=560,H=300;
  cv.width=W*dpr;cv.height=H*dpr;cv.style.width='100%';
  const g=cv.getContext('2d');g.scale(dpr,dpr);
  g.fillStyle='#0d141c';g.fillRect(0,0,W,H);
  const L=54,R=12,T=26,B=30,pw=W-L-R,ph=H-T-B;
  const x0=angs[0],x1=angs[angs.length-1];
  const X=v=>L+(v-x0)/(x1-x0)*pw, Y=v=>T+ph/2-(v/yr)*(ph/2);
  g.font='10px '+FONT;
  g.textAlign='right';g.textBaseline='middle';
  for(let v=-yr;v<=yr+1e-9;v+=ystep){
    const y=Y(v);
    g.strokeStyle=GRID;g.beginPath();g.moveTo(L,y);g.lineTo(W-R,y);g.stroke();
    g.fillStyle=DIM;g.fillText(String(Math.round(v)),L-6,y);
  }
  g.textAlign='center';g.textBaseline='top';
  for(let v=Math.ceil(x0/xstep)*xstep;v<=x1+1e-9;v+=xstep){
    g.fillStyle=DIM;g.fillText(v+'\u00b0',X(v),T+ph+6);
  }
  g.strokeStyle=DIM;g.beginPath();g.moveTo(L,Y(0));g.lineTo(W-R,Y(0));g.stroke();
  series.forEach((s,si)=>{
    g.strokeStyle=cols[si];g.lineWidth=1.6;g.beginPath();
    let pen=false;
    for(let i=0;i<s.length;i++){
      if(s[i]==null){pen=false;continue;}
      const x=X(angs[i]),y=Y(s[i]);
      if(!pen){g.moveTo(x,y);pen=true;}else g.lineTo(x,y);
    }
    g.stroke();
  });
  g.fillStyle=INK;g.textAlign='left';g.textBaseline='top';
  g.fillText(title,L,T-18);
  g.textAlign='right';
  let lx=W-R-8;
  for(let si=names.length-1;si>=0;si--){
    const nm=names[si],w=g.measureText(nm).width;
    lx-=w+30;
    g.strokeStyle=cols[si];g.lineWidth=2;g.beginPath();
    g.moveTo(lx,T-12);g.lineTo(lx+16,T-12);g.stroke();
    g.fillStyle=INK;g.fillText(nm,lx+20,T-18);
    lx-=6;
  }
  return {L:L,R:R,W:W,pw:pw,n:angs.length};
}
function drawCuePlots(mitd,mild,sitd,sild,hitd,hild){
  const showM=$('cc_m').checked, showS=$('cc_s').checked;
  const showH=HAS_HUMAN&&$('cc_h').checked&&hitd&&hild;
  const names=[],cols=[],itdSeries=[],ildSeries=[];
  if(showM){names.push(LBL.meta);cols.push(COLM);itdSeries.push(mitd);ildSeries.push(mild);}
  if(showH){names.push(LBL.human);cols.push(COLH);itdSeries.push(hitd);ildSeries.push(hild);}
  if(showS){names.push(LBL.soni);cols.push(COLS);itdSeries.push(sitd);ildSeries.push(sild);}
  CUEG=drawLine($('itdplot'),'ITD (\u00b5s)',sAngles,itdSeries,ITD_R,250,sType==='az'?60:30,names,cols);
  drawLine($('ildplot'),'ILD (dB)',sAngles,ildSeries,ILD_R,5,sType==='az'?60:30,names,cols);
}
let shCueScheduled=false;
function refreshCues(){
  // Same sweep as heatmaps; all three datasets use barycentric interpolation.
  const n=sAngles.length;
  const mitd=new Array(n),mild=new Array(n),sitd=new Array(n),sild=new Array(n);
  const hitd=new Array(n),hild=new Array(n);
  for(let i=0;i<n;i++){
    const az=toDataAz(sType==='az'?sAngles[i]:sFixAz), el=sType==='az'?sFixEl:sAngles[i];
    const qm=query(IM,az,el),qs=query(IS,az,el),qh=HAS_HUMAN?query(IH,az,el):null;
    mitd[i]=qm?(qm.itd+(compITD?ITD_COMP:0)):null; mild[i]=qm?qm.ild:null;
    sitd[i]=qs?qs.itd:null; sild[i]=qs?qs.ild:null;
    hitd[i]=qh?qh.itd:null; hild[i]=qh?qh.ild:null;
  }
  drawCuePlots(mitd,mild,sitd,sild,hitd,hild);
  CUE={angs:sAngles,mitd:mitd,mild:mild,sitd:sitd,sild:sild,hitd:hitd,hild:hild};
  $('cueinfo').innerHTML=(sType==='az'?'Azimuth sweep at elevation ':'Elevation sweep at azimuth ')+
    (sType==='az'?sFixEl:sFixAz)+'&deg;';
}
function hoverCues(e){
  if(!CUE||!CUEG)return;
  const r=e.target.getBoundingClientRect(),px=e.clientX-r.left;
  const x=px/r.width*CUEG.W;
  let i=Math.round((x-CUEG.L)/CUEG.pw*(CUEG.n-1));
  i=Math.max(0,Math.min(CUEG.n-1,i));
  const alab=sType==='az'?'az':'el',ang=CUE.angs[i];
  const f1=v=>v==null?'<span class="warn">n/a</span>':'<b>'+v.toFixed(0)+'</b>';
  const f2=v=>v==null?'<span class="warn">n/a</span>':'<b>'+v.toFixed(1)+'</b>';
  const parts=[];
  if($('cc_m').checked)parts.push(LBL.meta+' '+f1(CUE.mitd[i])+' &micro;s');
  if(HAS_HUMAN&&$('cc_h').checked&&CUE.hitd)parts.push(LBL.human+' '+f1(CUE.hitd[i])+' &micro;s');
  if($('cc_s').checked)parts.push(LBL.soni+' '+f1(CUE.sitd[i])+' &micro;s');
  const ilds=[];
  if($('cc_m').checked)ilds.push(LBL.meta+' '+f2(CUE.mild[i])+' dB');
  if(HAS_HUMAN&&$('cc_h').checked&&CUE.hild)ilds.push(LBL.human+' '+f2(CUE.hild[i])+' dB');
  if($('cc_s').checked)ilds.push(LBL.soni+' '+f2(CUE.sild[i])+' dB');
  $('cueread').innerHTML=alab+' '+ang+'&deg; &middot; ITD: '+parts.join(', ')+' &middot; ILD: '+ilds.join(', ');
}
$('itdplot').addEventListener('mousemove',hoverCues);
$('ildplot').addEventListener('mousemove',hoverCues);

window.addEventListener('resize',()=>{if(!$('panel1').classList.contains('hidden'))refresh1();else refreshSlice();});
syncFixSlider();
$('lmatchv').textContent='(+'+LVL.toFixed(1)+' dB)';
if(HAS_HUMAN){
  $('chk_h').style.display='';$('th_h').style.display='';
  $('itd_h').style.display='';$('ild_h').style.display='';
  $('ds_h').style.display='';
  const pb=$('psrc_h');if(pb){pb.style.display='';pb.dataset.d='human';pb.textContent=LBL.human;}
  $('cchk_h').style.display='';
  $('c_h').checked=true;$('cc_h').checked=true;
  $('shfooter').innerHTML='P0006: SONICOM measured human HRTF (Windowed 48 kHz, with ITD) — reference/benchmark against KEMAR EarSim and BTE.';
}else{
  $('chk_h').style.display='none';$('th_h').style.display='none';
  $('itd_h').style.display='none';$('ild_h').style.display='none';
  $('ds_h').style.display='none';
  const pb=$('psrc_h');if(pb)pb.style.display='none';
  $('cchk_h').style.display='none';
  $('shfooter').style.display='none';
}
document.querySelectorAll('.shhint').forEach(el=>{el.style.display='none';});
refresh1();

})().catch(err => {
  console.error(err);
  const d = document.createElement('div');
  d.style.cssText = 'margin:1em;padding:1em;border:1px solid #a33;background:#fee;color:#900;font-family:sans-serif';
  d.textContent = 'Could not load page data (' + err.message + '). ' +
    'Serve this folder over HTTP — e.g. run `python3 -m http.server` here — or use GitHub Pages.';
  document.body.prepend(d);
});
