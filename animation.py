# SCAL Studio - made by Merouane Hamdani for educational purposes only.
# Released under the MIT License (see LICENSE).
"""animation.py - self-contained HTML/JS (canvas) animations embedded with st.components.v1.html.

Tab 1  Waterflood front : Buckley-Leverett (Welge) displacement in a core plug, live saturation
                          profile, fractional-flow curve with tangent, recovery vs PVI.
Tab 2  Hysteresis loop  : saturation cycle walking along drainage and scanning curves.
All computation for the front is done in JS so the viscosity ratio slider is instantaneous.
"""
from __future__ import annotations

import json

HTML = r"""
<!DOCTYPE html><html><head><meta charset="utf-8">
<style>
 :root{--navy:#051C2C;--blue:#2251FF;--cyan:#00A9F4;--grey:#6B7C93;--line:#D5DEEA;--bg:#FFFFFF;--soft:#F3F6FB;--oil:#B9822B}
 *{box-sizing:border-box} body{margin:0;font-family:Arial,Helvetica,sans-serif;color:var(--navy);background:var(--bg)}
 .tabs{display:flex;gap:0;border-bottom:2px solid var(--navy);margin-bottom:10px}
 .tab{padding:8px 16px;cursor:pointer;font-size:13px;color:var(--grey);border:0;background:none;font-weight:600}
 .tab.on{color:var(--bg);background:var(--navy)}
 .panel{display:none}.panel.on{display:block}
 .row{display:flex;gap:14px;flex-wrap:wrap;align-items:center;margin:6px 0 10px}
 .row label{font-size:12px;color:var(--grey)} .row b{color:var(--navy);font-size:12px}
 button.p{background:var(--blue);color:#fff;border:0;padding:7px 16px;font-weight:700;cursor:pointer;font-size:12px;letter-spacing:.3px}
 button.p:hover{background:#1a3fd6}
 input[type=range]{accent-color:var(--blue);vertical-align:middle;width:150px}
 .grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
 canvas{width:100%;background:var(--bg);border:1px solid var(--line)}
 .cap{font-size:11px;color:var(--grey);margin:3px 0 0}
 .kpi{display:flex;gap:10px;margin:8px 0}
 .kpi div{flex:1;background:var(--soft);border-left:3px solid var(--blue);padding:6px 10px;font-size:11px;color:var(--grey)}
 .kpi div span{display:block;font-family:Georgia,serif;font-size:19px;color:var(--navy)}
</style></head><body>
<div class="tabs"><button class="tab on" data-t="a">Waterflood front</button><button class="tab" data-t="b">Hysteresis loop</button></div>

<div class="panel on" id="a">
 <div class="row">
  <button class="p" id="play">&#9654; Play</button>
  <label>Pore volumes injected <b id="pvv">0.00</b></label><input type="range" id="pv" min="0" max="3" step="0.005" value="0">
  <label>Speed</label><input type="range" id="spd" min="0.2" max="4" step="0.1" value="1" style="width:90px">
  <label>&mu;<sub>o</sub>/&mu;<sub>w</sub> <b id="mrv">10.0</b></label><input type="range" id="mr" min="0.5" max="100" step="0.5" value="10" style="width:120px">
 </div>
 <canvas id="core" height="92"></canvas>
 <p class="cap">Core plug (injection on the left). Blue = water, amber = oil. Saturation front from Buckley&ndash;Leverett / Welge theory, no capillary effects.</p>
 <div class="kpi"><div>Breakthrough<span id="kbt">-</span></div><div>RF at breakthrough<span id="krfb">-</span></div><div>RF now<span id="krf">-</span></div><div>Water cut now<span id="kwc">-</span></div></div>
 <div class="grid"><div><canvas id="prof" height="230"></canvas><p class="cap">Water saturation profile S<sub>w</sub>(x<sub>D</sub>)</p></div>
 <div><canvas id="ffl" height="230"></canvas><p class="cap">Fractional flow with Welge tangent</p></div></div>
 <div style="margin-top:12px"><canvas id="rec" height="200"></canvas><p class="cap">Oil recovery factor vs pore volumes injected</p></div>
</div>

<div class="panel" id="b">
 <div class="row">
  <button class="p" id="hplay">&#9654; Play cycle</button>
  <label>Reversal saturation <b id="hsv">-</b></label><input type="range" id="hs" min="0" max="1" step="1" value="0">
  <label id="hlabel" style="color:var(--navy);font-weight:700"></label>
 </div>
 <div class="grid"><div><canvas id="hk" height="300"></canvas><p class="cap">Non-wetting phase relative permeability: drainage (blue), boundary imbibition (grey), scanning (cyan)</p></div>
 <div><canvas id="hsat" height="300"></canvas><p class="cap">Saturation history driving the cycle</p></div></div>
 <div class="kpi"><div>Current S<sub>nw</sub><span id="hcs">-</span></div><div>Current kr<span id="hck">-</span></div><div>Trapped S<sub>nw</sub> (this reversal)<span id="hct">-</span></div></div>
</div>

<p class="cap" style="margin-top:14px">SCAL Studio &middot; made by Merouane Hamdani &middot; educational use only &middot; MIT License</p>
<script>
const D = __PAYLOAD__;
const C = {navy:'#051C2C',blue:'#2251FF',cyan:'#00A9F4',grey:'#6B7C93',line:'#D5DEEA',oil:'#B9822B'};
const $ = id => document.getElementById(id);
function fit(c){const r=window.devicePixelRatio||1,h=+(c.dataset.logicalHeight||c.getAttribute('height'));c.dataset.logicalHeight=h;c.style.height=h+'px';const w=c.clientWidth;c.width=w*r;c.height=h*r;const x=c.getContext('2d');x.setTransform(r,0,0,r,0,0);x.W=w;x.H=h;return x;}
function interp(xs,ys,x){const n=xs.length;if(x<=xs[0])return ys[0];if(x>=xs[n-1])return ys[n-1];let lo=0,hi=n-1;while(hi-lo>1){const m=(lo+hi)>>1;if(xs[m]<=x)lo=m;else hi=m;}return ys[lo]+(ys[hi]-ys[lo])*(x-xs[lo])/(xs[hi]-xs[lo]);}

/* ---------- generic axes plot ---------- */
function axes(x,o){const m={l:44,r:12,t:12,b:34},fmt=(v,span)=>v.toFixed(Math.max(2,Math.min(4,Math.ceil(-Math.log10(span/5))+1)));x.clearRect(0,0,x.W,x.H);x.m=m;x.sx=v=>m.l+(v-o.x0)/(o.x1-o.x0)*(x.W-m.l-m.r);x.sy=v=>x.H-m.b-(v-o.y0)/(o.y1-o.y0)*(x.H-m.t-m.b);
 x.font='11px Arial';x.fillStyle=C.grey;x.strokeStyle=C.line;x.lineWidth=1;
 for(let i=0;i<=5;i++){const vx=o.x0+(o.x1-o.x0)*i/5,vy=o.y0+(o.y1-o.y0)*i/5;
  x.beginPath();x.moveTo(x.sx(vx),x.sy(o.y0));x.lineTo(x.sx(vx),x.sy(o.y1));x.stroke();
  x.beginPath();x.moveTo(x.sx(o.x0),x.sy(vy));x.lineTo(x.sx(o.x1),x.sy(vy));x.stroke();
  x.textAlign='center';x.fillText(fmt(vx,o.x1-o.x0),x.sx(vx),x.H-m.b+14);x.textAlign='right';x.fillText(fmt(vy,o.y1-o.y0),m.l-5,x.sy(vy)+4);}
 x.strokeStyle=C.navy;x.lineWidth=1.4;x.beginPath();x.moveTo(m.l,m.t);x.lineTo(m.l,x.H-m.b);x.lineTo(x.W-m.r,x.H-m.b);x.stroke();
 x.fillStyle=C.navy;x.textAlign='center';x.fillText(o.xl||'',x.W/2+10,x.H-4);x.save();x.translate(11,x.H/2);x.rotate(-Math.PI/2);x.fillText(o.yl||'',0,0);x.restore();}
function line(x,xs,ys,col,w,dash){x.strokeStyle=col;x.lineWidth=w||2;x.setLineDash(dash||[]);x.beginPath();xs.forEach((v,i)=>{const px=x.sx(v),py=x.sy(ys[i]);i?x.lineTo(px,py):x.moveTo(px,py);});x.stroke();x.setLineDash([]);}
function dot(x,vx,vy,col,r){x.fillStyle=col;x.beginPath();x.arc(x.sx(vx),x.sy(vy),r||5,0,6.2832);x.fill();x.strokeStyle='#fff';x.lineWidth=1.5;x.stroke();}

/* ---------- Waterflood (Welge) ---------- */
const wo=D.wo, swi=wo.sw[0];
let smax=wo.sw[0]; wo.kro.forEach((k,i)=>{if(k>1e-9)smax=Math.max(smax,wo.sw[i]);});
const N=1500, SW=[...Array(N)].map((_,i)=>swi+(smax-swi)*i/(N-1));
const KRW=SW.map(s=>interp(wo.sw,wo.krw,s)), KRO=SW.map(s=>interp(wo.sw,wo.kro,s));
let FW,DFW,sf,fwf,slope,idx,tbt,rev;
function solve(){const mr=+$('mr').value; const muw=1, muo=mr;
 FW=SW.map((s,i)=>{const a=KRW[i]/muw,b=KRO[i]/muo;return a+b>0?a/(a+b):0;});
 DFW=FW.map((f,i)=>i==0?(FW[1]-FW[0])/(SW[1]-SW[0]):i==N-1?(FW[N-1]-FW[N-2])/(SW[N-1]-SW[N-2]):(FW[i+1]-FW[i-1])/(SW[i+1]-SW[i-1]));
 let best=-1;idx=1;for(let i=1;i<N;i++){const sl=(FW[i]-FW[0])/(SW[i]-swi);if(sl>best){best=sl;idx=i;}}
 slope=best;sf=SW[idx];fwf=FW[idx];tbt=1/slope;
 rev={d:DFW.slice(idx).reverse(),s:SW.slice(idx).reverse(),f:FW.slice(idx).reverse()};}
function state(t){ // returns {sw2 (outlet), fw2, np (PV of oil produced)}
 if(t<=tbt) return {sw2:swi,fw2:FW[0],np:t};
 const sw2=interp(rev.d,rev.s,1/t), fw2=interp(SW,FW,sw2);
 return {sw2,fw2,np:Math.min(sw2+(1-fw2)*t-swi,smax-swi)};}
function profile(xd,t){ if(t<=0) return swi; const v=xd/t; if(v>slope) return swi; return interp(rev.d,rev.s,v);} 
function drawCore(t){const x=fit($('core')),W=x.W,H=x.H,m=10,w=W-2*m,h=H-26;
 const cols=200; for(let i=0;i<cols;i++){const xd=(i+.5)/cols,sw=profile(xd,t),f=(sw-swi)/(smax-swi);
  const r=Math.round(185*(1-f)+34*f),g=Math.round(130*(1-f)+81*f),b=Math.round(43*(1-f)+255*f);x.fillStyle=`rgb(${r},${g},${b})`;x.fillRect(m+i*w/cols,12,w/cols+1,h);}
 x.strokeStyle=C.navy;x.lineWidth=2;x.strokeRect(m,12,w,h);
 x.fillStyle=C.navy;x.font='11px Arial';x.textAlign='left';x.fillText('Injector',m,H-4);x.textAlign='right';x.fillText('Producer',W-m,H-4);
 if(t>0&&t<=tbt){const px=m+Math.min(slope*t,1)*w;x.strokeStyle='#fff';x.setLineDash([4,3]);x.beginPath();x.moveTo(px,12);x.lineTo(px,12+h);x.stroke();x.setLineDash([]);}}
function drawProf(t){const x=fit($('prof'));axes(x,{x0:0,x1:1,y0:swi,y1:smax,xl:'Dimensionless distance x_D',yl:'Water saturation Sw'});
 const xs=[],ys=[];for(let i=0;i<=300;i++){const xd=i/300;xs.push(xd);ys.push(profile(xd,t));}line(x,xs,ys,C.blue,2.6);
 x.fillStyle=C.grey;x.font='11px Arial';x.textAlign='left';x.fillText('Sw,front = '+sf.toFixed(3),x.sx(.55),x.sy(smax)+16);}
function drawFF(t){const x=fit($('ffl'));axes(x,{x0:swi,x1:smax,y0:0,y1:1,xl:'Water saturation Sw',yl:'Fractional flow fw'});
 line(x,SW,FW,C.navy,2.4);line(x,[swi,swi+1/slope],[FW[0],FW[0]+1],C.cyan,1.6,[5,4]);dot(x,sf,fwf,C.cyan,5);
 const s=state(t);dot(x,s.sw2,s.fw2,C.blue,6);}
function drawRec(t){const x=fit($('rec')),pm=+$('pv').max;const rfu=(smax-swi)/(1-swi);
 axes(x,{x0:0,x1:pm,y0:0,y1:Math.min(1,rfu*1.15),xl:'Pore volumes injected',yl:'Recovery factor'});
 const xs=[],ys=[];for(let i=0;i<=300;i++){const tt=pm*i/300;xs.push(tt);ys.push(state(tt).np/(1-swi));}line(x,xs,ys,C.navy,2.4);
 x.strokeStyle=C.grey;x.setLineDash([4,4]);x.beginPath();x.moveTo(x.sx(0),x.sy(rfu));x.lineTo(x.sx(pm),x.sy(rfu));x.stroke();x.setLineDash([]);
 x.fillStyle=C.grey;x.textAlign='left';x.fillText('Displacement efficiency limit '+(rfu*100).toFixed(0)+'%',x.sx(pm*.35),x.sy(rfu)-5);
 dot(x,t,state(t).np/(1-swi),C.blue,6);x.strokeStyle=C.cyan;x.setLineDash([3,3]);x.beginPath();x.moveTo(x.sx(tbt),x.sy(0));x.lineTo(x.sx(tbt),x.sy(rfu));x.stroke();x.setLineDash([]);}
let T=0,playing=false,last=0;
function renderA(){const t=T;$('pvv').textContent=t.toFixed(2);$('pv').value=t;$('mrv').textContent=(+$('mr').value).toFixed(1);
 drawCore(t);drawProf(t);drawFF(t);drawRec(t);const s=state(t);
 $('kbt').textContent=tbt.toFixed(2)+' PV';$('krfb').textContent=(tbt/(1-swi)*100).toFixed(0)+'%';$('krf').textContent=(s.np/(1-swi)*100).toFixed(1)+'%';$('kwc').textContent=(s.fw2*100).toFixed(0)+'%';}
function loop(ts){if(!playing)return;const dt=(ts-last)/1000;last=ts;T+=dt*0.35*(+$('spd').value);if(T>=+$('pv').max){T=+$('pv').max;playing=false;$('play').innerHTML='&#9654; Play';}renderA();requestAnimationFrame(loop);}
$('play').onclick=()=>{if(T>=+$('pv').max)T=0;playing=!playing;$('play').innerHTML=playing?'&#10074;&#10074; Pause':'&#9654; Play';if(playing){last=performance.now();requestAnimationFrame(loop);}};
$('pv').oninput=e=>{T=+e.target.value;renderA();};$('mr').oninput=()=>{solve();renderA();};$('spd').oninput=()=>{};
$('mr').value=D.mu_ratio;

/* ---------- Hysteresis ---------- */
const H=D.hyst; let hi=0,hp=0,hplaying=false,hlast=0;
if(H){$('hs').max=H.scans.length-1;$('hlabel').textContent=H.label;}
function path(i){const sc=H.scans[i],dr=H.drain,pts=[];
 for(let k=0;k<dr.s.length;k++){if(dr.s[k]>=H.s_crit&&dr.s[k]<=sc.s_hy)pts.push([dr.s[k],dr.kr[k]]);}
 pts.push([sc.s_hy,interp(dr.s,dr.kr,sc.s_hy)]);
 for(let k=sc.s.length-1;k>=0;k--)pts.push([sc.s[k],sc.kr[k]]);return pts;}
function drawH(){if(!H)return;const x=fit($('hk')),y=fit($('hsat'));$('hsv').textContent=H.scans[hi].s_hy.toFixed(3);
 const P=path(hi);if(P.length<2)return;
 const n=P.length-1,f=Math.max(0,Math.min(1,hp))*n,k=Math.min(Math.floor(f),n-1),a=P[k],b=P[k+1],u=f-k;
 const cs=a[0]+(b[0]-a[0])*u,ck=a[1]+(b[1]-a[1])*u;
 const sx=P.map(p=>p[0]),ky=P.map(p=>p[1]),x0=Math.min(...sx),x1=Math.max(...sx),y1=Math.max(...ky,1e-6)*1.15;
 const visible=c=>c.s.map((s,i)=>[s,c.kr[i]]).filter(p=>p[0]>=x0&&p[0]<=x1);
 axes(x,{x0,x1,y0:0,y1,xl:'Non-wetting saturation',yl:'Relative permeability'});
 const dr=visible(H.drain),bd=visible(H.bound),sc=visible(H.scans[hi]);
 line(x,dr.map(p=>p[0]),dr.map(p=>p[1]),C.blue,2.6);line(x,bd.map(p=>p[0]),bd.map(p=>p[1]),C.grey,2,[6,4]);
 line(x,sc.map(p=>p[0]),sc.map(p=>p[1]),C.cyan,2.8);
 dot(x,cs,ck,C.navy,6.5);
 const pad=Math.max((x1-x0)*0.1,0.01);
 axes(y,{x0:0,x1:1,y0:Math.max(0,x0-pad),y1:Math.min(H.s_max*1.05,x1+pad),xl:'Cycle progress',yl:'Non-wetting saturation'});
 const xs=P.map((p,i)=>i/n),ys=P.map(p=>p[0]);line(y,xs,ys,C.blue,2.4);dot(y,hp,cs,C.navy,6.5);
 $('hcs').textContent=cs.toFixed(3);$('hck').textContent=ck.toFixed(4);$('hct').textContent=H.scans[hi].s_trap.toFixed(3);}
function hloop(ts){if(!hplaying)return;const dt=(ts-hlast)/1000;hlast=ts;hp+=dt*0.18;if(hp>=1){hp=1;hplaying=false;$('hplay').innerHTML='&#9654; Play cycle';}drawH();requestAnimationFrame(hloop);}
$('hplay').onclick=()=>{if(hp>=1)hp=0;hplaying=!hplaying;$('hplay').innerHTML=hplaying?'&#10074;&#10074; Pause':'&#9654; Play cycle';if(hplaying){hlast=performance.now();requestAnimationFrame(hloop);}};
$('hs').oninput=e=>{hi=+e.target.value;hp=0;drawH();};

/* ---------- tabs ---------- */
document.querySelectorAll('.tab').forEach(b=>b.onclick=()=>{document.querySelectorAll('.tab').forEach(t=>t.classList.remove('on'));b.classList.add('on');
 document.querySelectorAll('.panel').forEach(p=>p.classList.remove('on'));$(b.dataset.t).classList.add('on');if(b.dataset.t=='a')renderA();else drawH();});
window.addEventListener('resize',()=>{renderA();if($('b').classList.contains('on'))drawH();});
solve();renderA();
</script></body></html>
"""


def build_animation_html(wo: dict, mu_ratio: float, hyst: dict | None = None) -> str:
    """wo: {sw,krw,kro} lists. hyst: {label,s_crit,s_max,drain{s,kr},bound{s,kr},scans[{s_hy,s_trap,s,kr}]}"""
    payload = {"wo": wo, "mu_ratio": float(mu_ratio), "hyst": hyst}
    return HTML.replace("__PAYLOAD__", json.dumps(payload))
