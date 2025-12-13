var 
    obs_capu obs_unemp obs_spread obs_stockret obs_impinfl obs_dg obs_r10
    ewma epinfma  
    zcapf rkf kf pkf cf invef yf labf wf rrf
    mc zcap rk k pk c inve y lab pinf w r 
    a b g qs ms spinf sw 
    kpf kp;

varexo ea eb eg eqs em epinf ew;

parameters 
    curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
    czcap csadjcost ctou csigma chabb ccs cinvs cfc cindw cprobw cindp cprobp csigl 
    clandaw crdpi crpi crdy cry crr 
    crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow 
    ctrend cg
    scale_capu scale_unemp scale_spread scale_stockret scale_impinfl scale_dg scale_r10;

ctou = 0.025;
clandaw = 1.5;
cg = 0.18;
curvp = 10;
curvw = 10;

calfa = 0.24;
cfc = 1.5;
cgy = 0.51;

csadjcost = 6.0;
csigl = 2.0;
czcap = 0.5;

csigma = 1.5;
chabb = 0.7;
cprobw = 0.66;
cprobp = 0.66;
cindw = 0.5;
cindp = 0.5;
crpi = 1.7;
crr = 0.8;
cry = 0.1;
crdy = 0.1;

crhoa = 0.95;
crhob = 0.85;
crhog = 0.1;
crhols = 0.99;
crhoqs = 0.3;
crhoas = 1; 
crhoms = 0.5;
crhopinf = 0.7;
crhow = 0.95;
cmap = 0.5;
cmaw = 0.5;

constelab = 0;
constepinf = 0.7;
constebeta = 0.25;
ctrend = 0.4;

ccs = 0; cinvs = 0; crdpi = 0;

scale_capu = 1.0;
scale_unemp = 1.0;
scale_spread = 1.0;
scale_stockret = 1.0;
scale_impinfl = 1.0;
scale_dg = 1.0;
scale_r10 = 1.0;

model(linear); 

#cpie = 1 + constepinf/100;
#cgamma = 1 + ctrend/100;
#cbeta = 1/(1 + constebeta/100);
#clandap = cfc;
#cbetabar = cbeta*cgamma^(-csigma);
#cr = cpie/(cbeta*cgamma^(-csigma));
#crk = (cbeta^(-1))*(cgamma^csigma) - (1-ctou);
#cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
#cikbar = (1-(1-ctou)/cgamma);
#cik = (1-(1-ctou)/cgamma)*cgamma;
#clk = ((1-calfa)/calfa)*(crk/cw);
#cky = cfc*(clk)^(calfa-1);
#ciy = cik*cky;
#ccy = 1-cg-cik*cky;
#crkky = crk*cky;
#cwhlc = (1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
#cwly = 1-crk*cky;
#conster = (cr-1)*100;

0*(1-calfa)*a + 1*a = calfa*rkf + (1-calfa)*(wf);
zcapf = (1/(czcap/(1-czcap)))* rkf;
rkf = (wf) + labf - kf;
kf = kpf(-1) + zcapf;
invef = (1/(1+cbetabar*cgamma))*(invef(-1) + cbetabar*cgamma*invef(1) + (1/(cgamma^2*csadjcost))*pkf) + qs;
pkf = -rrf - 0*b + (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rkf(1) + ((1-ctou)/(crk+(1-ctou)))*pkf(1);
cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) + ((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b;
yf = ccy*cf + ciy*invef + g + crkky*zcapf;
yf = cfc*(calfa*kf + (1-calfa)*labf + a);
wf = csigl*labf + (1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1);
kpf = (1-cikbar)*kpf(-1) + (cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs;

mc = calfa*rk + (1-calfa)*(w) - 1*a - 0*(1-calfa)*a;
zcap = (1/(czcap/(1-czcap)))* rk;
rk = w + lab - k;
k = kp(-1) + zcap;
inve = (1/(1+cbetabar*cgamma))*(inve(-1) + cbetabar*cgamma*inve(1) + (1/(cgamma^2*csadjcost))*pk) + qs;
pk = -r + pinf(1) - 0*b + (1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) + ((1-ctou)/(crk+(1-ctou)))*pk(1);
c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) + ((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) + b;
y = ccy*c + ciy*inve + g + crkky*zcap;
y = cfc*(calfa*k + (1-calfa)*lab + a);
pinf = (1/(1+cbetabar*cgamma*cindp))*(cbetabar*cgamma*pinf(1) + cindp*pinf(-1) + ((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)) + spinf;
w = (1/(1+cbetabar*cgamma))*w(-1) + (cbetabar*cgamma/(1+cbetabar*cgamma))*w(1) + (cindw/(1+cbetabar*cgamma))*pinf(-1) - (1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf + (cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1) + (1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*(csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) - w) + 1*sw;
r = crpi*(1-crr)*pinf + cry*(1-crr)*(y-yf) + crdy*(y-yf-(y(-1)-yf(-1))) + crr*r(-1) + ms;

a = crhoa*a(-1) + ea;
b = crhob*b(-1) + eb;
g = crhog*(g(-1)) + eg + cgy*ea;
qs = crhoqs*qs(-1) + eqs;
ms = crhoms*ms(-1) + em;
spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
epinfma = epinf;
sw = crhow*sw(-1) + ewma - cmaw*ewma(-1);
ewma = ew;
kp = (1-cikbar)*kp(-1) + cikbar*inve + cikbar*cgamma^2*csadjcost*qs;


obs_capu = scale_capu * zcap;

obs_unemp = -scale_unemp * lab + constelab;

obs_spread = scale_spread * b;

obs_stockret = scale_stockret * (pk - pk(-1) + inve - inve(-1));

obs_impinfl = scale_impinfl * pinf;

obs_dg = scale_dg * (g - g(-1));

obs_r10 = scale_r10 * r;

end; 

steady_state_model;
obs_capu = 0;
obs_unemp = constelab;
obs_spread = 0;
obs_stockret = 0;
obs_impinfl = 0;
obs_dg = 0;
obs_r10 = 0;
end;

shocks;
var ea; stderr 0.5;
var eb; stderr 0.2;
var eg; stderr 0.5;
var eqs; stderr 0.5;
var em; stderr 0.2;
var epinf; stderr 0.5;
var ew; stderr 0.2;
end;

estimated_params;

stderr ea, 0.5, 0.01, 5, INV_GAMMA_PDF, 0.5, 2;
stderr eb, 0.2, 0.01, 3, INV_GAMMA_PDF, 0.2, 2;
stderr eg, 0.5, 0.01, 3, INV_GAMMA_PDF, 0.5, 2;
stderr eqs, 0.5, 0.01, 5, INV_GAMMA_PDF, 0.5, 2;
stderr em, 0.2, 0.01, 3, INV_GAMMA_PDF, 0.2, 2;
stderr epinf, 0.5, 0.01, 5, INV_GAMMA_PDF, 0.5, 2;
stderr ew, 0.2, 0.01, 3, INV_GAMMA_PDF, 0.2, 2;

crhoa, 0.9, 0.5, 0.999, BETA_PDF, 0.85, 0.1;

crhob, 0.8, 0.3, 0.999, BETA_PDF, 0.8, 0.1;

crhog, 0.2, 0.01, 0.999, BETA_PDF, 0.3, 0.2;

crhoqs, 0.3, 0.01, 0.999, BETA_PDF, 0.4, 0.2;

crhoms, 0.3, 0.01, 0.999, BETA_PDF, 0.4, 0.2;

crhopinf, 0.6, 0.01, 0.999, BETA_PDF, 0.6, 0.15;

crhow, 0.9, 0.5, 0.999, BETA_PDF, 0.85, 0.1;

cmap, 0.5, 0.01, 0.999, BETA_PDF, 0.5, 0.2;
cmaw, 0.5, 0.01, 0.999, BETA_PDF, 0.5, 0.2;


csigma, 1.5, 0.5, 4, NORMAL_PDF, 1.5, 0.5;

chabb, 0.65, 0.1, 0.95, BETA_PDF, 0.65, 0.1;

cprobw, 0.66, 0.3, 0.95, BETA_PDF, 0.66, 0.1;

cprobp, 0.66, 0.3, 0.95, BETA_PDF, 0.66, 0.1;

cindw, 0.5, 0.01, 0.99, BETA_PDF, 0.5, 0.15;
cindp, 0.5, 0.01, 0.99, BETA_PDF, 0.5, 0.15;

crpi, 1.7, 1.0, 3.0, NORMAL_PDF, 1.7, 0.3;

crr, 0.75, 0.5, 0.97, BETA_PDF, 0.75, 0.1;

cry, 0.1, 0.001, 0.5, NORMAL_PDF, 0.1, 0.05;
crdy, 0.1, 0.001, 0.5, NORMAL_PDF, 0.1, 0.05;

scale_capu, 1.0, 0.1, 10, NORMAL_PDF, 1.0, 0.5;
scale_unemp, 1.0, 0.1, 10, NORMAL_PDF, 1.0, 0.5;
scale_spread, 1.0, 0.1, 10, NORMAL_PDF, 1.0, 0.5;
scale_stockret, 1.0, 0.1, 20, NORMAL_PDF, 2.0, 1.0;
scale_impinfl, 1.0, 0.1, 30, NORMAL_PDF, 5.0, 2.0;
scale_dg, 1.0, 0.1, 10, NORMAL_PDF, 1.0, 0.5;
scale_r10, 1.0, 0.1, 10, NORMAL_PDF, 1.0, 0.5;

constelab, 0, -5, 5, NORMAL_PDF, 0, 1;

end;

estimated_params_init(use_calibration);
end;

varobs obs_capu obs_unemp obs_spread obs_stockret obs_impinfl obs_dg obs_r10;

estimation(
    datafile = data_alternative,
    first_obs = 1,
    prefilter = 0,
    lik_init = 2,
    mh_replic = 0,
    mode_compute = 4,
    mode_check
) obs_capu obs_unemp obs_spread obs_stockret obs_impinfl obs_dg obs_r10;
