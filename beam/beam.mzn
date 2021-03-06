% Use this editor as a MiniZinc scratch book


enum MATERIALS = {AL7010, AL2198, TI64};
array[MATERIALS] of float: E = [71000.0, 76000.0, 110000.0];
array[MATERIALS] of float: nu = [0.33, 0.33, 0.33];
array[MATERIALS] of float: ftu = [515.0, 495.0, 900.0];
array[MATERIALS] of float: fty = [455.0, 430.0, 800.0];
array[MATERIALS] of float: fcy = [440.0, 415.0, 780.0];
array[MATERIALS] of float: fsu = [295.0, 270.0, 520.0];
array[MATERIALS] of float: rho = [2.82e-6, 2.69e-6, 4.5e-6];
array[MATERIALS] of float: cmat = [5.0, 10.0, 60.0];
array[MATERIALS] of float: cprod = [5.0, 10.0, 60.0];

float: L;

var MATERIALS: mat;
var 20.0..300.0: h;
var 1.6..10.0: tw;
var 1.6..50.0: blf;
var 1.6..10.0: tlf;
var 1.6..50.0: buf;
var 1.6..10.0: tuf;

var float: a_uf;
constraint a_uf = buf * tuf;
var float: a_lf;
constraint a_lf = blf * tlf;
var float: a_w;
constraint a_w = (h - tlf - tuf)*tw;
var float: area;
constraint area = a_w + a_lf + a_uf;
var float: cg_w;
constraint cg_w = (tlf + h - tuf)/2.0;
var float: cg_uf;
constraint cg_uf = h - tuf/2.0;
var float: cg_lf;
constraint cg_lf = tlf/2.0;
var float: cg;
constraint cg = (cg_uf*a_uf + cg_w*a_w + cg_lf*a_lf) / area;
var float: iyy;
constraint iyy = buf*tuf^3/12.0 + a_uf*(cg_uf - cg)^2 + blf*tlf^3/12.0 + a_lf*(cg_lf - cg)^2 + tw*(h - tuf - tlf)^3/12.0 + a_w*(cg_w - cg)^2;

% var float: it = 
var float: mass;
constraint mass = rho[mat] * L * area;

solve minimize mass;
