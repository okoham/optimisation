% Use this editor as a MiniZinc scratch book
%include "alldifferent.mzn";

float: PI = 3.141592653589793;

float: El;
float: Et;
float: Glt;
float: nult;
float: tply;

El = 131000;
Et = 11700;
Glt = 4820;
nult = 0.29; 
tply = 0.184;

% Transverse isotropic material
float: nutl = nult * Et/El;
float: Qll = El / (1 - nult*nutl);
float: Qtt = Et / (1 - nult*nutl);
float: Qlt = nult * Et / (1 - nult*nutl);
float: Qss = Glt;

function float: radians(int: theta) = 
  PI*theta/180.0;

% TODO: indices 1..3 oder 0..2 ?
% FIXME: hilft das eigentlich, oder kann ich das Ding nachher nicht zusammenbauen??
/*
function array[1..3, 1..3] of float: q(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));
    float: m2 = m*m;
    float: n2 = n*n;
    float: m3 = m2*m;
    float: n3 = n2*n;
    float: m4 = m3*m;
    float: n4 = n3*n;
    float: Q11 = m4*Qll + n4*Qtt + m2*n2*(2*Qlt + 4*Qss);
    float: Q12 = m2*n2*(Qll + Qtt - 4*Qss) + (m4 + n4)*Qlt;
    float: Q16 = m3*n*Qll + m*n3*Qtt + (m*n3 - m3*n)*Qlt + 2*(m*n3 - m3*n)*Qss;
    float: Q22 = n4*Qll + m4*Qtt + m2*n2*(2*Qlt + 4*Qss);
    float: Q26 = m*n3*Qll + m3*n*Qtt + (m3*n - m*n3)*Qlt + 2*(m3*n - m*n3)*Qss;
    float: Q66 = m2*n2*(Qll + Qtt - 2*Qlt) + (m2 - n2)*Qss;
    array[1..3, 1..3] of float: q_ = [| Q11, Q12, Q16 | Q12, Q22, Q26 | Q16, Q26, Q66 |];
  } in q_;
*/


function float: q11(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));  
    float: q = pow(m, 4)*Qll + pow(n, 4)*Qtt + 2*pow(m, 2)*pow(n, 2)*Qlt + 4*pow(m, 2)*pow(n, 2)*Qss;
  } in q;
  
%    q12 = m**2*n**2*qxx + m**2*n**2*qyy + (m**4 + n**4)*qxy - 4*m**2*n**2*qss
function float: q12(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));  
    float: q = m*m*n*n*Qll + m*m*n*n*Qtt + (m*m*m*m + n*n*n*n)*Qlt - 4*(m*m*n*n)*Qss;
  } in q;
  
%    q16 = m**3*n*qxx - m*n**3*qyy + (m*n**3 - m**3*n)*qxy + 2*(m*n**3 - m**3*n)*qss
function float: q16(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));  
    float: q = m*m*m*n*Qll + m*n*n*n*Qtt + (m*n*n*n - m*m*m*n)*Qlt + 2*(m*n*n*n - m*m*m*n)*Qss;
  } in q;
  
% qt22 = n**4*qxx + m**4*qyy + 2*m**2*n**2*qxy + 4*m**2*n**2*qss
function float: q22(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));  
    float: q = pow(n, 4)*Qll + pow(m, 4)*Qtt + 2*pow(m, 2)*pow(n, 2)*Qlt + 4*pow(m, 2)*pow(n, 2)*Qss;
  } in q;

%    q26 = m*n**3*qxx - m**3*n*qyy + (m**3*n - m*n**3)*qxy + 2*(m**3*n - m*n**3)*qss
function float: q26(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));  
    float: q = m*n*n*n*Qll + m*m*m*n*Qtt + (m*m*m*n - m*n*n*n)*Qlt + 2*(m*m*m*n - m*n*n*n)*Qss;
  } in q;

    
%    q66 = m**2*n**2*qxx + m**2*n**2*qyy - 2*m**2*n**2*qxy + (m**2-n**2)*qss
function float: q66(int: theta) = 
  let {
    float: m = cos(radians(theta));
    float: n = sin(radians(theta));  
    float: q = m*m*n*n*(Qll + Qtt - 2*Qlt) + (m*m - n*n)*Qss;
  } in q;  

% M ... number of orientations
int: M = 4;
set of int: IDX_A = 1..M; % 0, 90, 45, -45
array[IDX_A] of int: ORIENTATIONS = [0, 90, 45, -45];

% N ... number of layers
% COIN-BC is significantly faster than GECODE!!
% however, for 50 plies or more, it slows down.
% use time limit?
% model only half of the laminate, because of symmetry?
int: N = 30;
set of int: LAYERS = 1..N;

array[LAYERS] of var IDX_A: stacking;

array[IDX_A] of float: Q11 = [ q11(a) | a in ORIENTATIONS ];
array[IDX_A] of float: Q12 = [ q12(a) | a in ORIENTATIONS ];
array[IDX_A] of float: Q16 = [ q16(a) | a in ORIENTATIONS ];
array[IDX_A] of float: Q22 = [ q22(a) | a in ORIENTATIONS ];
array[IDX_A] of float: Q26 = [ q26(a) | a in ORIENTATIONS ];
array[IDX_A] of float: Q66 = [ q66(a) | a in ORIENTATIONS ];


% Dicke des Laminats
float: tlam = N*tply;

% z koordinaten der trennebenen zwischen den schichten, von mittelebene aus gemessen
% Annahme, dass Mittellinie die neutrale Faser darstellt
array[1..N+1] of float: z = [ tply*(k - N/2) | k in 1..N+1 ];

% ply counts
%array[1..M] of var int: counts = [sum(a in stacking)(a == i) | i in 1..M];

% Dehnsteifigkeiten
var float: A11 = sum (j in 1..M) (count(stacking, j) * Q11[j] * tply);
var float: A12 = sum (j in 1..M) (count(stacking, j) * Q12[j] * tply);
var float: A16 = sum (j in 1..M) (count(stacking, j) * Q16[j] * tply);
var float: A22 = sum (j in 1..M) (count(stacking, j) * Q22[j] * tply);
var float: A26 = sum (j in 1..M) (count(stacking, j) * Q26[j] * tply);
var float: A66 = sum (j in 1..M) (count(stacking, j) * Q66[j] * tply);

% Biegesteifigkeiten. 
% Annahme, dass Mittelebene die neutrale Ebene darstellt. Ist gültig,
% solange Laminate symmetrisch (und D16 = D26 = 0). Sonst reduzierte 
% Biegesteifigkeiten verwenden.
var float: D11 = sum (k in 1..N) (Q11[stacking[k]] * (pow(z[k+1], 3) - pow(z[k], 3))/3.0); 
var float: D12 = sum (k in 1..N) (Q12[stacking[k]] * (pow(z[k+1], 3) - pow(z[k], 3))/3.0); 
var float: D16 = sum (k in 1..N) (Q16[stacking[k]] * (pow(z[k+1], 3) - pow(z[k], 3))/3.0); 
var float: D22 = sum (k in 1..N) (Q22[stacking[k]] * (pow(z[k+1], 3) - pow(z[k], 3))/3.0); 
var float: D26 = sum (k in 1..N) (Q26[stacking[k]] * (pow(z[k+1], 3) - pow(z[k], 3))/3.0); 
var float: D66 = sum (k in 1..N) (Q66[stacking[k]] * (pow(z[k+1], 3) - pow(z[k], 3))/3.0); 

% Koppelsteifigkeiten
var float: B11 = sum (k in 1..N) (Q11[stacking[k]] * (pow(z[k+1], 2) - pow(z[k], 2))/2.0); 
var float: B12 = sum (k in 1..N) (Q12[stacking[k]] * (pow(z[k+1], 2) - pow(z[k], 2))/2.0); 
var float: B16 = sum (k in 1..N) (Q16[stacking[k]] * (pow(z[k+1], 2) - pow(z[k], 2))/2.0); 
var float: B22 = sum (k in 1..N) (Q22[stacking[k]] * (pow(z[k+1], 2) - pow(z[k], 2))/2.0); 
var float: B26 = sum (k in 1..N) (Q26[stacking[k]] * (pow(z[k+1], 2) - pow(z[k], 2))/2.0); 
var float: B66 = sum (k in 1..N) (Q66[stacking[k]] * (pow(z[k+1], 2) - pow(z[k], 2))/2.0); 


% ------------------
% Design Constraints
% ------------------

% symmetry
constraint forall (i in 1..(N div 2)) (stacking[i] == stacking[N+1-i]);

% TODO: symmetry exception - a +/- pair may occur in the middle of the laminate, if laminate is thick enough

% balanced: Anzahl 45° = Anzahl -45°
% TODO: Konvention machen - "balanced" gilt für +/- Lagen, nicht für 0/90
% TODO: kann man das verallgemeinern, z.B. auf Laminate mit 
%constraint counts[3] == counts[4];

%constraint count(alpha in stacking)(alpha = 3) = count(alpha in stacking)(alpha = 4);
constraint count(stacking, 3) = count(stacking, 4);

% at least 10% in each direction
float: MINRATIO = 0.1;
int: Nmin = ceil(MINRATIO * N);
%constraint forall (i in 1..M) (counts[i] >= ceil(MINRATIO * N));
constraint forall (i in 1..M) (count_leq(stacking, i, Nmin));

% +/-45 on outer plies
% TODO: auch 45/90/-45 zulassen
constraint stacking[1] in {3};
constraint stacking[2] in {4};
constraint stacking[N-1] in {4};
constraint stacking[N] in {3};

% max 2 consecutive
% TODO: verallgemeinern, max K gleiche aufeinanderfolgende Richtungen
constraint forall (k in 1..N-2) (not ((stacking[k] = stacking[k+1]) /\ (stacking[k+1] = stacking[k+2])));

% ------------------------------------------------------------
% compression buckling, all edges ss, kassapoglou eq. 6.7
% N0 = (pi/a)^2 (D11 m^2 + 2(D12 + 2D66) AR^2 + D22 AR^2/m^4)
% AR = length a / width b
% look for the m (integer) that minimises N0
% ------------------------------------------------------------

float: l = 500;
float: b = 200;
float: AR = l/b;
int: MMax = 5; % max. Wellenzahl
var float: ncrit = min (m in 1..MMax) (pow(PI/l, 2)*(D11*m*m + 2*(D12 + 2*D66)*pow(AR, 2) + D22*pow(AR, 4)/(m*m)));

solve maximize ncrit;

output [ "stacking = ", show([ ORIENTATIONS[a] | a in stacking ]), "\n" ] ;
%output [show(counts), "\n"];
output [ "t = ", show(tlam), "\n" ] ;
output [ "ncrit = ", show(ncrit), "\n" ];


% ---------------------------------
% problem: match in-plane stiffness
% ---------------------------------

% constraint A11 > 10000.0;
% constraint A11 < 100000.0;
% constraint A22 > 10000.0;
% constraint A22 < 100000.0;

% solve satisfy;

%output [ "a11 = ", show(A11), "\n" ] ;
