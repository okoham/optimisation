% Use this editor as a MiniZinc scratch book
include "alldifferent.mzn";
include "count.mzn";

set of int: werte = 1..12;

int: Nmax = 8;
set of int: LAYERS = 1..Nmax;
array[LAYERS] of var bool: active;
array[LAYERS] of var werte: x;
%var int: nonzero;


constraint sum([x[i] | i in LAYERS where active[i] = true]) = 28;
constraint sum([x[i] | i in LAYERS where active[i] = false]) = 27;
%constraint count(i in LAYERS)(active[i]) = nonzero;
%constraint active[2] = false;
constraint all_different(x);
constraint forall(i in 1..Nmax-1)(x[i] < x[i+1]);
%constraint count(active, true) > 3;


solve satisfy;

output [show(active), "\n"];
output [show(x), "\n"];
output [show(array2set([x[i] | i in LAYERS where active[i]])), "\n"];
output [show([x[i] | i in LAYERS where active[i] = true]), "\n"];
output [show(reverse([x[i] | i in LAYERS where active[i] = true])), "\n"];
output [show(reverse([x[i] | i in {1, 3, 5, 7}])), "\n"];

%output [show(index_set([active[i] | i in LAYERS where active[i] = true])), "\n"];