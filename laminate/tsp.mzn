% Use this editor as a MiniZinc scratch book

/*
# City X Y
Berlin 5251 1340
Frankfurt 5011 864
Leipzig 5133 1237
Heidelberg 4941 867
Karlsruhe 4901 840
Hamburg 5356 998
Bayreuth 4993 1159
Trier 4974 668
Hannover 5237 972
Stuttgart 4874 909
Passau 4856 1344
Augsburg 4833 1089
Koblenz 5033 759
Dortmund 5148 741
Bochum 5145 728
Duisburg 5142 679
Wuppertal 5124 715
Essen 5145 701
Jena 5093 1158

set V := { read "tsp.dat" as "<1s>" comment "#" };
set E := { <i , j> in V ∗ V with i < j };
set P[] := powerset (V);
set K   := indexset (P);

param px[V] := read "tsp.dat" as "<1s> 2n" comment "#";
param py[V] := read "tsp.dat" as "<1s> 3n" comment "#";

defnumb dist(a, b) := sqrt((px[a]−px[b])^2 + (py[a]−py[b])^ 2);

var x[E] binary;

minimize cost: sum <i,j> in E: dist(i,j) ∗ x[i,j];

subto two_connected: forall <v> in V do
  (sum <v,j> in E : x[v,j]) + (sum <i , v> in E: x[i,v]) == 2;

subto no_subtour: 
  forall <k> in K with
    card (P[k]) > 2 and card (P[k]) < card(V) − 2 do
      sum <i,j> in E with <i> in P[k] and <j> in P[k]:
        x[i,j] <= card (P[k]) − 1;
*/

include "alldifferent.mzn";

enum cities = {Berlin, Frankfurt, Leipzig, Heidelberg, Karlsruhe, Hamburg, Bayreuth, Trier,
               Hannover, Stuttgart, Passau, Augsburg, Koblenz, Dortmund, Bochum, Duisburg, 
               Wuppertal, Essen, Jena};
               
int: N = length(cities);
               
array[1..N, 1..2] of int: coordinates = [| 5251, 1340
                                       | 5011, 864
                                       | 5133, 1237
                                       | 4941, 867
                                       | 4901, 840
                                       | 5356, 998
                                       | 4993, 1159
                                       | 4974, 668
                                       | 5237, 972
                                       | 4874, 909
                                       | 4856, 1344
                                       | 4833, 1089
                                       | 5033, 759
                                       | 5148, 741
                                       | 5145, 728
                                       | 5142, 679
                                       | 5124, 715
                                       | 5145, 701
                                       | 5093, 1158 |];
                                       
array[1..N, 1..N] of int: distances;

                                       
array[1..N] of var 1..N: route;

constraint alldifferent(route);

var int: laenge = sum(i in 1..N-1)((coordinates[route[i+1],1] - coordinates[route[i],1])^2 + (coordinates[route[i+1],2] - coordinates[route[i],2])^2) + (coordinates[route[1],1] - coordinates[route[N],1])^2 + (coordinates[route[1],2] - coordinates[route[N],2])^2;

solve minimize laenge;

output ["\(laenge)", "\n"];