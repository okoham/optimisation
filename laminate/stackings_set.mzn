/*
Basic idea: one array that defines the global sequence of layer orientations, another array that 
defines if a layer is active in a given cell.

Layer  Angle  Cell:  123 ... M
1      45            ---------
2      -45           ---------   
3      0             - --- - -
4      90            -----  
.      ..               ---
.      ..            
N      45            ---------


Decision variables: 

- the angle per layer
- which layer is active in which cell?

Constraints:

1. for each pair of adjacent cells c1 and c2 with MinPlyCount n1 and n2, 
   min(n1, n2) layers have to be continuous. Ply k is continuous if 
   besetzt[c1, k] /\ besetzt[c2, k].
2. Symmetric  
3. Balanced
4. min 10% in any direction
5. Max nc consecutive plies of same orientation
6. outer plies are +/- 45
7. outer plies are continuous, i.e. active in each cell
8. for each cell, and each orientation: actual counts per orientation >= MinPlyCount


Complexity of the problem: how many possible designs are there?
*/


% Stacking
int: MaxPlies = 30; % TODO: kann ich das als decision variable definieren?
int: NOrient = 4; 

set of int: Cells = 1..NCells;

% array[1..NOrient] of int: angles = [ 0, 90, 45, -45 ];
% could use an enum???
set of int: Angles = 1..NOrient;
set of int: Plies = 1..MaxPlies;
% FIXME: hift das?
%set of int: AnglesWithDummy = 0..NOrient; % use 0 as default: ply not used


% Vektor der Orientierungen
array[Plies] of var Angles: stacking;

% MinPlyCounts
array[Cells, Angles] of int: MinPlyCount = 
[| 4, 4, 2, 2 |
   2, 2, 2, 2 |
   3, 2, 2, 2 |
   3, 2, 2, 2 |
   5, 5, 4, 4 |
   3, 4, 2, 2 |
   4, 4, 4, 4 |
   4, 4, 4, 4 |
   3, 4, 2, 2 |];
int: NCells = 9;

/*
% adjancency matrix. upper right triangle only
array[Cells, Cells] of bool: adjacent = 
[| false, true,  false, true,  false, false, false, false, false |
   false, false, true,  false, true,  false, false, false, false |
   false, false, false, false, false, true,  false, false, false |
   false, false, false, false, true,  false, true,  false, false |
   false, false, false, false, false, true,  false, true,  false |
   false, false, false, false, false, false, false, false, true  |
   false, false, false, false, false, false, false, true,  false |
   false, false, false, false, false, false, false, false, true  |
   false, false, false, false, false, false, false, false, false |];

 
% alternativ: paarweise
array[1..12, 1..2] of Cells: adjacent_pairs = 
[| 1, 2 | 2, 3 | 1, 4 | 2, 5 | 3, 6 | 4, 5 | 5, 6 | 4, 7 | 5, 8 | 6, 9 | 7, 8 | 8, 9 |];
*/

% TODO: add some assertions on input data

% array[Cells, Plies] of var bool: besetzt;


%function var int: count_orientation(int: c, int: a) = 
%  sum([1 | i in Plies where stacking[i] = a /\ besetzt[c, i]]);
  
% set based version
function var int: count_orientation_s(int: c, int: a) = 
  sum([1 | i in Plies where stacking[i] = a /\ i in used[c]]);
 
  
% at least MinPlyCount
% TODO: should be >= ???
% TODO: use at_least
%constraint forall (c in Cells, a in Angles) (count_orientation(c, a) >= MinPlyCount[c, a]);
% set based version
constraint forall (c in Cells, a in Angles) (count_orientation_s(c, a) = MinPlyCount[c, a]);
 
% outer plies +/-
% TODO: zusätzliches Constraint, so dass alle komplett leeren Lagen sich aussen sammeln
constraint stacking[1] = 3;
constraint stacking[2] = 4;
constraint stacking[MaxPlies - 1] = 4;
constraint stacking[MaxPlies] = 3;

% das hier funktioniert!
array[Cells] of var set of Plies: used; 
constraint forall (c in Cells) ({1, 2, MaxPlies, MaxPlies-1} subset used[c]);

/*
% constraint 4: min 10% in each direction
function var int: plycount(int: cell) = 
  sum([1 | i in Plies where besetzt[cell, i]]);
  
  
% FIXME: ceil funktioniert nicht mit decision variables?
float: MINRATIO = 0.1;
% dict: nplies -> min plies per orientation
array[Plies] of Plies: nmin = [ceil(MINRATIO * i) | i in Plies ];
array[Cells] of var Plies: thickness = [sum([1 | i in Plies where besetzt[c, i]]) | c in Cells]; %[plycount(c) | c in Cells];
constraint forall (c in Cells, a in Angles) (count_orientation(c, a) >= nmin[thickness[c]]);
*/

% objects covered by the alternatives

% balanced
predicate balanced(int: c) = 
  count_orientation_s(c, 3) = count_orientation_s(c, 4);
  
constraint forall (c in Cells) (balanced(c));

% symmetry
/*
predicate symmetric(int: c) = 
    let {
        array[Plies] of var opt Angles: a = [x | x in used[c]];
        var Plies: n = length(a);
   } in forall (i in 1..n div 2) (a[i] = a[n+1-i]);

constraint forall (c in Cells) (symmetric(c));
*/

% solver; chuffed - schnell, coin-bc - geht so, gecode - langsam

%solve satisfy;
solve minimize sum(c in Cells)(card(used[c])); %sum(c in Cells)(plycount(c));

% output
%output [ "MinPlyCount = ", show(MinPlyCount), "\n" ];
%output [ "stacking = ", show(stacking), "\n" ];
%output [ "besetzt = ", show(bool2int(besetzt)), "\n" ];
output [ "used = ", show(used), "\n" ];
% output [ "plycounts = ", show([sum([1 | i in Plies where besetzt[c, i]]) | c in Cells]), "\n" ];
% output [show([stacking[i] | i in Plies where besetzt[1, i]]), "\n"];
% output [show([stacking[i] | i in Plies where besetzt[2, i]]), "\n"];
% output [show([stacking[i] | i in Plies where besetzt[3, i]]), "\n"];
% output [show([stacking[i] | i in Plies where besetzt[4, i]]), "\n"];
% output [show([stacking[i] | i in Plies where besetzt[5, i]]), "\n"];
%output [show(reverse([stacking[i] | i in Plies where besetzt[1, i]])), "\n"];
%output [show(reverse([x[i] | i in LAYERS where active[i] = true])), "\n"];
