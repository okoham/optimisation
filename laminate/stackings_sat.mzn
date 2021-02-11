/*
Given: a structure consisting of M cells, the minimum number of layers required per cell,
and the information with cells are connected.

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

1. for each cell, and each orientation: actual counts per orientation >= MinPlyCount
2. Symmetric  
3. Balanced, i.e. number of 45° plies equals number of -45° plies
4. min 10% in any direction
5. Max nc consecutive plies of same orientation
6. outer plies are +/- 45
7. outer plies are continuous, i.e. active in each cell
8. +/- always adjacent

9. for each pair of adjacent cells c1 and c2 with MinPlyCount n1 and n2, 
   min(n1, n2) layers have to be continuous. Ply k is continuous if 
   besetzt[c1, k] /\ besetzt[c2, k].

Complexity of the problem: how many possible designs are there?
*/


% Stacking
int: MaxPlies = 21; 
set of int: Plies = 1..MaxPlies;

int: NumDirections = 4; 
% array[1..NumDirections] of int: angles = [ 0, 90, 45, -45 ];
set of int: Directions = 1..NumDirections;

int: NCells = 10;
set of int: Cells = 1..NCells;

% MinPlyCounts
array[Cells, Directions] of int: MinPlyCount = 
[| 4, 4, 2, 2 |
   2, 2, 2, 2 |
   3, 2, 2, 2 |
   3, 2, 4, 4 |
   5, 5, 4, 4 |
   3, 4, 2, 2 |
   4, 4, 4, 4 |
   4, 4, 4, 4 |
   3, 4, 2, 2 |
   2, 2, 4, 4 |];

% Adjacency information
int: NumAdjacent = 15;
array[1..NumAdjacent, 1..2] of Cells: adjacent_pairs = 
  [| 1, 2 | 2, 3 | 1, 4 | 2, 5 | 3, 6 | 4, 5 | 5, 6 | 4, 7 | 5, 8 | 6, 9 | 7, 8 | 8, 9 | 1, 10 | 2, 10 | 3, 10 |];

% Main decision variable: global stacking sequence, and activation table
%array[Plies] of var Directions: global_sequence;
%array[Cells, Plies] of var bool: active;

array[Cells, Plies] of var bool: active;
array[Plies, Directions] of var bool: X; % global_sequence

predicate atmostone(array[int] of var bool:x) =
  forall(i,j in index_set(x) where i < j)((not x[i] \/ not x[j]));

predicate exactlyone(array[int] of var bool:x) =
  atmostone(x) /\ exists(x);

% there has to be one direction per ply
constraint forall(p in Plies)(exactlyone(X[p, Directions]));

% Ein paar Hilfsvariable
array[Cells] of Plies: MinThickness = [sum(row(MinPlyCount, c)) | c in Cells];

function var int: counts_per_angle(Cells: c, Directions: a) = sum(p in Plies)(bool2int(active[c, p] /\ X[p, a]));
function var int: thickness(Cells: c) = sum(p in Plies)(bool2int(active[c, p]));

% ---------------------------------------
% Constraint 1: MinPlyCount per direction
% ---------------------------------------
constraint forall (c in Cells, a in Directions) (counts_per_angle(c, a) >= MinPlyCount[c, a]);


% ----------------------
% Constraint 2: symmetry
% ----------------------
predicate symmetric(int: c) = 
    let {
        array[Plies] of var opt Directions: seq = [index2num(X[p, Directions]) | p in Plies where active[c, p]];
        array[Plies] of var opt Directions: revseq = reverse(seq);
        var Plies: n = length(seq);
   } in 
   forall (i in 1..n div 2) (seq[i] == revseq[i]); 
   
constraint forall (c in Cells) (symmetric(c));

% ----------------------
% Constraint 3: balanced
% ----------------------

constraint forall (c in Cells) (
  counts_per_angle(c, 3) = counts_per_angle(c, 4)
);

/* 
----------------------------------------
Constraint 4: min 10% in each direction
-----------------------------------------

Note ceil function does not work with var.
Therefore, create an auxiliary array nmin that maps from total ply count 
to minimum ply count per orientation
*/

float: MINRATIO = 0.1;
array[Plies] of Plies: nmin = [ceil(MINRATIO * n) | n in Plies ];
constraint forall (c in Cells, a in Directions) (counts_per_angle(c, a) >= nmin[thickness(c)]);


% --------------------------------------------------
% Constraint 5: not more than Ncmax contiguous plies
% --------------------------------------------------
% FIXME: noch falsch???


constraint forall (c in Cells) (
    let {
        array[Plies] of var opt int: seq = [p | p in Plies where active[c, p]];
        var Plies: n = length(seq);
   } in 
   forall (k in 1..n-2, a in Directions) (
     sum(i in k..k+2)(bool2int(X[i, a])) < 3
   )
);


% -----------------------------
% Constraint 6: outer plies +/-
% -----------------------------
constraint X[1, 3];
constraint X[2, 4];
constraint X[MaxPlies-1, 4];
constraint X[MaxPlies, 3];


% ----------------------------------------------------
% Constraint 7: 2 outer plies are present in all cells
% ----------------------------------------------------
set of Plies: outer_plies = {1, 2, MaxPlies-1, MaxPlies};
constraint forall (c in Cells, p in outer_plies) (active[c, p]);


% --------------------------------------------
% Constraint 8: +/- should be grouped together
% --------------------------------------------

predicate grouped(Cells: c) = 
    let {
        array[Plies] of var opt Plies: seq = [p | p in Plies where active[c, p]];
        var Plies: n = length(seq);
   } in 
   forall (p in 2..n-1) (
     (X[p, 3] -> (X[p-1, 4] \/ X[p+1, 4])) /\ 
     (X[p, 4] -> (X[p+1, 3] \/ X[p+1, 3]))
   ); 

constraint forall(c in Cells)(grouped(c));

/* ------------------------------------------------
 Constraint 9: continuity between neighboring cells
 --------------------------------------------------
 
 all plies, except dropped ones, must be continuous
*/
constraint forall (k in 1..NumAdjacent) (
  let {
    Cells: c1 = adjacent_pairs[k, 1]; 
    Cells: c2 = adjacent_pairs[k, 2];
    var 0..MaxPlies: common = sum(p in Plies)(bool2int(active[c1, p] /\ active[c2, p]));
  } in 
  common >= min(thickness(c1), thickness(c2))
);


% ---------------
% solver 
% ---------------
% chuffed - schnell, coin-bc - geht so, gecode - langsam

var int: vtotal = sum(c in Cells)(thickness(c));
% solve satisfy;
solve minimize vtotal; 

% -------
% output
% -------

function var int: index2num(array[Directions] of var bool: x) =
  let {array[Directions] of var int: b = bool2int(x)} in 
  sum(i in Directions)(i * b[i]);

output ["\(c)  " ++ show([index2num(row(X, p)) | p in Plies where active[c, p]]) ++ "\n" | c in Cells];
output ["Desired volume: \(sum(MinThickness)), realised volume: \(vtotal)\n"];
