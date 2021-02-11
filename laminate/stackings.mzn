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

int: NCells = 9;
set of int: Cells = 1..NCells;

% MinPlyCounts
array[Cells, Directions] of int: MinPlyCount = 
[| 4, 4, 2, 2 |
   2, 2, 2, 2 |
   3, 2, 2, 2 |
   3, 2, 2, 4 |
   5, 5, 4, 4 |
   3, 4, 2, 2 |
   4, 4, 4, 4 |
   4, 4, 4, 4 |
   3, 4, 2, 2 |];

% Adjacency information
int: NumAdjacent = 12;
array[1..NumAdjacent, 1..2] of Cells: adjacent_pairs = 
  [| 1, 2 | 2, 3 | 1, 4 | 2, 5 | 3, 6 | 4, 5 | 5, 6 | 4, 7 | 5, 8 | 6, 9 | 7, 8 | 8, 9 |];

% Main decision variable: global stacking sequence, and activation table
array[Plies] of var Directions: global_sequence;
array[Cells, Plies] of var bool: active;

% Ein paar Hilfsvariable
array[Cells] of Plies: MinThickness = [sum(row(MinPlyCount, c)) | c in Cells];
array[Cells] of var Plies: thickness = [count(i in Plies)(active[c, i]) | c in Cells];
array[Cells, Directions] of var Plies: counts_per_angle;
constraint forall(c in Cells, a in Directions)(
    counts_per_angle[c, a] = count(i in Plies)(global_sequence[i] = a /\ active[c, i])
);


% ----------------------------
% Constraint 1: MinPlyCount
% ----------------------------
constraint forall (c in Cells, a in Directions) (counts_per_angle[c, a] >= MinPlyCount[c, a]);


% ----------------------
% Constraint 2: symmetry
% ----------------------
predicate symmetric(int: c) = 
    let {
        array[Plies] of var opt Directions: a = [global_sequence[j] | j in Plies where active[c, j]];
        var int: n = length(a);
   } in 
   forall (i in 1..n div 2) (a[i] = a[n+1-i]);
   
constraint forall (c in Cells) (symmetric(c));

% ----------------------
% Constraint 3: balanced
% ----------------------
predicate balanced(int: c) = 
  counts_per_angle[c, 3] = counts_per_angle[c, 4];

constraint forall (c in Cells) (balanced(c));


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
constraint forall (c in Cells, a in Directions) (counts_per_angle[c, a] >= nmin[thickness[c]]);


% --------------------------------------------------
% Constraint 5: not more than Ncmax contiguous plies
% --------------------------------------------------
% FIXME: noch falsch!

constraint forall (c in Cells) (
  let {
    array[Plies] of var opt Directions: seq = [global_sequence[p] | p in Plies where active[c, p]];
    var Plies: n = length(seq);
  } in 
  forall (a in Directions, i in 1..n-2) (
       sum([1 | j in i..i+2 where seq[j] == a]) < 3
  ) 
);

% -----------------------------
% Constraint 6: outer plies +/-
% -----------------------------
constraint global_sequence[1] = 3;
constraint global_sequence[2] = 4;
constraint global_sequence[MaxPlies - 1] = 4;
constraint global_sequence[MaxPlies] = 3;

% ----------------------------------------------------
% Constraint 7: 2 outer plies are present in all cells
% ----------------------------------------------------
set of Plies: outer_plies = {1, 2, MaxPlies-1, MaxPlies};
constraint forall (c in Cells, p in outer_plies) (active[c, p]);


% --------------------------------------------
% Constraint 8: +/- should be grouped together
% --------------------------------------------
% whenever there is a 45 ply, there should be an adjacent -45 ply
constraint forall (c in Cells) (
  let {
    array[Plies] of var opt Directions: seq = [global_sequence[p] | p in Plies where active[c, p]];
    var Plies: n = length(seq);
  } in 
  forall (i in 2..n-1) (
     (seq[i] == 3 -> (seq[i-1] == 4 \/ seq[i+1] == 4)) /\ 
     (seq[i] == 4 -> (seq[i-1] == 3 \/ seq[i+1] == 3))
   ) 
%   forall (i in 2..n-2 where seq[i] == 3) (seq[i-1] == 4 \/ seq[i+1] == 4) /\
%   forall (i in 2..n-2 where seq[i] == 4) (seq[i-1] == 3 \/ seq[i+1] == 3)
);

/* ------------------------------------------------
 Constraint 9: continuity between neighboring cells
 --------------------------------------------------
 
 all plies, except dropped ones, must be continuous
*/
constraint forall (k in 1..NumAdjacent) (
  let {
    Cells: c1 = adjacent_pairs[k, 1]; 
    Cells: c2 = adjacent_pairs[k, 2];
    var set of Plies: s1 = {global_sequence[j] | j in Plies where active[c1, j]};
    var set of Plies: s2 = {global_sequence[j] | j in Plies where active[c2, j]};
  } in 
  card(s1 intersect s2) >= min(card(s1), card(s2))
);


% ---------------
% solver 
% ---------------
% chuffed - schnell, coin-bc - geht so, gecode - langsam

%solve satisfy;
solve minimize sum(thickness); 

% -------
% output
% -------

%array[0..1] of string: lala = [" ", "X"];
%output [ "Global global_sequence sequence = ", show(global_sequence), "\n" ];
%output [ "Min plycounts    = ", show([MinThickness[c] | c in Cells]), "\n" ];
%output [ "Actual plycounts = ", show([thickness[c] | c in Cells]), "\n" ];
output [show_int(2, c) ++ ": " ++ show([global_sequence[i] | i in Plies where active[c, i]]) ++ "\n" | c in Cells] ++ [ "\n" ];
% output [show_int(3, p) ++ ": \(global_sequence[p]) " ++ show([bool2int(active[c, p]) | c in Cells]) ++ "\n" | p in Plies] ++ ["\n"];
% output ["Ply counts\n"] ++ ["\(a):     " ++ show([counts_per_angle[c, a] | c in Cells]) ++ "\n" | a in Directions];
output ["Desired volume: " ++ show_int(4, sum(MinThickness)) ++ ", Achieved volume: " ++ show_int(4, sum(thickness)) ++ "\n"];
