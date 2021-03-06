/*
Given: a structure consisting of M cells, the minimum number of layers required per cell,
and the information with cells are connected.

Basic idea: use the array based set representation of minizinc course week 2,
bounded cardinality set

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
array[Plies] of var Directions: GO;

% active: sorted descending
array[Cells, Plies] of var 0..MaxPlies: active;
constraint forall(c in Cells, i in 1..MaxPlies-1)
                 (active[c,i] >= (active[c,i] != 0) + active[c,i+1]);

% Ein paar Hilfsvariable
array[Cells] of Plies: MinThickness = [sum(row(MinPlyCount, c)) | c in Cells];
array[Cells] of var Plies: thickness =  [sum(i in Plies)(active[c,i] > 0) | c in Cells]; 
array[Cells, Directions] of var Plies: counts_per_angle;
constraint forall(c in Cells, a in Directions)
                 (counts_per_angle[c, a] = sum(i in Plies)(active[c,i] > 0 /\ GO[active[c,i]] = a));


% ----------------------------
% Constraint 1: MinPlyCount
% ----------------------------
constraint forall (c in Cells, a in Directions)
                  (counts_per_angle[c, a] >= MinPlyCount[c, a]);


% ----------------------
% Constraint 2: symmetry
% ----------------------

constraint forall(c in Cells) (
           let {
             var 4..MaxPlies: t = thickness[c];
           } in 
           forall(i in 1..t div 2) 
                 (GO[active[c,i]] = GO[active[c, t+1-i]]));

% ----------------------
% Constraint 3: balanced
% ----------------------
% predicate balanced(int: c) = 
%   counts_per_angle[c, 3] = counts_per_angle[c, 4];

constraint forall (c in Cells) (counts_per_angle[c, 3] = counts_per_angle[c, 4]);


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
constraint forall (c in Cells, a in Directions) 
                  (counts_per_angle[c, a] >= nmin[thickness[c]]);
                  
% --------------------------------------------------
% Constraint 5: not more than Ncmax contiguous plies
% --------------------------------------------------
% FIXME: noch falsch?

int: MaxCont = 2;

constraint forall(c in Cells, i in 1..thickness[c] - MaxCont) (
           sum(j in i..i+MaxCont)(GO[active[c,i]] = GO[active[c,j]]) <= MaxCont
);

% -----------------------------
% Constraint 6: outer plies +/-
% -----------------------------
constraint GO[1] = 3;
constraint GO[2] = 4;
constraint GO[MaxPlies - 1] = 4;
constraint GO[MaxPlies] = 3;

% ----------------------------------------------------
% Constraint 7: 2 outer plies are present in all cells
% ----------------------------------------------------
constraint forall (c in Cells) (active[c, 1] = MaxPlies);
constraint forall (c in Cells) (active[c, 2] = MaxPlies-1);
constraint forall (c in Cells) (active[c, thickness[c]-1] = 2);
constraint forall (c in Cells) (active[c, thickness[c]] = 1);


% --------------------------------------------
% Constraint 8: +/- should be grouped together
% --------------------------------------------
% whenever there is a 45 ply, there should be an adjacent -45 ply


constraint forall (c in Cells) (
  forall (i in 2..thickness[c]-1) (
     (GO[active[c,i]] == 3 -> (GO[active[c,i-1]] == 4 \/ GO[active[c,i+1]] == 4)) /\ 
     (GO[active[c,i]] == 4 -> (GO[active[c,i-1]] == 3 \/ GO[active[c,i+1]] == 3))
   ) 
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
    var set of Plies: s1 = {j | j in Plies where active[c1, j] > 0};
    var set of Plies: s2 = {j | j in Plies where active[c2, j] > 0};
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
%output [ "Global GO sequence = ", show(GO), "\n" ];
%output [ "Min plycounts    = ", show([MinThickness[c] | c in Cells]), "\n" ];
%output [ "Actual plycounts = ", show([thickness[c] | c in Cells]), "\n" ];
output [show_int(2, c) ++ ": " ++ show([GO[active[c,i]] | i in 1..thickness[c]]) ++ "\n" | c in Cells] ++ [ "\n" ] ++ ["Lower Bound: " ++ show_int(4, sum(MinThickness)) ++ ", Achieved volume: " ++ show_int(4, sum(thickness)) ++ "\n"]; 
% ++ [show(row(active, c))  ++"\n" | c in Cells] ++ ["\n\n" ++ show(GO)];
% output [show_int(3, p) ++ ": \(GO[p]) " ++ show([bool2int(active[c, p]) | c in Cells]) ++ "\n" | p in Plies] ++ ["\n"];
% output ["Ply counts\n"] ++ ["\(a):     " ++ show([counts_per_angle[c, a] | c in Cells]) ++ "\n" | a in Directions];
