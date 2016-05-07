(* :Title: BackwardAlgorithm, Version 1.0, 2016*)

(* :Author: Kacper Pluta, kacper.pluta@esiee.fr
	    Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France 
*)

(* :Summary:

  A package which implements algorithm which allows to find intersection of lattices of preimages of
  elements of a cyclic group in the remainder range. In details, it allows to predict where in the
  transformed space non-injective and non-surjective points occur under 2D digitized rigid motions.
  The algorithm was introduced in:

@inproceedings{PLUTA:DGCI:2016,
  author = {Pluta, Kacper and Romon, Pascal and Kenmochi, Yukiko and Passat, Nicolas},
  booktitle = {{Discrete Geometry for Computer Imagery}},
  pages = {359–371},
  title = {{Bijective rigid motions of the 2D Cartesian grid}},
  url = {http://link.springer.com/chapter/10.1007/978-3-319-32360-2_28; 
      https://hal-upec-upem.archives-ouvertes.fr/hal-01275598/file/article.pdf},
  year = {2016}
}

*)

(* :Context: BackwardAlgorithm` *)

(* :Package Version: 1.0 *)

(* : Copyright: Kacper Pluta, 2016
  All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL Kacper Pluta BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*)


(* :History:
  Created by Kacper Pluta at University Paris-Est, 2016
*)

(* :Keywords: 

  2D digitized rigid motions, non-injective, non-surjective
*)


(* :Mathematica Version: 10.0.0 *)

(* :Limitations: 

  For this moment this package supports only some particular generators of primitive Pythagorean
  triples. Also translation is not supported.

*)


BeginPackage["BackwardAlgorithm`"]

GetNonInjectiveRegion::usage = "GetNonInjectiveRegion[index, p, q] returns Rectange[]";
GetNonSurjectiveRegion::usage = "GetNonSurjectiveRegion[index, p, q] returns Rectange[]";
FindGroupMembersInRegions::usage = "FindGroupMembersInRegions[p, q, F] returns a set of integer
points";
IntersectionSetLatticesNonInjective::usage = "IntersectionSetLatticesNonInjective[p, q, set] returns
subset of the set.";
IntersectionSetLatticesNonSurjective::usage = "IntersectionSetLatticesNonSurjective[p, q, set] returns
subset of the set.";

Begin["`Private`"]

(*

Procedure: GetNonInjectiveFrame

Summary:
  This procedure tests if generators of primitive Pythagorean triple are correct. If not it throws
  exceptions.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.

*)
TestPythagoreanTripleGenerators[p_, q_] := (
  If[index < 1 || index > 4, Throw["Frame index out of the range:[1,4]"]];
  If[GCD[p, q] != 1 || EvenQ[p - q], Throw["Pythagorean rotation generates are incorrect!"]];
)

(*

Procedure: GetNonInjectiveRegion

Summary:
  This procedure returns a non-injective region in the remainder range.

Parameters:
  index -- an index of given non-injective region where 1 - up, 2 - right, 3 - down and 4 left. For
  details look into a paper cited in :Summary: of this package.  
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.

Output:
  A non-injective region of a given index.

*)
GetNonInjectiveRegion[index_, p_, q_] := Module[{a, b, c},
  TestPythagoreanTripleGenerators[p, q];
  (*primitive Pythagorean triple*)
  b = 2 p q; a = p^2 - q^2; c = p^2 + q^2;
  (*up*)
  If[index == 1, Return[Rectangle[{b/c - 1/2, -1/2}, {1/2, 1/2 - a/c}]]];
  (*right*)
  If[index == 2, Return[Rectangle[{-1/2, -1/2}, {1/2 - a/c, 1/2 - b/c}]]];
  (*down*)
  If[index == 3, Return[Rectangle[{-1/2, a/c - 1/2}, {1/2 - b/c, 1/2}]]];
  (*left*)
  If[index == 4, Return[Rectangle[{a/c - 1/2, b/c - 1/2}, {1/2, 1/2}]]]
]; (* end of GetNonInjectiveRegion *)


(*

Procedure: GetNonSurjectiveRegion

Summary:
  This procedure returns a non-surjective region in the remainder range.

Parameters:
  index -- an index of given non-surjective region where 1 - up, 2 - right, 3 - down and 4 left. For
  details look into a paper cited in :Summary: of this package.  
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.

Output:
  A non-surjective region of a given index.

*)
GetNonSurjectiveRegion[index_, p_, q_] := Module[{a, b, c},
  TestPythagoreanTripleGenerators[p, q];
  b = 2 p q; a = p^2 - q^2; c = p^2 + q^2;
  (*up*)
  If[index == 1, Return[Rectangle[{1/2 - a/c, 3/2 - a/c - b/c}, {b/c - 1/2, 1/2}]]];
  (*right*) 
  If[index == 2, Return[Rectangle[{3/2 - a/c - b/c, 1/2 - b/c}, {1/2, a/c - 1/2}]]];
  (*down*) 
  If[index == 3, Return[Rectangle[{1/2 - b/c, -1/2}, {a/c - 1/2, (a/c + b/c) - 3/2}]]];
  (*left*)
  If[index == 4, Return[Rectangle[{-1/2, 1/2 - a/c}, {a/c + b/c - 3/2, b/c - 1/2}]]];
]; (* end of GetNonSurjectiveRegion *)


(*

Procedure: FindGroupMembersInRegions

Summary:
  This procedure returns coordinates of elements of a group induces \
  on remainder range. These coordinates are integers. See more in a paper cited in :Summary: section
  of this package.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  F -- a non-surjective or non-injective frame obtained by GetNonInjectiveRegion or
  GetNonSurjectiveRegion.

Output:
  A set of integer points.

*)
FindGroupMembersInRegions[p_, q_, F_] := Module[{c, s},
  TestPythagoreanTripleGenerators[p, q];
   c = p^2 + q^2;
   s = Solve[({p/c, q/c}*u + {-q/c, p/c}*v) \[Element] F, {u, v}, Integers]; 
   Return[Transpose[{s[[All, 1, 2]], s[[All, 2, 2]]}]];
]; (* end of FindGroupMembersInRegions *)


(*

Procedure: BezoutCoefficients

Summary:
  This procedure returns Bézout coefficients used to generate lattices of preimages of elements of
  cyclic group induced on the remainder range. See more in a paper cited in :Summary: section of
  this package.


Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.

Output:
  Two pairs of integers such that p^2 u + q^2 v = 1 (first) and a x + b y = 1 (second).

*)
BezoutCoefficients[p_, q_] := Module[{a, b, c, first, second},
  TestPythagoreanTripleGenerators[p, q];
  b = 2*p*q; a = p^2 - q^2; c = p^2 + q^2; n = c;
  first = ExtendedGCD[p^2, q^2][[2]];
  second = ExtendedGCD[a, b][[2]];
  Return[{first, second}];
]; (* end of BezoutCoefficients *)


(*

Procedure: LatticePoint

Summary:
  This procedure calculate a point of the lattice of preimegaes for an element of a cyclic group
  induced on the remainder range. See more in a paper cited in :Summary: section of this package.

Parameters:
  element -- corrdinates of an element in a cyclic group
  x -- index in (a, -b) direction
  y -- index in c (coeffs[[2]]) direction
  coeffs -- Bézout coefficients (see procedure BezoutCoefficients)
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.

Output:
  A point such that its image in the remainder range is a given element of a cyclic group.

*)
LatticePoint[element_, x_, y_, coeffs_, p_, q_] := Module[{a, b, c},
   TestPythagoreanTripleGenerators[p, q];
   b = 2 p*q; a = p^2 - q^2; c = p^2 + q^2;
   If[ p^2 coeffs[[1]][[1]] + q^2 coeffs[[1]][[2]] != 1, Throw["The first Bézout coefficients are
     inncorect!"] ];
   If[ a coeffs[[2]][[1]] + b coeffs[[2]][[2]] != 1, Throw["The second Bézout coefficients are
     inncorect!"] ];
   Return[( element * p * ( coeffs[[1]][[1]] - coeffs[[1]][[2]] ) / 2 + x * {a, -b} + c * y *
   coeffs[[2]])];
]; (* end of LatticePoint *)


(*

Procedure: IntersectionSetLatticesNonInjective

Summary:
  This procedure returns integer points in a given set such that they have non-injective mapping in
  the target space.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  set -- a rectangular region.

Output:
  A set of integer points of the set such that their mapping under 2D digitized rigid motion is
  non-injective.

*)
IntersectionSetLatticesNonInjective[p_, q_, set_] := Module[{a, b, c, FI, MI, nonInj, s, k, j, coeffs},
   TestPythagoreanTripleGenerators[p, q];
   b = 2 * p * q; a = p^2 - q^2; c = p^2 + q^2;
   coeffs = BezoutCoefficients[p, q];
   FI = { GetNonInjectiveRegion[1, p, q], GetNonInjectiveRegion[2, p, q],
          GetNonInjectiveRegion[3, p, q], GetNonInjectiveRegion[4, p, q] };
   MI = {};
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[1]]]];
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[2]]]];
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[3]]]];
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[4]]]];
   nonInj = {};
   Do[
      s = Solve[( k * p * ( coeffs[[1]][[1]] - coeffs[[1]][[2]] ) / 2 +  x * {a, -b} + c * y *
      coeffs[[2]]) \[Element] set, {x, y}, Integers];
      s = Transpose[{s[[All, 1, 2]], s[[All, 2, 2]]}];
      If[s == {}, Continue[]];
      AppendTo[nonInj, {}];
      Do[nonInj[[-1]] = Join[nonInj[[-1]], {LatticePoint[k, j[[1]], j[[2]], coeffs, p, q]}], {j, s}];
  , {k, MI}]; (* end of outer Do[] *)
  Return[nonInj];
]; (* end of IntersectionSetLatticesNonInjective *)


(*

Procedure: IntersectionSetLatticesNonSurjective

Summary:
  This procedure returns integer points in a given set such that in a neighbourhood of their images
  in the target space exist points without preimage.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  set -- a rectangular region.

Output:
  A set of integer points of the set such that in a nieghbourhood of their images in the target
  space exist non-surjective points.

*)
IntersectionSetLatticesNonSurjective[p_, q_, set_] := Module[{a, b, c, FI, MI, nonInj, s, k, j},
   TestPythagoreanTripleGenerators[p, q];
   b = 2 * p * q; a = p^2 - q^2; c = p^2 + q^2;
   coeffs = BezoutCoefficients[p, q];
   FI = { GetNonSurjectiveRegion[1, p, q], GetNonSurjectiveRegion[2, p, q],
          GetNonSurjectiveRegion[3, p, q], GetNonSurjectiveRegion[4, p, q] };
   MI = {};
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[1]]]];
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[2]]]];
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[3]]]];
   MI = Join[MI, FindGroupMembersInRegions[p, q, FI[[4]]]];
   nonInj = {};
   Do[
      s = Solve[( k * p * ( coeffs[[1]][[1]] - coeffs[[1]][[2]] ) / 2 +  x * {a, -b} + c * y *
      coeffs[[2]]) \[Element] set, {x, y}, Integers];
      s = Transpose[{s[[All, 1, 2]], s[[All, 2, 2]]}];
      If[s == {}, Continue[]];
      AppendTo[nonInj, {}];
      Do[nonInj[[-1]] = Join[nonInj[[-1]], {LatticePoint[k, j[[1]], j[[2]], coeffs, p, q]}], {j, s}];
  , {k, MI}]; (* end of outer Do[] *)
  Return[nonInj];
]; (* end of IntersectionSetLatticesNonSurjective *)

End[]
EndPackage[]

