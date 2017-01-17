(* :Title: ForwardAlgorithm, Version 1.0, 2016*)

(* :Author: Kacper Pluta, kacper.pluta@esiee.fr
	    Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France 
*)

(* :Summary:

  A package which implements an algorithm which allows to verify if a digitized rigid motions
  constrained to a finite set S is bijective.
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

(* :Context: ForwardAlgorithm` *)

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

  2D digitized rigid motions, non-injective
*)


(* :Mathematica Version: 10.0.0 *)


BeginPackage["ForwardAlgorithm`"]


CheckInjectivity::usage = "CheckInjectivity[p, q, t, set] returns a subset of S";

Begin["`Private`"]

(*

Procedure: TestPythagoreanTripleGenerators

Summary:
  This procedure tests if generators of primitive Pythagorean triple are correct. If not it throws
  exceptions.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1

*)
TestPythagoreanTripleGenerators[p_, q_] := (
  If[GCD[p, q] != 1 || EvenQ[p - q], Throw["Pythagorean rotation generates are incorrect!"]];
)


(*

Procedure: TestTranslationVector

 Summary:
  This procedure tests if translation vector is valid.

Parameters:
  t -- 2D vector with rational elements

*)
TestTranslationVector[t_] := (
   If[Length[t] != 2, Throw["Translation vector has wrong dimension!"]]; 
   If[NotElement[Head[t[[1]]], {Rational, Integer}] || NotElement[Head[t[[2]]], {Rational,
   Integer}], Throw["Translation vector has irrational elements!"]];
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

Procedure: RemainderMap

Summary:
  This procedure implements the remainder map.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  t -- 2D translation vector
  x -- integer point

Output:
  A non-injective region of a given index.

*)
RemainderMap[p_, q_, t_, x_] := Module[{a, b, c, U},
  TestPythagoreanTripleGenerators[p, q];
  TestTranslationVector[t];
  TestTranslationVector[x];
  (*primitive Pythagorean triple*)
  b = 2 p q; a = p^2 - q^2; c = p^2 + q^2;
  U = {{a/c, -b/c},{b/c, a/c}}.x + t;
  Return[U - Round[U]];
];


(*

Procedure: CheckInjectivity

Summary:
  This procedure returns a subset of input set such that each point of such subset induce
  non-injectivty.

Parameters:
  p -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  q -- generator of primitive Pythagorean triple such that p - q is odd and gcd of p and q is 1.
  t -- 2D translation vector
  set -- a finite set of integer points

Output:
  A non-injective region of a given index.

*)
CheckInjectivity[p_, q_, t_, set_] := Module[{xp, B, F1, F2},
  B = {};
  F1 = GetNonInjectiveRegion[1, p, q];
  F2 = GetNonInjectiveRegion[2, p, q];
  
  Do[
    xp = RemainderMap[p, q, t, x];
    If[xp \[Element] F1 && MemberQ[set, x + {0,1}], AppendTo[B, x]; AppendTo[B, x + {0,1}]];
    If[xp \[Element] F2 && MemberQ[set, x + {1,0}], AppendTo[B, x]; AppendTo[B, x + {1,0}]];
  ,{x, set}];
  Return[B];
]; (* end of CheckInjectivity *)

End[]
EndPackage[]

