(* :Title: ExtendedForwardAlgorithm, Version 1.1, 2017*)

(* :Author: Kacper Pluta, kacper.pluta@esiee.fr
	    Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France 
*)

(* :Summary:

  A package which implements algorithms which allow to verify if a digitized rigid motions
  constrained to a finite set S is bijective.
  The algorithm was introduced in:

@Article{
  Pluta2017,
  author="Pluta, Kacper and Romon, Pascal and Kenmochi, Yukiko and Passat, Nicolas",
  title="Bijective Digitized Rigid Motions on Subsets of the Plane",
  journal="Journal of Mathematical Imaging and Vision",
  year="2017",
  pages="1--22",
  issn="1573-7683",
  doi="10.1007/s10851-017-0706-8",
  url="http://dx.doi.org/10.1007/s10851-017-0706-8"
}

*)

(* :Context: ExtendedForwardAlgorithm` *)

(* :Package Version: 1.1 *)

(* : Copyright: Kacper Pluta, 2017
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
  An implementation of ExtendedForwardAlgorithm has been added, 2017
  Created by Kacper Pluta at University Paris-Est, 2016
*)

(* :Keywords: 

  2D digitized rigid motions, non-injective
*)


(* :Mathematica Version: 10.0.0 *)


BeginPackage["ExtendedForwardAlgorithm`"]


CheckInjectivity::usage = "CheckInjectivity[p, q, t, set] returns a subset of S";
CheckInjectivityRange::usage = "CheckInjectivityRange[p, q, t, set] returns a set of hinge angles";

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
  If[GCD[p, q] != 1 || EvenQ[p - q], Throw["The Pythagorean rotation generates are incorrect!"]];
)


(*

Procedure: TestTranslationVector

 Summary:
  This procedure tests if a translation vector is valid.

Parameters:
  t -- 2D vector with rational elements.

*)
TestTranslationVector[t_] := (
   If[Length[t] != 2, Throw["The translation vector has wrong dimension!"]]; 
   If[NotElement[Head[t[[1]]], {Rational, Integer}] || NotElement[Head[t[[2]]], {Rational,
   Integer}], Throw["The translation vector has irrational elements!"]];
)


(*

Procedure: TestHingeAngle

 Summary:
  This procedure tests if a hinge angle is valid.

Parameters:
  h -- 4D vector with integer elements. The format is {p_1, p_2, k, s}.

*)
TestHingeAngle[h_] := (
  If[Length[h] != 4, Throw["The hinge angle is given in inccrect form: wrong length!"]];
  If[DeleteDuplicates[Head /@ h] != {Integer}, 
    Throw["The hinge angle is given in inccrect form: non-integer elements!"]
  ];
)


(*

Procedure: CompareHingePythagorean

 Summary:
  This procedure compares hinge and Pythagorean angles.

Parameters:
  h -- 4D vector with integer elements. The format is {p_1, p_2, k, s}.
  p -- 3D vector with integer elements. The format is {b, a, c}
  t -- 2D vector with rational elements.

Comments:
  It can throw an exception if some parameters are not valid!

Output:
  True if the Pythagorean angle is bigger than the hinge angle and false otherwise. 
*)
CompareHingePythagorean[h_, p_, t_] := Module[{u, s, l},
   TestHingeAngle[h]; TestTranslationVector[t];
   If[h[[4]] == 1,
    s = p[[3]] h[[2]] ((2 h[[3]] + 1) Denominator[t[[2]]] - 2 Numerator[t[[2]]]) - 2 p[[1]] Denominator[t[[2]]] (h[[1]]^2 + h[[2]]^2);
    l = Sqrt[h[[1]]^2 + h[[2]]^2 - (h[[3]] - t[[2]] + 1/2)^2];
    u = -2 p[[3]] l h[[1]] Denominator[t[[2]]];,
    s = p[[3]] h[[1]] ((2 h[[3]] + 1) Denominator[t[[1]]] - 2 Numerator[t[[1]]]) - 2 p[[1]] Denominator[t[[1]]] (h[[1]]^2 + h[[2]]^2);
    l = Sqrt[h[[1]]^2 + h[[2]]^2 - (h[[3]] - t[[1]] + 1/2)^2];
    u = -2 p[[3]] l h[[1]] Denominator[t[[1]]];
   ];
   If[Im[l] != 0, Throw["On the hinge angle is not valid: a half-grid index is too high!"];];
   If[NonNegative[s], Return[True], Return[TrueQ[s^2 < u^2]];];
];


(*

Procedure: CompareHingeHinge

 Summary:
  This procedure compares two hinge angles.

Parameters:
  h -- 4D vector with integer elements. The format is {p_1, p_2, k, s}.
  g -- 4D vector with integer elements. The format is {p_1, p_2, k, s}.
  t -- 2D vector with rational elements.

Comments:
  It can throw an exception if some parameters are not valid!

Output:
  True if the second angle is bigger than the first one and false otherwise. 
*)
CompareHingeHinge[h_, g_, t_] := Module[{A, B, C, D, l, ll, den},
   TestHingeAngle[h]; TestHingeAngle[g]; TestTranslationVector[t];
   If[h == g, Throw["The angles are the same!"];];
   (*Set the variales with respect to the critical lines*)
   
   If[h[[4]] == 1 && g[[4]] == 1,
    A = (g[[1]]^2 + g[[2]]^2) h[[2]] (2 h[[3]] - t[[2]] + 1);
    B = (h[[1]]^2 + h[[2]]^2) g[[2]] (2 g[[3]] - t[[2]] + 1);
    C = (h[[1]]^2 + h[[2]]^2) g[[1]];
    D = (g[[1]]^2 + g[[2]]^2) h[[1]];
    l = Sqrt[h[[1]]^2 + h[[2]]^2 - (h[[3]] - t[[2]] + 1/2)^2];
    ll = Sqrt[g[[1]]^2 + g[[2]]^2 - (g[[3]] - t[[2]] + 1/2)^2];
    den = Denominator[t[[2]]];
   ];
   If[h[[4]] == 0 && g[[4]] == 1,
    A = (g[[1]]^2 + g[[2]]^2) h[[1]] (2 h[[3]] - t[[1]] + 1);
    B = (h[[1]]^2 + h[[2]]^2) g[[2]] (2 g[[3]] - t[[2]] + 1);
    C = (h[[1]]^2 + h[[2]]^2) g[[1]];
    D = (g[[1]]^2 + g[[2]]^2) h[[2]];
    l = Sqrt[h[[1]]^2 + h[[2]]^2 - (h[[3]] - t[[1]] + 1/2)^2];
    ll = Sqrt[g[[1]]^2 + g[[2]]^2 - (g[[3]] - t[[2]] + 1/2)^2];
    den = Denominator[t[[2]]] * Denominator[t[[1]]];
   ];
   If[h[[4]] == 1 && g[[4]] == 0,
    A = (g[[1]]^2 + g[[2]]^2) h[[2]] (2 h[[3]] - t[[2]] + 1);
    B = (h[[1]]^2 + h[[2]]^2) g[[1]] (2 g[[3]] - t[[1]] + 1);
    C = (h[[1]]^2 + h[[2]]^2) g[[2]];
    D = (g[[1]]^2 + g[[2]]^2) h[[1]];
    l = Sqrt[h[[1]]^2 + h[[2]]^2 - (h[[3]] - t[[2]] + 1/2)^2];
    ll = Sqrt[g[[1]]^2 + g[[2]]^2 - (g[[3]] - t[[1]] + 1/2)^2];
    den = Denominator[t[[2]]] * Denominator[t[[1]]];
   ];
   If[h[[4]] == 0 && g[[4]] == 0,
    A = (g[[1]]^2 + g[[2]]^2) h[[1]] (2 h[[3]] - t[[1]] + 1);
    B = (h[[1]]^2 + h[[2]]^2) g[[1]] (2 g[[3]] - t[[1]] + 1);
    C = (h[[1]]^2 + h[[2]]^2) g[[2]];
    D = (g[[1]]^2 + g[[2]]^2) h[[2]];
    l = Sqrt[h[[1]]^2 + h[[2]]^2 - (h[[3]] - t[[1]] + 1/2)^2];
    ll = Sqrt[g[[1]]^2 + g[[2]]^2 - (g[[3]] - t[[1]] + 1/2)^2];
    den = Denominator[t[[1]]];
   ];
   If[Im[l] != 0 || Im[ll] != 0, Throw["On of the angles is not valid: a half-grid index is too high!"];];
   If[Positive[den (A - B)] && Negative[2 den (C ll - D l)], Return[True];];
   If[Negative[den (A - B)] && Positive[2 den (C ll - D l)], Return[False];];
   If[Positive[den (A - B)] && Positive[2 den (C ll - D l)],
   If[Positive[den^2 (A - B)^2 - 4 den^2 (D^2 l^2 + C^2 ll^2)], 
     Return[True],
     Return[TrueQ[(den^2 (A - B)^2 - 4 den^2 (D^2 l^2 + C^2 ll^2))^2 < 64 den^4 C^2 D^2 l^2 ll^2]];];
   ];
   If[Negative[den (A - B)] && Negative[2 den (C ll - D l)],
   If[Positive[4 den^2 (C ll - D l)^2 - den^2 (A^2 + B^2)], 
     Return[True],
     Return[TrueQ[4 den^4 A^2 B^2 > den^4 (4 (C ll - D l)^2 - (A^2 + B^2))^2]];];
   ];
]; (* end of CompareHingeHinge *)


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



(*

Procedure: Get4Neighborhod

 Summary:
  This procedure compares hinge and Pythagorean angles.

Parameters:
  set -- a finite set of integer points
  x -- a finite set of integer points

Comments:
  It can throw an exception if some parameters are not valid!

Output:
  True if the Pythagorean angle is bigger than the hinge angle and false otherwise. 
*)
Get4Neighborhod[set_, x_] := Module[{N},
 If[NotElement[x, set], Throw["The point is not in the set!"];];
 N = {x};
 If[Element[x + {1, 0}, set], AppendTo[set, x]; ];
 If[Element[x + {-1, 0}, set], AppendTo[set, x]; ];
 If[Element[x + {0, 1}, set], AppendTo[set, x]; ];
 If[Element[x + {0, -1}, set], AppendTo[set, x]; ];
 Return[N];
]; 



FindInjectiveRanges[p_, q_, t_, set_] := Module[{C},
 Do[
  If[Length[CheckInjectivity[p_, q_, t_, Get4Neighborhod[set, x]]] != 0, Return[{}]]; 
  
 , {x, set}
 ];
];


End[]
EndPackage[]

