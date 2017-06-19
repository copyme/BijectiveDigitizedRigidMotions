(* :Title: QuaternionCertification, Version 1.1, 2017*)

(* :Author: Kacper Pluta, kacper.pluta@esiee.fr
	    Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France 
*)

(* :Summary:

  A package which implements algorithm which allows to verify if a given Lipschitz quaternion induce
  bijective digitized rotations. The algorithm was introduced in:

@Inbook{Pluta:CTIC:2016,
  author="Pluta, Kacper and Romon, Pascal and Kenmochi, Yukiko and Passat, Nicolas",
  editor="Bac, Alexandra and Mari, Jean-Luc",
  title="Bijectivity Certification of 3D Digitized Rotations",
  bookTitle="Computational Topology in Image Context: 6th International Workshop, CTIC 2016,
             Marseille, France, June 15-17, 2016, Proceedings", 
  year="2016",
  publisher="Springer International Publishing",
  address="Cham",
  pages="30--41",
  isbn="978-3-319-39441-1",
  doi="10.1007/978-3-319-39441-1_4",
  url="http://dx.doi.org/10.1007/978-3-319-39441-1_4"
}

*)

(* :Context: QuaternionCertification` *)

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
  Created by Kacper Pluta at University Paris-Est, 2016
*)

(* :Keywords: 

  3D digitized rotations, certification, quaternions, Lipschitz quaternions, bijectivity
*)


(* :Mathematica Version: 11.0.0 *)

BeginPackage["QuaternionCertification`"]

 CertifyQuaternion::usage = "CertifyQuaternion[q] returns True or False.";

Begin["`Private`"]

(*

Procedure: OpenUnitCube

Summary:
  Open cube which will work with quaternion multiplication

*)
OpenUnitCube[] := ImplicitRegion[ w == 0 && -1/2 < x < 1/2 && -1/2 < y < 1/2 && -1/2 < z < 1/2, {w,
x, y, z}];


(*

Procedure: QuaternionRgionMultiplication

Summary:
  This procedure allows to apply Hamilton product to regions.

Parameters:
  v -- a quaternion represented as a set of cardinality 4.
  c -- a quaternion represented as a set of cardinality 4.

Output:
  A non-surjective region of a given index.

*)
QuaternionRegionMultiplication[v_, c_] := (
   Return[ { 
     v[[1]] * c[[1]] - v[[2]] * c[[2]] - v[[3]] * c[[3]] - v[[4]] * c[[4]],
     v[[1]] * c[[2]] + v[[2]] * c[[1]] + v[[3]] * c[[4]] - v[[4]] * c[[3]],
     v[[1]] * c[[3]] - v[[2]] * c[[4]] + v[[3]] * c[[1]] + v[[4]] * c[[2]],
     v[[1]] * c[[4]] + v[[2]] * c[[3]] - v[[3]] * c[[2]] + v[[4]] * c[[1]] 
  } ];
); (* end of QuaternionRgionMultiplication *)


(*

Procedure: SystemA

Summary:
  This procedure builds the matrix A from the system (see paper cited in :Summary: section of this
  package).

Parameters:
  a, b, c, d -- a quaternion represented by four integers.

Output:
  The matrix A.

*)
SystemA[a_, b_, c_, d_] := 
  Return[ {
   { b, c, d, -b, -c, -d }, { a, -d, c, -a, -d, c },
   { d, a, -b, d, -a, -b }, { -c, b, a, -c, b, -a } 
  } ];
(* end of SystemA *)


(*

Procedure: CheckMultiplicity

Summary:
  This procedure checks of elements of the second list are multiples of the respective elements of
  the first list

Parameters:
  Two lists of integers which are of the same length.

Output:
  True if all the elements of the second lists are integer multiples of the elements of the first
  list.

*)
CheckMultiplicity[ListA_, ListB_] := Module[{result},
  If[Length[ListA] != Length[ListB], Return[False],
    result = True;
    Do[ result = result && IntegerQ[ ListB[[x]] / ListA[[x]] ], {x, 1, Length[ListA]} ];
    Return[result];
  ];
];
(* end of CheckMultiplicity *)


(*

Procedure: HasIntegerSolutions

Summary:
  This procedure checks if y is a solution to the system Ax = y (see the paper cited in :Summary:)
  To check the feasibility of the system Ax = y we use a similar condition to the one give by
  A.  Schrijver in Theory of linear and integer programming (page 51)

Parameters:
  A Smith Normal Form of the matrix A and a four dimensional integer point y

Output:
  True if Ax = y has integer solutions.

*)
HasIntegerSolutions[SNF_, y_] := Module[{yp, diag},
  diag = Diagonal[ SNF[[2]] ];
  yp = SNF[[1]].y;
  If[ yp[[-1]] != 0, Return[False], Return[ CheckMultiplicity[ diag[[;;-2]], yp[[;;-2]] ] ] ];
];
(* end of HasIntegerSolutions *)


(*

Procedure: CertifyQuaternion

Summary:
  Check if a given quaternion generates 3D bijective digitized rotations.

Parameters:
  q -- a Lipschitz quaternion.

Output:
  True if q leads to 3D bijective digitized rotation and False otherwise.

*)
CertifyQuaternion[q_] := Module[ {p, pp, cube1, cube2, bounds, y, h, i, j, k, SNF},
  If[MemberQ[Map[IntegerQ, q], False], Throw["This is not a Lipschitz quaternion."]];
  p = q / GCD[ q[[1]], q[[2]], q[[3]], q[[4]] ];
  pp = { p[[1]], -p[[2]], -p[[3]], -p[[4]] };
  cube1 = TransformedRegion[OpenUnitCube[], QuaternionRegionMultiplication[ pp, { Indexed[#, 1],
  Indexed[#, 2], Indexed[#, 3], Indexed[#, 4] } ] &];
  cube2 = TransformedRegion[OpenUnitCube[], QuaternionRegionMultiplication[ { Indexed[#, 1],
  Indexed[#, 2], Indexed[#, 3], Indexed[#, 4] }, pp] & ]; 
  bounds = Round[RegionBounds[cube1]];
  SNF = SmithDecomposition[ SystemA[ p[[1]], p[[2]], p[[3]], p[[4]] ] ];
  (*Compute GCDs for minors of the all valid dimensions*)
  Do[
     y = {h, i, j, k};
     (*If zero the y times q goes to R^3*)
     If[ h * p[[1]] - i * p[[2]] - j * p[[3]] - k * p[[4]] != 0, Continue[]];
     If[y \[NotElement] cube1 || !HasIntegerSolutions[ SNF, y ], Continue[]];
     If[y \[NotElement] cube2, Goto[end]];
   , {h, bounds[[1]][[1]], bounds[[1]][[2]]},
   {i, bounds[[2]][[1]], bounds[[2]][[2]]},
   {j, bounds[[3]][[1]], bounds[[3]][[2]]},
   {k, bounds[[4]][[1]], bounds[[4]][[2]]}
  ];
  (* Calling just Return[] in the If[] will take us at the end where is the second Return and output
    will be True. *) 
  Return[True];
  Label[end];
  Return[False];
]; (* end of CertifyQuaternion *)

End[]
EndPackage[]

