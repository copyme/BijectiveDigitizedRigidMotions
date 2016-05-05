(* :Title: QuaternionCertification, Version 1.0, 2016*)

(* :Author: Kacper Pluta, kacper.pluta@esiee.fr
	    Laboratoire d'Informatique Gaspard-Monge - LIGM, A3SI, France 
*)

(* :Summary:

  A package which implements algorithm which allows to verify if a given Lipschitz quaternion induce
  bijective digitized rotations. The algorithm was introduced in:

@inproceedings{Pluta:CTIC:2016,
	author = {Pluta, K. and Romon, P. and Kenmochi, Y. and Passat, N.},
	booktitle = {{CTIC}},
	publisher = {Springer},
	series = {{Lecture Notes in Computer Science}},
	title = {{Bijectivity certification of 3D digitized rotations}},
	year = {Forthcoming 2016}
}

*)

(* :Context: QuaternionCertification` *)

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

  3D digitized rotations, certification, quaternions, Lipschitz quaternions, bijectivity
*)


(* :Mathematica Version: 10.0.0 *)

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

Procedure: QuaternionRgionMultiplication

Summary:
  This procedure builds a system [A|B] (see paper cited in :Summary: section of this package).

Parameters:
  a, b, c, d -- a quaternion represented by four integers.
  bb1, bb2, bb3, bb4 -- a four dimensional integer point.

Output:
  The augment matrix [A|B].

*)
SystemAB[a_, b_, c_, d_, bb1_, bb2_, bb3_, bb4_] := 
  Return[ {
   { b, c, d, -b, -c, -d, bb1 }, { a, -d, c, -a, -d, c, bb2 },
   { d, a, -b, d, -a, -b, bb3 }, { -c, b, a, -c, b, -a, bb4 } 
  } ];
(* end of SystemAB *)


(*

Procedure: CertifyQuaternion

Summary:
  Check if a given quaternion generates 3D bijective digitized rotations.

Parameters:
  q -- a Lipschitz quaternion.

Output:
  True if q leads to 3D bijective digitized rotation and False otherwise.

*)
CertifyQuaternion[q_] := Module[ {p, pp, cube1, cube2, bounds, y, h, i, j, k, pivot47},
  If[MemberQ[Map[IntegerQ, q], False], Throw["This is not a Lipschitz quaternion."]];
  p = q / GCD[ q[[1]], q[[2]], q[[3]], q[[4]] ];
  pp = { p[[1]], -p[[2]], -p[[3]], -p[[4]] };
  cube1 = TransformedRegion[OpenUnitCube[], QuaternionRegionMultiplication[ pp, { Indexed[#, 1],
  Indexed[#, 2], Indexed[#, 3], Indexed[#, 4] } ] &];
  cube2 = TransformedRegion[OpenUnitCube[], QuaternionRegionMultiplication[ { Indexed[#, 1],
  Indexed[#, 2], Indexed[#, 3], Indexed[#, 4] }, pp] & ]; 
  bounds = Round[RegionBounds[cube1]];  
  Do[
     y = {h, i, j, k};
     (*If zero the y times q goes to R^3*)
     If[ h * q[[1]] - i * q[[2]] - j * q[[3]] - k * q[[4]] != 0, Continue[]];
     pivot47 = HermiteDecomposition[ SystemAB[q[[1]], q[[2]], q[[3]], q[[4]], y[[1]], y[[2]], y[[3]],
                                              y[[4]]]][[2]][[4]][[7]]; 
     (*System has solution iff the last pivot is zero.*)
     If[y \[NotElement] cube1 || pivot47 != 0, Continue[]];
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

