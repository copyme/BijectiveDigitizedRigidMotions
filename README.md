Description
===========
BijectiveDigitzedRigidMotions is a set of Mathematica packages for which the main objective is to
implement algorithms used in the study of bjective digitized rigid motions.  These algorithms were
described/introduced in: [Kacper Pluta , Pascal Romon, Yukiko Kenmochi, Nicolas Passat,Bijective
Rigid Motions of the 2D Cartesian Grid, DGCI2016, Springer
2016](http://link.springer.com/chapter/10.1007/978-3-319-32360-2_28) ([free
preprint](https://hal.archives-ouvertes.fr/hal-01275598v2)); [Kacper Pluta, Pascal Romon, Yukiko
Kenmochi, Nicolas Passat. Bijectivity certification of 3D digitized rotations, CTIC2016, Springer
2016](http://link.springer.com/chapter/10.1007%2F978-3-319-39441-1_4) ([free
print](https://hal.archives-ouvertes.fr/hal-01315226v1)).

Quick Install
=============

1. Download or clone the repository.
2. [Install the packages](http://support.wolfram.com/kb/5648)

Examples
================

To check if a 2D digitized rigid motion is bijective while restricted to a finite digital set S
while using ForwardAlgorithm:

```
Needs["ForwardAlgorithm`"];
CheckInjectivity[3, 4, {1/3, 1/2}, S]
```

To check if a 2D digitized rigid motion is bijective while restricted to a finite digital set S while
using BackwardAlgorithm:

```
Needs["BackwardAlgorithm`"];
IntersectionSetLatticesNonInjective[3, 4, {1/3, 1/2}, S]
```

To check if a 3D digitized rotation given by a Lipschitz quaternion is bijective:

```
Needs["QuaternionCertification`"];
IntersectionSetLatticesNonInjectiveCertifyQuaternion[{3,0,0,1}]
```

Additional information
================

The algorithms were implemented in Mathematica 10, therefore, Mathematica 10 or higher is required.

