/*
  Copyright (c) Amos G. Anderson 2007
  Distributed under GNU general public license (GPL)
  No guarantee or warantee regarding usability or stability is expressed or implied.
  nitroamos@gmail.com
*/

#ifndef QMCSURFER_H
#define QMCSURFER_H

#include <string>

#include "Array2D.h"
#include "QMCWalkerData.h"
#include "QMCFunctions.h"
#include <iostream>
#include <iomanip>

using namespace std;

/**
   This set of functions will "surf" a wavefunction. At the moment, I can't be
   bothered to document this more carefully. Basically, you have a set of positions
   including the locations of the nuclei, and the locations of the electrons.
   You can start or stop at any of those positions, or go in a circle around one of
   the positions.

   When asked to input a value, the bracketed value is the one that will accepted
   as default if you only hit return.

   It will print out some energy terms, psi. Using GCC/Linux, the output is color coded;
   it will change when you cross a node.

   Just use a regular input file, except add a "use_surfer" flag. Obviously, you want to
   run this in interactive mode, only one processor, etc.
*/

class QMCSurfer
{
 private:
  QMCFunctions * QMF;
  Array2D<double> R;
  QMCWalkerData walkerData;

 public:
  QMCSurfer();
  ~QMCSurfer();

  int mainMenu(QMCFunctions * QMF, int iteration,
	       Array2D<double> newR);

  void interparticleDistanceMatrix();
  void equipotentialSurface();
  void scanEnergies(int, int, int, int, double, double, double);
  void surfaceExplorer();
  void grid3D();
};

#endif
