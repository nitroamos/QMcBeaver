//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#ifndef SORTEDPARAMETERSCORPAIRLIST_H
#define SORTEDPARAMETERSCORPAIRLIST_H

#include <list>

#include "ParameterScorePair.h"

using namespace std;

/**
  A sorted list of ParameterScorePair objects where the objects are ordered
  in an increasing order.
  */

class SortedParameterScorePairList
{
public:
  /**
    Creates an empty instance of this class.
    */
  SortedParameterScorePairList();

  /**
    Createsan instance of this class which is equal to another instance.

    @param SPSL this object to which this one will be made equal.
    */
  SortedParameterScorePairList(SortedParameterScorePairList & SPSL);

  /**
    Gets the number of elements in this list.
   
    @return number of elements in this list.
    */
  int size();

  /**
    Adds a new ParameterScorePair to this list.

    @param PSP new element to add to this list.
    */
  void add(const ParameterScorePair & PSP);

  /**
    Gets the ith element.
    
    @param i index of the element to return.
    @return the ith element of the list.
    */
  ParameterScorePair get(int i);

  /**
    Remove all elements from this list.
    */
  void clear();

  /**
    Sets two objects equal to one another.

    @param SPSL object to set this object equal to.
    */
  void operator=(const SortedParameterScorePairList & SPSL);

private:
  list <ParameterScorePair> PSPList;
};

#endif
