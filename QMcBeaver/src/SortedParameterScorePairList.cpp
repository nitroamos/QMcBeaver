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

#include "SortedParameterScorePairList.h"

SortedParameterScorePairList::SortedParameterScorePairList()
{
}

SortedParameterScorePairList::SortedParameterScorePairList(SortedParameterScorePairList & SPSL)
{
  PSPList = SPSL.PSPList;
}

int SortedParameterScorePairList::size()
{
  return PSPList.size();
}

void SortedParameterScorePairList::add(const ParameterScorePair & PSP)
{
  PSPList.push_back(PSP);
  PSPList.sort();
}

ParameterScorePair SortedParameterScorePairList::get(int i)
{
  list <ParameterScorePair>::iterator iter;

  int index = 0;
  for(iter=PSPList.begin();iter!=PSPList.end();++iter)
    {
      if( index == i )
	{
	  return *iter;
	}
      index++;
    }

  cerr << "ERROR: accessing element beyond the range of "
       << "SortedParameterScorePairList" << endl;
  exit(0);
  return *iter;
}

void SortedParameterScorePairList::clear()
{
  PSPList.clear();
}

void SortedParameterScorePairList::operator=(const 
					     SortedParameterScorePairList & 
					     SPSL)
{
  PSPList = SPSL.PSPList;
}





