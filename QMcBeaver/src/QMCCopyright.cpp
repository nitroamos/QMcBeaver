//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2000-2.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "QMCCopyright.h"

ostream & operator<<(ostream & strm, QMCCopyright & rhs)
{
  string text = "";

  // ascii beaver

  text += "             ___\n";
  text += "           .=\"   \"=._.---.\n";
  text += "         .\"         c \' Y\'`p\n";
  text += "        /   ,       `.  w_/\n";
  text += "        |   \'-.   /     /\n";
  text += "  _,..._|      )_-\\ \\_=.\\\n";
  text += " `-....-\'`------)))`=-\'\"`\'\"\n";


  // text

  text += "\n";
  text += "                   QMcBeaver\n";
  text += "\n";
  text += "               Version: " +
    StringManipulation::longToString(VERSION) + "\n";
  text += "\n";
  text += "\n";
  text += "                Constructed by\n";
  text += "\n";
  text += "        David Randall \"Chip\" Kent IV\n";
  text += "                    and\n";
  text += "            Michael Todd Feldmann\n";
  text += "\n";
  text += " Copyright 2000-2003.  All rights reserved.\n";
  text += "\n";
  text += "drkent@users.sourceforge.net feldmann@cacr.caltech.edu\n";
  text += "--------------------------------------------\n";
  text += "\n";

  strm << text;
  return strm;
}
