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

  // Caltech copyright info
  text += "Copyright 2003 California Institute of Technology.\n"; 
  text += "\n";
  text += "To contact the original authors, write to: \n";
  text += "drkent@users.sourceforge.net\n";
  text += "or\n";
  text += "mtfeldmann@users.sourceforge.net   \n";
  text += "\n";
  text += "\n";
  text += "This program is free software; you can redistribute it and/or modify\n";
  text += "it under the terms of the GNU General Public License as published by\n";
  text += "the Free Software Foundation and so long as the above copyright\n";
  text += "notice, this paragraph and the following three paragraphs appear in\n";
  text += "all copies.  \n";
  text += "\n";
  text += "This program is distributed in the hope that it will be useful, but\n";
  text += "WITHOUT ANY WARRANTY; without even the implied warranty of\n";
  text += "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. In no event \n";
  text += "shall California Institute of Technology or the authors be liable \n";
  text += "to any party for direct, indirect, special, incidental or\n";
  text += "consequential damages, including lost profits, arising out of the use\n";
  text += "of this software and its documentation, even if the California\n";
  text += "Institute of Technology or the authors have been advised of the\n";
  text += "possibility of such damage. Lastly, the California Institute of\n";
  text += "Technology and the authors have no obligations to provide maintenance,\n";
  text += "support, updates, enhancements or modifications.   \n";
  text += "\n";
  text += "To receive a copy of the GNU General Public License, go to:\n";
  text += "\n";
  text += "http://www.gnu.org/licenses/gpl.txt\n";
  text += "or write to:\n";
  text += "The Free Software Foundation, Inc.\n";
  text += "59 Temple Place, Suite 330\n";
  text += "Boston, MA 02111--1307 USA\n";
  text += "\n";
  text += "--------------------------------------------\n";
  text += "\n";

  strm << text;
  return strm;
}
