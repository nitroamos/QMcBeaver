//            QMcBeaver
//
//         Constructed by 
//
//     Michael Todd Feldmann 
//              and 
//   David Randall "Chip" Kent IV
//
// Copyright 2002.  All rights reserved.
//
// drkent@users.sourceforge.net mtfeldmann@users.sourceforge.net

#include "XMLElement.h"

XMLElement::XMLElement()
{
  initialize(false, true);
}
    

XMLElement::XMLElement(map<string,string> * entities)
{
  initialize(entities, false, true);
}


XMLElement::XMLElement(bool skipLeadingWhitespace)
{
  initialize(skipLeadingWhitespace, true);
}


XMLElement::XMLElement(map<string,string> * entities,
		       bool skipLeadingWhitespace)
{
  initialize(entities, skipLeadingWhitespace, true);
}


XMLElement::XMLElement(map<string,string> * entities,
		       bool skipLeadingWhitespace,
		       bool fillBasicConversionTable)
{
   initialize(entities, skipLeadingWhitespace, fillBasicConversionTable); 
}


void XMLElement::initialize(bool skipLeadingWhitespace,
			    bool fillBasicConversionTable)
{
   initialize(&entitiesInstance, skipLeadingWhitespace, 
	      fillBasicConversionTable);
}

void XMLElement::initialize(map<string,string> * entities,
			    bool skipLeadingWhitespace,
			    bool fillBasicConversionTable)
{
  this->ignoreWhitespace = skipLeadingWhitespace;
  this->contents = "";
  this->entities = entities;
  this->lineNr = 0;

  if (fillBasicConversionTable) 
    {
      string key = "amp";
      string val = "&";
      (*this->entities)[key] = val; 

      key = "quot";
      val = '"';
      (*this->entities)[key] = val; 

      key = "apos";
      val = '\'';
      (*this->entities)[key] = val; 

      key = "lt";
      val = '<';
      (*this->entities)[key] = val; 

      key = "gt";
      val = '>';
      (*this->entities)[key] = val; 
    }
}

int XMLElement::countChildren()
{
  return (int)this->children.size();
}

void XMLElement::addChild(XMLElement & child)
{
  this->children.push_back(child);
}


void XMLElement::setAttribute(string & name, string & value)
{
  this->attributes[name] = value; 
}


void XMLElement::setAttribute(string & name, int value)
{
  this->attributes[name] = StringManipulation::intToString(value);
}


void XMLElement::setAttribute(string & name, double value)
{
  this->attributes[name] = StringManipulation::doubleToString(value);
}

void XMLElement::parse(string & file)
{
  ifstream reader(file.c_str());
  this->parse(reader);
  reader.close();
}


void XMLElement::parse(istream & reader)
{
  this->parse(reader, /*startingLineNr*/ 1);
}

void XMLElement::parse(istream & reader, int startingLineNr)
{
  this->charReadTooMuch = '\0';
  this->reader = &reader;
  this->parserLineNr = startingLineNr;
  
  for (;;) 
    {
      unsigned char ch = this->scanWhitespace();

      if (ch != '<') 
	{
	  throw this->expectedInput("<");
	}

      ch = this->readChar();
    
      if ((ch == '!') || (ch == '?')) 
	{
	  this->skipSpecialTag(0);
	} 
      else 
	{
	  this->unreadChar(ch);
	  this->scanElement(*this);
	  return;
	}
    }
}

void XMLElement::removeChild(XMLElement & child)
{
  for(list<XMLElement>::iterator it=this->children.begin();
      it != this->children.end(); it++)
    {
      if( *it == child )
	{
	  it = this->children.erase(it);
	}
    }
}


void XMLElement::removeAttribute(string & name)
{
  this->attributes.erase(name);
}

void XMLElement::setContent(string & content)
{
  this->contents = content;
}

string XMLElement::getContent()
{
  return this->contents;
}

string XMLElement::getName()
{
  return this->name;
}

void XMLElement::setName(string & name)
{
  this->name = name;
}

int XMLElement::getLineNr()
{
  return this->lineNr;
}

list<XMLElement> * XMLElement::getChildren()
{
  return &children;
}



string XMLElement::getStringAttribute(string & name)
{
  string noValue = "";
  return this->getStringAttribute(name, noValue);
}


string XMLElement::getStringAttribute(string & name,
				      string & defaultValue)
{
  string value = this->attributes[name];

  if (value.empty()) 
    {
      value = defaultValue;
    }

  return value;
}


int XMLElement::getIntAttribute(string & name)
{
  return this->getIntAttribute(name, 0);
}


int XMLElement::getIntAttribute(string & name,
				int defaultValue)
{
  string value = this->attributes[name];

  if (value.empty()) 
    {
      return defaultValue;
    } 
  else 
    {
      return StringManipulation::stringToInt(value);
    }
}


double XMLElement::getDoubleAttribute(string & name)
{
  return this->getDoubleAttribute(name, 0.);
}


double XMLElement::getDoubleAttribute(string & name,
				      double defaultValue)
{
  string value = this->attributes[name];

  if (value.empty()) 
    {
      return defaultValue;
    } 
  else 
    {
      return StringManipulation::stringToDouble(value);
    }
}


bool XMLElement::getBooleanAttribute(string & name,
				     string & trueValue,
				     string & falseValue,
				     bool defaultValue)
{
  string value = this->attributes[name];
 
  if (value.empty()) 
    {
      return defaultValue;
    } 
  else if (value == trueValue) 
    {
      return true;
    } 
  else if (value == falseValue) 
    {
      return false;
    } 
  else 
    {
      throw this->invalidValue(name, (string) value);
    }
}



void XMLElement::scanIdentifier(string & result)
{
  for (;;) 
    {
      unsigned char ch = this->readChar();
    
      if (((ch < 'A') || (ch > 'Z')) && ((ch < 'a') || (ch > 'z'))
	  && ((ch < '0') || (ch > '9')) && (ch != '_') && (ch != '.')
	  && (ch != ':') && (ch != '-')) //UNICODE && (ch <= '\u007E')) 
	{
	  this->unreadChar(ch);
	  return;
	}
      result += ch;
    }
}


unsigned char XMLElement::scanWhitespace()
{
  for (;;) 
    {
      unsigned char ch = this->readChar();
    
      switch (ch) 
	{
	case ' ':
	case '\t':
	case '\n':
	case '\r':
	  break;
	default:
	  return ch;
	}
    }
}


unsigned char XMLElement::scanWhitespace(string & result)
{
  for (;;) 
    {
      unsigned char ch = this->readChar();
    
      switch (ch) 
	{
	case ' ':
	case '\t':
	case '\n':
	  result += ch;
	case '\r':
	  break;
	default:
	  return ch;
	}
    }
}


void XMLElement::scanString(string & str)
{
  unsigned char delimiter = this->readChar();

  if ((delimiter != '\'') && (delimiter != '"')) 
    {
      throw this->expectedInput("' or \"");
    }

  for (;;) 
    {
      unsigned char ch = this->readChar();
    
      if (ch == delimiter) 
	{
	  return;
	} 
      else if (ch == '&') 
	{
	  this->resolveEntity(str);
	} 
      else 
	{
	  str += ch;
	}
    }
}


void XMLElement::scanPCData(string & data)
{
  for (;;) 
    {
      unsigned char ch = this->readChar();
    
      if (ch == '<') 
	{
	  ch = this->readChar();
	  
	  if (ch == '!') 
	    {
	      this->checkCDATA(data);
	    } 
	  else 
	    {
	      this->unreadChar(ch);
	      return;
	    }
	} 
      else if (ch == '&') 
	{
	  this->resolveEntity(data);
	} 
      else 
	{
	  data += ch; 
	}
    }
}


bool XMLElement::checkCDATA(string & buf)
{
  unsigned char ch = this->readChar();
  
  if (ch != '[') 
    {
      this->unreadChar(ch);
      this->skipSpecialTag(0);
      return false;
    } 
  else if (! this->checkLiteral("CDATA[")) 
    {
      this->skipSpecialTag(1); // one [ has already been read
      return false;
    } 
  else 
    {
      int delimiterCharsSkipped = 0;
      
      while (delimiterCharsSkipped < 3) 
	{
	  ch = this->readChar();

	  switch (ch) 
	    {
	    case ']':
	      if (delimiterCharsSkipped < 2) 
		{
		  delimiterCharsSkipped += 1;
		} 
	      else 
		{
		  buf += ']'; 
		  buf += ']'; 
		  delimiterCharsSkipped = 0;
		}
	      break;
	    case '>':
	      if (delimiterCharsSkipped < 2) 
		{
		  for (int i = 0; i < delimiterCharsSkipped; i++) 
		    {
		      buf += ']';
		    }

		  delimiterCharsSkipped = 0;
		  buf += '>'; 
		} 
	      else 
		{
		  delimiterCharsSkipped = 3;
		}

	      break;
	    default:
	      for (int i = 0; i < delimiterCharsSkipped; i += 1) 
		{
		  buf += ']'; 
		}

	      buf += ch; 
	      delimiterCharsSkipped = 0;
	    }
	}
      return true;
    }
}


void XMLElement::scanElement(XMLElement & elt)
{
  string name;
  this->scanIdentifier(name);
  elt.setName(name);
  unsigned char ch = this->scanWhitespace();

  while ((ch != '>') && (ch != '/')) 
    {
      string key;
      this->unreadChar(ch);
      this->scanIdentifier(key);
      ch = this->scanWhitespace();
      
      if (ch != '=') 
	{
	  throw this->expectedInput("=");
	}

      this->unreadChar(this->scanWhitespace());
      string value;
      this->scanString(value);
      elt.setAttribute(key, value);
      ch = this->scanWhitespace();
    }

  if (ch == '/') 
    {
      ch = this->readChar();
      if (ch != '>') 
	{
	  throw this->expectedInput(">");
	}
      return;
    }

  string buf;
  ch = this->scanWhitespace(buf);
  
  if (ch != '<') 
    {
      this->unreadChar(ch);
      this->scanPCData(buf);
    } 
  else 
    {
      for (;;) 
	{
	  ch = this->readChar();
	  if (ch == '!') 
	    {
	      if (this->checkCDATA(buf)) 
		{
		  this->scanPCData(buf);
		  break;
		} 
	      else 
		{
		  ch = this->scanWhitespace(buf);
		  if (ch != '<') 
		    {
		      this->unreadChar(ch);
		      this->scanPCData(buf);
		      break;
		    }
		}
	    } 
	  else 
	    {
	      buf.resize(0);
	      break;
	    }
	}
    }

  if (buf.length() == 0) 
    {
      while (ch != '/') 
	{
	  if (ch == '!') 
	    {
	      ch = this->readChar();
	
	      if (ch != '-') 
		{
		  throw this->expectedInput("Comment or Element");
		}

	      ch = this->readChar();
	      
	      if (ch != '-') 
		{
		  throw this->expectedInput("Comment or Element");
		}

	      this->skipComment();
	    } 
	  else 
	    {
	      this->unreadChar(ch);
	      XMLElement child = this->createAnotherElement();
	      this->scanElement(child);
	      elt.addChild(child);
	    }

	  ch = this->scanWhitespace();
	  
	  if (ch != '<') 
	    {
	      throw this->expectedInput("<");
	    }

	  ch = this->readChar();
	}
      this->unreadChar(ch);
    } 
  else 
    {
      //      if (this->ignoreWhitespace) 
      //	{
      //	  elt.setContent(buf.toString().trim());
      //	} 
      //      else 
      //	{
      elt.setContent(buf);
      //	}
    }

  ch = this->readChar();
  
  if (ch != '/') 
    {
      throw this->expectedInput("/");
    }

  this->unreadChar(this->scanWhitespace());

  if (! this->checkLiteral(name)) 
    {
      throw this->expectedInput(name);
    }

  if (this->scanWhitespace() != '>') 
    {
      throw this->expectedInput(">");
    }
}



void XMLElement::skipComment()
{
  int dashesToRead = 2;
  
  while (dashesToRead > 0) 
    {
      unsigned char ch = this->readChar();
    
      if (ch == '-') 
	{
	  dashesToRead -= 1;
	} 
      else 
	{
	  dashesToRead = 2;
	}
    }
  
  if (this->readChar() != '>') 
    {
      throw this->expectedInput(">");
    }
}



void XMLElement::skipSpecialTag(int bracketLevel)
{
  int tagLevel = 1; // <
  unsigned char stringDelimiter = '\0';
  
  if (bracketLevel == 0) 
    {
      unsigned char ch = this->readChar();

      if (ch == '[') 
	{
	  bracketLevel += 1;
	} 
      else 
	{
	  if (ch == '-') 
	    {
	      ch = this->readChar();

	      if (ch == '[') 
		{
		  bracketLevel += 1;
		} 
	      else if (ch == ']') 
		{
		  bracketLevel -= 1;
		} 
	      else if (ch == '-') 
		{
		  this->skipComment();
		  return;
		}
	    }
	}
  }

  while (tagLevel > 0) 
    {
      unsigned char ch = this->readChar();
      
      if (stringDelimiter == '\0') 
	{
	  if ((ch == '"') || (ch == '\'')) 
	    {
	      stringDelimiter = ch;
	    } 
	  else if (bracketLevel <= 0) 
	    {
	      if (ch == '<') 
		{
		  tagLevel += 1;
		} 
	      else if (ch == '>') 
		{
		  tagLevel -= 1;
		}
	    }
	  
	  if (ch == '[') 
	    {
	      bracketLevel += 1;
	    } 
	  else if (ch == ']') 
	    {
	      bracketLevel -= 1;
	    }
	} 
      else 
	{
	  if (ch == stringDelimiter) 
	    {
	      stringDelimiter = '\0';
	    }
	}
    }
}



bool XMLElement::checkLiteral(string literal)
{
  int length = (int)literal.length();

  for (int i = 0; i < length; i += 1) 
    {
      if (this->readChar() != literal.at(i)) 
	{
	  return false;
	}
    }
  return true;
}


void XMLElement::resolveEntity(string & buf)
{
  unsigned char ch = '\0';
  string key;
  
  for (;;) 
    {
      ch = this->readChar();
      
      if (ch == ';') 
	{
	  break;
	}
      
      key += ch; 
    }

  if (key.at(0) == '#') 
    {
      if (key.at(1) == 'x') 
	{
	  string numberStr = key.substr(2);
	  int intRep = StringManipulation::hexstringToInt(numberStr);
	  
	  if(intRep < 0 || intRep > 256)
	    {
	      throw unicodeError(intRep);
	    }

	  ch = (unsigned char) intRep;
	} 
      else 
	{
	  string numberStr = key.substr(1);
	  int intRep = StringManipulation::stringToInt(numberStr);

	  if(intRep < 0 || intRep > 256)
	    {
	      throw unicodeError(intRep);
	    }

	  ch = (unsigned char) intRep;
	}
      
      buf += ch;
    } 
  else 
    {
      string value = (*this->entities)[key]; 
      if (value.empty()) 
	{
	  throw this->unknownEntity(key);
	}
      buf.append(value);
    }
}


unsigned char XMLElement::readChar()
{
  if (this->charReadTooMuch != '\0') 
    {
      unsigned char ch = this->charReadTooMuch;
      this->charReadTooMuch = '\0';
      return ch;
    } 
  else 
    {
      int i = (*this->reader).get();
      if (i < 0) 
	{
	  throw this->unexpectedEndOfData();
	} 
      else if (i == 10) 
	{
	  this->parserLineNr += 1;
	  return '\n';
	} 
      else 
	{
	  return (unsigned char) i;
	}
    }
}


void XMLElement::unreadChar(unsigned char ch)
{
  this->charReadTooMuch = ch;
}


XMLElement XMLElement::createAnotherElement()
{
  return XMLElement(this->entities, this->ignoreWhitespace, false);
}


void XMLElement::singleLineWriter(ostream& writer)
{
  if (this->name.empty()) 
    {
      this->writeEncoded(writer, this->contents);
      return;
    }

  writer << '<';
  writer << this->name;

  if (! this->attributes.empty()) 
    {
      for(map<string,string>::iterator it=attributes.begin(); 
	  it != attributes.end(); ++it)
	{
	  writer << ' ';
	  string key   = it->first;
	  string value = it->second;
	  writer << key << '=' << '"';
	  this->writeEncoded(writer,value);
	  writer << '"';
	}
    }

  if (!this->contents.empty()) 
    {
      writer << '>';
      this->writeEncoded(writer, this->contents);
      writer << '<'; 
      writer << '/';
      writer << this->name;
      writer << '>';
    } 
  else if (this->children.empty()) 
    {
      writer << '/'; 
      writer << '>';
    } 
  else 
    {
      writer << '>';

      for(list<XMLElement>::iterator it=children.begin(); 
	  it != children.end(); ++it)
	{
	  it->singleLineWriter(writer);
	}

      writer << '<'; 
      writer << '/';
      writer << this->name;
      writer << '>';
  }
}

void XMLElement::indent(ostream & writer, int depth)
{
  for(int i=0; i<depth; i++)
    {
      writer << "    ";
    }
}

void XMLElement::write(string & file)
{
  ofstream writer(file.c_str());
  this->prettyWriter(writer);
  writer.close();
}

void XMLElement::prettyWriter(ostream& writer)
{
  this->prettyWriter(writer, /*depth*/ 0);
}

void XMLElement::prettyWriter(ostream& writer,int depth)
{
  if (this->name.empty()) 
    {
      this->writeEncoded(writer, this->contents);
      return;
    }

  this->indent(writer,depth);
  writer << '<';
  writer << this->name;

  if (! this->attributes.empty()) 
    {
      for(map<string,string>::iterator it=attributes.begin(); 
	  it != attributes.end(); ++it)
	{
	  writer << ' ';
	  string key   = it->first;
	  string value = it->second;
	  writer << key << '=' << '"';
	  this->writeEncoded(writer,value);
	  writer << '"';
	}
    }

  if (!this->contents.empty()) 
    {
      writer << '>';
      this->writeEncoded(writer, this->contents);
      writer << '<'; 
      writer << '/';
      writer << this->name;
      writer << '>';
      writer << '\n';
    } 
  else if (this->children.empty()) 
    {
      writer << '/'; 
      writer << '>';
      writer << '\n';
    } 
  else 
    {
      writer << '>';
      writer << '\n';

      for(list<XMLElement>::iterator it=children.begin(); 
	  it != children.end(); ++it)
	{
	  it->prettyWriter(writer,depth+1);
	}

      this->indent(writer,depth);
      writer << '<'; 
      writer << '/';
      writer << this->name;
      writer << '>';
      writer << '\n';
  }
}

void XMLElement::writeEncoded(ostream & writer,
			      string & str)
{
  for (unsigned int i = 0; i < str.length(); i += 1) 
    {
      unsigned char ch = str.at(i);

      switch (ch) 
	{
	case '<':
	  writer << "&lt;"; 
	  break;
	case '>':
	  writer << "&gt;";
	  break;
	case '&':
	  writer << "&amp;";
	  break;
	case '"':
	  writer << "&quot;";
	  break;
	case '\'':
	  writer << "&apos;";
	  break;
	default:
	  int unicode = (int) ch;

	  if ((unicode < 32) || (unicode > 126)) 
	    {
	      writer << "&#x" << StringManipulation::intToHexString(unicode);
	      writer << ";";
	    } 
	  else 
	    {
	      writer << ch;
	    }
	}
    }
}


void XMLElement::operator=(XMLElement & rhs)
{
  this->entities = rhs.entities;
  this->lineNr   = rhs.lineNr;
  this->ignoreWhitespace = rhs.ignoreWhitespace;
  this->name = rhs.name;
  this->contents = rhs.contents;

  this->attributes.clear();

  for(map<string,string>::iterator it=rhs.attributes.begin(); 
      it != rhs.attributes.end(); ++it)
    {
      this->attributes[it->first] = it->second;
    }

  this->children.clear();

  for(list<XMLElement>::iterator it=rhs.children.begin();
      it != rhs.children.end(); ++it)
    {
      XMLElement child = *it;
      this->children.push_back(child);
    }
}

bool XMLElement::operator==(XMLElement & rhs)
{
  // Check the easy fields
  if( this->name != rhs.name ) return false;
  if( this->contents != rhs.contents ) return false;
  
  // Check for equal attributes
  if( this->attributes.size() != rhs.attributes.size() ) return false;

  for(map<string,string>::iterator it=rhs.attributes.begin();
      it != rhs.attributes.end(); ++it)
    {
      if( this->attributes[it->first] != it->second ) return false;
    }

  // Check for equal children
  if( this->children.size() != rhs.children.size() ) return false;

  for(list<XMLElement>::iterator it1=rhs.children.begin();
      it1 != rhs.children.end(); ++it1)
    {
      bool childFound = false;

      for(list<XMLElement>::iterator it2=this->children.begin();
	  it2 != this->children.end(); ++it2)
	{
	  if( *it1 == *it2 ) childFound = true;
	}

      if( !childFound ) return false;
    }


  return true;
}

XMLParseException XMLElement::invalidValueSet(string name)
{
  string msg = "Invalid value set (entity name = \"" + name + "\")";
  return XMLParseException(this->getName(), this->parserLineNr, msg);
}


XMLParseException XMLElement::invalidValue(string name,
					 string value)
{
  string msg = "Attribute \"" + name + "\" does not contain a valid "
    + "value (\"" + value + "\")";
  return XMLParseException(this->getName(), this->parserLineNr, msg);
}


XMLParseException XMLElement::unexpectedEndOfData()
{
  string msg = "Unexpected end of data reached";
  return XMLParseException(this->getName(), this->parserLineNr, msg);
}


XMLParseException XMLElement::syntaxError(string context)
{
  string msg = "Syntax error while parsing " + context;
  return XMLParseException(this->getName(), this->parserLineNr, msg);
}

XMLParseException XMLElement::unicodeError(int value)
{
  string decUni = "#&" + StringManipulation::intToString(value) + ";";
  string hexUni = "#&x" + StringManipulation::intToHexString(value) + ";";

  string msg = "Unicode error while parsing " + decUni + " or " + hexUni
    + ". Only ASCII values (&#0; to &#256; or &#x0; to &#x100;) are "
    + "currently supported.";

  return XMLParseException(this->getName(), this->parserLineNr, msg);
}

XMLParseException XMLElement::expectedInput(string charSet)
{
  string msg = "Expected: " + charSet;
  return XMLParseException(this->getName(), this->parserLineNr, msg);
}


XMLParseException XMLElement::unknownEntity(string name)
{
  string msg = "Unknown or invalid entity: &" + name + ";";
  return XMLParseException(this->getName(), this->parserLineNr, msg);
}



