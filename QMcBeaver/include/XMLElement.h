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

#ifndef XMLElement_H
#define XMLElement_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <list>

#include "XMLParseException.h"

using namespace std;


/**
  XMLElement is a representation of an XML object. The object is able to parse
  and write XML code.
 */

class XMLElement
{
  /**
    The attributes given to the element.
   */
private: map<string,string> attributes;

/**
  Child elements of the element.
  */
private: list<XMLElement> children;

/**
  The name of the element.
  */
private: string name;


/**
  The #PCDATA content of the object.
  */
private: string contents;


/**
  Conversion table for "&amp;...;" entities. The keys are the entity names
  without the & and ; delimiters.
  */
private: map<string,string> * entities;

/**
  A concrete instance of the entities conversion table.  If the object
  is constructed without providing a value for entities, this will be
  used.  By doing it this way, some speed and memory usage are sacrificed,
  but it takes care of it's own memory deallocation.
  */
private: map<string,string> entitiesInstance;


/**
  The line number where the element starts.
  */
private: int lineNr;

/**
  The current line number in the source content.
  */
private: int parserLineNr;

/**
  Character read too much.
  This character provides push-back functionality to the input stream
  without having to use a pushback stream.
  If there is no such character, this field is '\0'.
  */
private: unsigned char charReadTooMuch;

/**
  The input stream provided by the caller of the parse method.
  */
private: istream * reader;


/**
  <code>true</code> if the leading and trailing whitespace of #PCDATA
  sections have to be ignored.
  */
private: bool ignoreWhitespace;



/**
  Creates and initializes a new XML element.
  A basic entity ("&amp;", etc.) conversion table is used and leading 
  whitespace is not skipped.
  */
public: XMLElement();
    
/**
  Creates and initializes a new XML element.
  A basic entity ("&amp;", etc.) conversion table and the provided entity
  conversion table are used and leading whitespace
  is not skipped.
  
  @param entities The entity conversion table.
  */
public: XMLElement(map<string,string> * entities);

/**
  Creates and initializes a new XML element.
  A basic entity ("&amp;", etc.) conversion table is used and skipping of 
  leading whitespace is controled by <code>skipLeadingWhitespace</code>.
  
  @param skipLeadingWhitespace <code>true</code> if leading and trailing 
  whitespace in PCDATA content has to be removed.
  */
public: XMLElement(bool skipLeadingWhitespace);

/**
  Creates and initializes a new XML element.
  A basic entity ("&amp;", etc.) conversion table and the provided entity
  conversion table are used and leading whitespace is controled by 
  <code>skipLeadingWhitespace</code>.
  
  @param entities The entity conversion table.
  @param skipLeadingWhitespace <code>true</code> if leading and trailing 
  whitespace in PCDATA content has to be removed.
  */
public: XMLElement(map<string,string> * entities, bool skipLeadingWhitespace);


/**
  Creates and initializes a new XML element.
  A basic entity ("&amp;", etc.) conversion table can be provided by setting
  <code>fillBasicConversionTable</code> to <code>true</code>. The 
  provided entity conversion table is used and skipping leading whitespace is 
  controled by <code>skipLeadingWhitespace</code>.  
  <P>
  This constructor should <I>only</I> be called from
  {@link #createAnotherElement() createAnotherElement}
  to create child elements.

  @param entities The entity conversion table.
  @param skipLeadingWhitespace <code>true</code> if leading and trailing 
  whitespace in PCDATA content has to be removed.
  @param fillBasicConversionTable <code>true</code> if the basic entities 
  need to be added to the entity list.
     */
private: XMLElement(map<string,string> * entities, bool skipLeadingWhitespace,
		    bool fillBasicConversionTable);

/**
  Initializes a new XML element.
  
  @param skipLeadingWhitespace
  <code>true</code> if leading and trailing whitespace in PCDATA
  content has to be removed.
  @param fillBasicConversionTable
  <code>true</code> if the basic entities need to be added to
  the entity list.
  */ 
private: void initialize(bool skipLeadingWhitespace, 
			 bool fillBasicConversionTable);

/**
  Initializes a new XML element.
  
  @param entities
  The entity conversion table.
  @param skipLeadingWhitespace
  <code>true</code> if leading and trailing whitespace in PCDATA
  content has to be removed.
  @param fillBasicConversionTable
  <code>true</code> if the basic entities need to be added to
  the entity list.
  */ 
private: void initialize(map<string,string> * entities,
			 bool skipLeadingWhitespace,
			 bool fillBasicConversionTable);

/**
  Returns the number of child elements of the element.

  @return number of child elements.
  */
public: int countChildren();

/**
  Adds a child element.
  
  @param child The child element to add.
  */
public: void addChild(XMLElement & child);

/**
  Adds or modifies an attribute.
  
  @param name The name of the attribute.
  @param value The value of the attribute.
  */
public: void setAttribute(string & name, string & value);

/**
  Adds or modifies an attribute.
  
  @param name The name of the attribute.
  @param value The value of the attribute.
  */
public: void setAttribute(string & name, int value);

/**
  Adds or modifies an attribute.
  
  @param name The name of the attribute.
  @param value The value of the attribute.
  */
public: void setAttribute(string & name, double value);
 
/**
  Reads one XML element from a file and parses it.
  
  @param file The file from which to retrieve the XML data.

  @throws XMLParseException If an error occured while parsing the read data.
  */
public: void parse(string & file);

   
/**
  Reads one XML element from a stream and parses it.
  
  @param reader The stream from which to retrieve the XML data.

  @throws XMLParseException If an error occured while parsing the read data.
  */
public: void parse(istream & reader);

/**
  Reads one XML element from a stream and parses it.
  
  @param reader The stream from which to retrieve the XML data.
  @param startingLineNr The line number of the first line in the data.

  @throws XMLParseException If an error occured while parsing the read data.
  */
private: void parse(istream & reader, int startingLineNr);

/**
  Removes a child element.
  
  @param child The child element to remove.
  */
public: void removeChild(XMLElement & child);

/**
  Returns the child elements as a Vector. It is safe to modify this
  Vector.

  @return The child elements of this element.
  */
public: list<XMLElement> * getChildren();


/**
  Returns an attribute of the element.
  If the attribute doesn't exist, an empty string is returned.
  
  @param name The name of the attribute.

  @return The value of the attribute.
  */
public: string getStringAttribute(string & name);

/**
  Returns an attribute of the element.
  If the attribute doesn't exist, <code>defaultValue</code> is returned.
  
  @param name The name of the attribute.
  @param defaultValue Key to use if the attribute is missing.
  
  @return The value of the attribute.
  */
public: string getStringAttribute(string & name, string & defaultValue);

/**
  Returns an attribute of the element.
  If the attribute doesn't exist, <code>0</code> is returned.
  
  @param name The name of the attribute.

  @return The value of the attribute.
  */
public: int getIntAttribute(string & name);

/**
  Returns an attribute of the element.
  If the attribute doesn't exist, <code>defaultValue</code> is returned.
  
  @param name The name of the attribute.
  @param defaultValue Key to use if the attribute is missing.

  @return The value of the attribute.
  */
public: int getIntAttribute(string & name, int defaultValue);

/**
  Returns an attribute of the element.
  If the attribute doesn't exist, <code>0.0</code> is returned.
  
  @param name The name of the attribute.

  @return The value of the attribute.
  */
public: double getDoubleAttribute(string & name);

/**
  Returns an attribute of the element.
  If the attribute doesn't exist, <code>defaultValue</code> is returned.
  
  @param name The name of the attribute.
  @param defaultValue Key to use if the attribute is missing.

  @return The value of the attribute.
  */
public: double getDoubleAttribute(string & name, double defaultValue);

/**
  Returns an attribute of the element.
  If the attribute doesn't exist, <code>defaultValue</code> is returned.
  If the value of the attribute is equal to <code>trueValue</code>,
  <code>true</code> is returned.
  If the value of the attribute is equal to <code>falseValue</code>,
  <code>false</code> is returned.
  If the value doesn't match <code>trueValue</code> or
  <code>falseValue</code>, an exception is thrown.
  
  @param name         The name of the attribute.
  @param trueValue    The value associated with <code>true</code>.
  @param falseValue   The value associated with <code>true</code>.
  @param defaultValue Value to use if the attribute is missing.

  @return The value of the attribute.

  @throws XMLParseException   If the value doesn't match 
  <code>trueValue</code> or <code>falseValue</code>.
  */
public: bool getBooleanAttribute(string & name, string & trueValue,
				 string & falseValue, bool defaultValue);

/**
  Removes an attribute.
  
  @param name The name of the attribute.
  */
public: void removeAttribute(string & name);

/**
  Changes the content string.
  
  @param content The new content string.
  */
public: void setContent(string & content);

/**
  Returns the PCDATA content of the object. If there is no such content,
  an empty string is returned.

  @return PCDATA content.
  */
public: string getContent();

/**
  Returns the name of the element.

  @return name of the element.
  */
public: string getName();

/**
  Changes the name of the element.
  
  @param name The new name.
  */
public: void setName(string &name);

/**
  Returns the line number in the source data on which the element is found.
  This method returns <code>0</code> there is no associated source data.

  @return Line number in the source data on which the element is found.
  */
public: int getLineNr();

/**
  Scans an identifier from the current stream.
  The scanned identifier is appended to <code>result</code>.
  
  @param result The buffer in which the scanned identifier will be put.
  */
private: void scanIdentifier(string & result);

/**
  This method scans an identifier from the current stream.
  
  @return the next character following the whitespace.
  */
private: unsigned char scanWhitespace();

/**
  This method scans an identifier from the current stream.
  The scanned whitespace is appended to <code>result</code>.
  
  @param result Buffer to which the scanned whitespace is appended.

  @return the next character following the whitespace.
  */
private: unsigned char scanWhitespace(string & result);

/**
  This method scans a delimited string from the current stream.
  The scanned string without delimiters is appended to <code>str</code>.

  @param str Buffer to which a scaned delimited string is added.
  */
private: void scanString(string & str);

/**
  Scans a #PCDATA element. CDATA sections and entities are resolved.
  The next &lt; char is skipped.
  The scanned data is appended to <code>data</code>.

  @param data Scaned PCDATA is added to this buffer.
  */
private: void scanPCData(string & data);

/**
  Scans an XML element.
  
  @param elt The element that will contain the result.
  */
private: void scanElement(XMLElement & elt);

/**
  Scans a special tag and if the tag is a CDATA section, append its
  content to <code>buf</code>.
  
  @param buf Buffer to add scaned data to.
  */
private: bool checkCDATA(string & buf);

/**
  Skips a comment.
  */
private: void skipComment();

/**
  Skips a special tag or comment.
  
  @param bracketLevel The number of open square brackets ([) that have
  already been read.
  */
private: void skipSpecialTag(int bracketLevel);

/**
  Scans the data for literal text.
  Scanning stops when a character does not match or after the complete
  text has been checked, whichever comes first.
  
  @param literal the literal to check.
  */
private: bool checkLiteral(string literal);

/**
  Resolves an entity. The name of the entity is read from the current stream.
  The value of the entity is appended to <code>buf</code>.
  
  @param buf Where to put the entity value.
  */
private: void resolveEntity(string & buf);

/**
  Reads a character from a stream..

  @return Character from the stream.
 */
private: unsigned char readChar();

/**
  Pushes a character back to the read-back buffer.
  
  @param ch The character to push back.
  */
private: void unreadChar(unsigned char ch);

/**
  Creates a new similar XML element.
  <P>
  You should override this method when subclassing XMLElement.

  @return a new XMLElement.
  */
private: XMLElement createAnotherElement();

/**
  Writes the XML element to an output stream as a single line.
  
  @param writer The stream to write the XML data to.
  */
public: void singleLineWriter(ostream & writer);

/**
  Provides the proper indentation spacing when pretty printing an element
  depth deep in the XML tree.

  @param writer The stream to write the XML data to.
  @param depth Depth to indent to.
  */
private: void indent(ostream & writer, int depth);

/**
  Writes the XML element to a file using a pretty format.
  
  @param file The file to write the XML data to.
  */
public: void write(string & file);


/**
  Writes the XML element to an output stream using a pretty format.
  
  @param writer The stream to write the XML data to.
  */
public: void prettyWriter(ostream & writer);

/**
  Writes the XML element to an output stream using a pretty format.  The 
  indentation begins depth deep.
  
  @param writer The stream to write the XML data to.
  @param depth Depth to begin indentation at.
  */
private: void prettyWriter(ostream & writer,int depth);


/**
  Writes a string encoded to a stream.
  
  @param writer The stream to write the XML data to.
  @param str The string to write encoded.
     */
private: void writeEncoded(ostream & writer, string & str);

/**
  Sets two objects equal to one another.

  @param rhs object to set this object equal to.
  */
public: void operator=(XMLElement & rhs);

/**
  Determines if two objects equal to one another.

  @param rhs object to determine if this one is equal to.

  @return <code>true</code> if both objects are equal and <code>false</code>
  otherwise.
  */
public: bool operator==(XMLElement & rhs);

/**
  Creates a parse exception for when an invalid valueset is given to
  a method.
  
  @param name The name of the entity.

  @return An exception to throw.
  */
private: XMLParseException invalidValueSet(string name);

/**
  Creates a parse exception for when an invalid value is given to a
  method.
  
  @param name  The name of the entity.
  @param value The value of the entity.

  @return An exception to throw.
  */
private: XMLParseException invalidValue(string name, string value);

/**
  Creates a parse exception for when the end of the data input has been
  reached.

  @return An exception to throw.
  */
private: XMLParseException unexpectedEndOfData();

/**
  Creates a parse exception for when a syntax error occured.
  
  @param context The context in which the error occured.

  @return An exception to throw.
  */
private: XMLParseException syntaxError(string context);

/**
  Creates a parse exception for when an unsupported unicode value 
  is encountered.
  
  @param value integer value of the unsupported unicode representation.

  @return An exception to throw.
  */ 
private: XMLParseException unicodeError(int value);

/**
  Creates a parse exception for when the next character read is not
  the character that was expected.
  
  @param charSet The set of characters (in human readable form) that was
  expected.

  @return An exception to throw.
  */
private: XMLParseException expectedInput(string charSet);

/**
  Creates a parse exception for when an entity could not be resolved.
  
  @param name The name of the entity.

  @return An exception to throw.
  */
private: XMLParseException unknownEntity(string name);
};

#endif
