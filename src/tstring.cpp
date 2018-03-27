/** @file tstring.cpp
 *
 *   Copyright (C) 2008 Samuel Neuenschwander <samuel.neuenschwander@unil.ch>
 *
 *   quantiNemo:
 *   quantiNemo is an individual-based, genetically explicit stochastic
 *   simulation program. It was developed to investigate the effects of
 *   selection, mutation, recombination, and drift on quantitative traits
 *   with varying architectures in structured populations connected by
 *   migration and located in a heterogeneous habitat.
 *
 *   quantiNemo is built on the evolutionary and population genetics
 *   programming framework NEMO (Guillaume and Rougemont, 2006, Bioinformatics).
 *
 *
 *   Licensing:
 *   This file is part of quantiNemo.
 *
 *   quantiNemo is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   quantiNemo is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with quantiNemo.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "tstring.h"
#include "functions.h"
#include <fstream>
using namespace std;

//------------------------------------------------------------------------------
/** read white space and comment until the end of the file, or the specified character "end".
 * returns FALSE if the "end "character is first reached
 * returns TRUE if valuable text is reached before the "end" character
 *
 * types of comments:
 * 	#: single line (until the end of the line)
 * 	#/ ... /# commented parts (may be among several lines)
 */
bool STRING::removeCommentAndSpace(istream & IN, unsigned int& line,
		int comment, char end, bool removed)
{
	char c;

	switch (comment) {
		case 0: // no comment
			while (IN.get(c)) {
				if (c == end) {
					if (isEOL(c, IN)) ++line;
					break;
				}
				if (isspace(c)) {
					if (isEOL(c, IN)) {
						++line;
						if (end == '\n') break;
					}
					continue;
				}
				if (c == '#') {      // start of a comment
					if (IN.peek() == '/') {
						IN.get();
						return removeCommentAndSpace(IN, line, 2, end);
					}
					else return removeCommentAndSpace(IN, line, 1, end);
				}

				// no comment
				IN.putback(c);  // put back the character and return
				if (removed) IN.putback(' '); // if it was a comment 2 control that there is at least a space
				return true;
			}
			return false;

		case 1: // within comment #
			while (IN.get(c)) {
				if (isEOL(c, IN)) {
					++line;
					if (end == '\n') return false; // internally the end of line is always '\n'
					return removeCommentAndSpace(IN, line, 0, end);
				}
			}
			return false;

		case 2: // within comment #/ ... /#
			while (IN.get(c)) {
				if (isEOL(c, IN)) ++line;
				if (c == '/' && IN.peek() == '#') {
					IN.get();     // remove #
					return removeCommentAndSpace(IN, line, 0, end, true);
				}
			}
			return false;
	}
	return false;     // will never be used, but the compiler prefers it...
}

//------------------------------------------------------------------------------
/** file information starts with the key word "[FILE_INFO]" and the information
 * is enclosed by {...}
 */
map<string, int> STRING::readFileInfo(istream & IN, unsigned int& line)
{
	char c;
	string key, keyWord = "[FILE_INFO]";
	int value, i;
	map<string, int> info;

	// check if it is really a file info
	for (i = 0; i < (int) keyWord.size(); ++i) {
		IN.get(c);
		if (c != keyWord[i]) {
			// it is not a file info: put all read characters back
			IN.putback(c);     // current character
			while (i >= 0) {     // all previous read characters
				IN.putback(keyWord[i]);
				--i;
			}
			return info;
		}
	}

	// get the opening bracket
	if (!removeCommentAndSpace(IN, line))
		error("File information is not complete!\n");
	IN.get(c);
	if (c != '{')
		error("File info could not be properly read: misplaced '{'!\n");
	if (!removeCommentAndSpace(IN, line))
		error("File information is not complete!\n");

	// read the info line by line
	while (IN.good() && !IN.eof()) {
		if (!removeCommentAndSpace(IN, line))
			error("File information is not complete1!\n");
		if (IN.peek() == '}') break;

		//read the parameter name:
		key = "";
		while (IN.get(c) && !isspace(c) && !isEOL(c, IN)) {
			switch (c) {
				case '}':
					error("File info could not be read properly!\n");
					break;
				case '#':
					IN.putback(c);
					if (!removeCommentAndSpace(IN, line))
						error("File information is not complete!\n");
					break;
				default:
					key += c;
					break;
			}
		}
		if (isEOL(c)) error("File info could not be properly read!\n");

		// read the argument (is an integer!)
		removeCommentAndSpace(IN, line);
		IN >> value;

		// save the file info to a map
		if (info.find(key) != info.end())
			warning("Column '%s' has previously been set!\n", key.c_str());
		info[key] = value;
	}

	// remove the }
	assert(IN.peek() == '}');
	IN.get();

	return info;
}

//------------------------------------------------------------------------------
/** read a matrix and return it as a string
 * a matrix is read until the number of opening and closing brackets is the same
 * comments are removed
 * no control if numbers are read, or text is between rows */
string STRING::readBrackets2String(istream & IN, unsigned int& linecnt,
		char dim)
{
	char c;
	int whiteSpace = -1; // -1: ignore coming ws, 0: last char was important, 1: ws needed
	string args;

	// read the first char
	IN.get(c);
	args += c;
	assert(c=='{' || c=='(' || c=='[');

	// read brackets
	while (IN.get(c) && IN.good() && !IN.eof()) {
		switch (c) {
			case '{':
				IN.putback(c);
				args += readBrackets2String(IN, linecnt, '}');
				continue;
			case '(':
				IN.putback(c);
				args += readBrackets2String(IN, linecnt, ')');
				continue;
			case '[':
				IN.putback(c);
				args += readBrackets2String(IN, linecnt, ']');
				continue;
		}

		if (c == dim) {
			args += c;
			return args;
		}

		if (isspace(c)) {
			if (isEOL(c, IN)) ++linecnt;
			if (whiteSpace != -1) whiteSpace = 1;
			continue;
		}

		if (c == '#') {
			IN.putback(c);
			if (!removeCommentAndSpace(IN, linecnt))
				throw "Could not read text between brackets!";
			if (whiteSpace != -1) whiteSpace = 1;
			continue;
		}

		// that character is needed, add it to the string
		if (whiteSpace == 1) args += " ";
		whiteSpace = 0;
		args += c;
	}

	throw "Could not read text between brackets: missing closing bracket!";
}

//------------------------------------------------------------------------------
/** read the stream  until one of the specified characters or any space occur,
 * parameter include: 0: found char is put back to the stream (default)
 *                    1: found char is removed and also not returned
 *                    2: found char is added to the returned string, but removed from the stream
 */
string STRING::find_first_of(istream & IN, string chars, bool stopAtSpace,
		unsigned int include)
{
	unsigned int i, size = (unsigned int)chars.size();
	string arg;
	char c;
	while (IN.get(c)) {
		if (stopAtSpace && isspace(c)) break;
		for (i = 0; i < size; ++i) {
			if (c == chars[i]) {
				break;
			}
		}
		arg += c;
	}
	switch (include) {
		case 0:
			IN.putback(c);
			break;
		case 1:
			break;
		case 2:
			arg += c;
			break;
	}
	return arg;
}

//------------------------------------------------------------------------------
/** read string until the specified character occurs. The starting and ending
 * character may be included or excluded and all comments on the way are removed.
 * The last character "item" is removed from the stream
 */
bool STRING::readUntilCharacter(istream & IN, string& output,
		unsigned int& linecnt, char item, bool includeFirst, bool includeLast)
{
	output = readUntilCharacter(IN, linecnt, item, includeFirst, includeLast);
	if (IN.eof() && output.empty()) return false;
	return true;
}

string STRING::readUntilCharacter(istream & IN, unsigned int& linecnt,
		char item, bool includeFirst, bool includeLast)
{
	char c;
	string args;
	int whiteSpace = 0;

	if (IN.get(c)) {
		if (includeFirst) args = c;

		// read brackets
		while (IN.get(c)) {
			if (c == item || (item == '\n' && isEOL(c, IN))) break;
			if (c == '#') {	// start of a comment
				IN.putback(c);
				if (!removeCommentAndSpace(IN, linecnt, 0, item)) break;
				if (!whiteSpace) args += " ";
			}
			else if (isspace(c)) { // reduce space to a single space
				if (!whiteSpace) {
					++whiteSpace;
					args += " ";
				}
			}
			else {
				whiteSpace = 0;
				args += c;
			}
		}

		//remove the trailing space
		int i;
		for (i = (int) args.size(); i > 0; --i) {
			if (!isspace(args[i - 1])) break;
		}
		args = args.substr(0, i);

		if (includeLast) args += item;
	}
	return args;
}

//------------------------------------------------------------------------------
/** function checks if the file exists */
bool STRING::file_exists(const string& name)
{
	ifstream ifs(name.c_str());
	if (!ifs) return false;
	ifs.close();
	return true;
}

//------------------------------------------------------------------------------
/** in this version the name of the file will first be tested and if the file
 * is not available it wil lbe searched at the given directory
 */
bool STRING::file_exists(string& name, const string& dir)
{
	if (file_exists(name)) return true;
	name = dir + name;
	return file_exists(name);
}

//------------------------------------------------------------------------------
/** Remove extension if present. */
string STRING::remove_file_extension(string file)
{
    const size_t period_idx = file.rfind('.');
    if (std::string::npos == period_idx) return file;
    return  file.substr(0, period_idx);
}

//------------------------------------------------------------------------------
/** Remove directory if present. */
string STRING::remove_file_path(string file)
{
    const size_t last_slash_idx = file.find_last_of("\\/");
    if (std::string::npos == last_slash_idx) return file;
    return file.substr(last_slash_idx + 1);
}

//------------------------------------------------------------------------------
/** Return the directory if present. */
string STRING::get_file_path(string file)
{
    const size_t last_slash_idx = file.find_last_of("\\/");
    if (std::string::npos == last_slash_idx) return "./";
    return file.substr(0,last_slash_idx);
}

//------------------------------------------------------------------------------
/** Return the file if present, without the path. */
string STRING::get_file_name(string file)
{
    return remove_file_path(file);
}

//------------------------------------------------------------------------------
/** Remove directory and extension if present. */
string STRING::get_file_basename(string file)
{
    return remove_file_extension(remove_file_path(file));
}
