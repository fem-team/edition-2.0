#pragma once

#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

class CPostOutputter
{
	ofstream OutputFile;

protected:

	CPostOutputter(string FileName);

	static CPostOutputter* _instance;

public:

	ofstream* GetOutputFile() 
		{ return &OutputFile; }

	static CPostOutputter* Instance(string FileName = " ");

	void OutputElementStress();

	template <typename T>
	CPostOutputter& operator<<(const T& item) 
	{
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	CPostOutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		op(OutputFile);
		return *this;
	}


};
