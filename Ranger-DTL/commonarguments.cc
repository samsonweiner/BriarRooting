/*
 * *   Copyright (C) 2017 Mukul S. Bansal (mukul.bansal@uconn.edu).
 * *   Based on code originally written by Andre Wehe and Mukul S. Bansal
 * *
 * *   This program is free software: you can redistribute it and/or modify
 * *   it under the terms of the GNU General Public License as published by
 * *   the Free Software Foundation, either version 3 of the License, or
 * *   (at your option) any later version.
 * *
 * *   This program is distributed in the hope that it will be useful,
 * *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 * *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * *   GNU General Public License for more details.
 * *
 * *   You should have received a copy of the GNU General Public License
 * *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */



bool helpFlag = (Argument::find("-h") != NULL) || (Argument::find("--help") != NULL);
quiet = (Argument::find("-q") != NULL) || (Argument::find("--quiet") != NULL);
bool versionFlag = (Argument::find("-v") != NULL) || (Argument::find("--version") != NULL);
bool summary = (Argument::find("-s") != NULL) || (Argument::find("--summary") != NULL);


// read in arguments for the Duplication, Transfer, and Loss costs.
const Argument *argD = Argument::find("-D");
if (argD != NULL)
{
	argD->convert(DCOST);
}
	

const Argument *argT = Argument::find("-T");
if (argT != NULL)
{
	argT->convert(TCOST);
}


const Argument *argL = Argument::find("-L");
if (argL != NULL)
{
	argL->convert(LCOST);
}


const Argument *argType = Argument::find("--type");
if (argType != NULL)
{
	argType->convert(TYPE);
	if ((TYPE > 2) || (TYPE < 0))
		EXCEPTION("Invalid argument for --type");
}



const Argument *argThr = Argument::find("--thr");
if (argThr != NULL)
{
	argThr->convert(THRESHOLD);
}
else
{
	THRESHOLD = 10;
}


const Argument *argAdd = Argument::find("--add");
if (argAdd != NULL)
{
	argAdd->convert(ADD);
}
else
{
	ADD = TCOST/2;
}



// choose input
istream *in;
ifstream infile;
{
	const Argument *arg1 = Argument::find("-i");
	const Argument *arg2 = Argument::find("--input");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg == NULL) in = &cin;
	else {
		string filename;
		arg->convert(filename);
		msgout << "Input file: " << filename << endl;
		infile.open(filename.c_str());
		if (!infile) EXCEPTION("input file " << filename << " not found");
		in = &infile;
	}
}
Input input(in);

// choose output
ostream *out;
ofstream outfile;
{
	const Argument *arg1 = Argument::find("-o");
	const Argument *arg2 = Argument::find("--output");
	const Argument *arg = arg1 != NULL ? arg1 : arg2;
	if (arg == NULL) out = &cout;
	else {
		string filename;
		arg->convert(filename);
		msgout << "Output file: " << filename << endl;
		outfile.open(filename.c_str());
		if (!outfile) EXCEPTION("Output file " << filename << " not created");
		out = &outfile;
	}
}

// output the values of the D T and L costs.
msgout << "Duplication cost: " << DCOST << ", Transfer cost: " << TCOST << ", Loss cost: " << LCOST << endl;

if (TYPE > 0)
	msgout << "Using Type " << TYPE << " variable transfer costs. Threshold: " << THRESHOLD << ", Additional transfer cost: " << ADD << endl;
 
// output version number
if (versionFlag) {
	cout << version << endl;
	return 0;
}

// complain about unknown arguments
vector<Argument*> unused = Argument::unusedArgs();
for (vector<Argument*>::iterator itr = unused.begin(); itr != unused.end(); itr++) {
	EXCEPTION("Unknown argument " << (*itr)->key << endl << "Try " << argsv[0] << " --help for a list of arguments");
}
