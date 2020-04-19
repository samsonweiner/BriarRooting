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



int DCOST = 2;  // duplication cost
int TCOST = 3;  // transfer cost
int LCOST = 1;  // losst cost
int TYPE = 0;   // type of transfer cost to be used (takes values 0, 1, or 2)
int THRESHOLD = 0;  // threshold for variable transfer costs
int ADD = 0;   // additional transfer cost for variable transfer costs


#include <signal.h>
#include "common.h"
#include "argument.h"
#include "input.h"
#include "node.h"
#include "tree.h"
#include "DTL-treeset.h"
#include "DTL-algorithm.h"

/*
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

boost::random::mt19937 gen;   // initialize boosts random number generator and make it globally accessible
*/



// long int getRandom(long int limit);



int main(int argc, char *argsv[]) {
	initialization();
	signal(SIGINT, interruptFunction);

	Argument::add(argc, argsv);
	
	
	
	// random seed given
	unsigned int randomseed = unsigned(time(NULL)) % 0xFFFF;
	{
		const Argument *arg = Argument::find("--seed");
		if (arg != NULL) arg->convert(randomseed);
		srand(randomseed);
	}
	
	#include "commonarguments.cc"

	// random seed
	msgout << "Random seed: " << randomseed << endl;
	
	
	// output help message
	if (helpFlag) {
		cout << "Usage: " << argsv[0] << " [ARGUMENT]" << endl;
		cout << endl;
		cout << "  -i, --input                 input file" << endl;
		cout << "  -o, --output                output file" << endl;
		cout << "  -D,                         Duplication cost (whole number only, default value 2)" << endl;
		cout << "  -T,                         Transfer cost (whole number only, default value 3)" << endl;
		cout << "  -L,                         Loss cost (whole number only, default value 1)" << endl;
		cout << "      --type  0|1|2           type of transfer cost to use" << endl;      
		cout << "                                0 - Fixed transfer cost [default]" << endl;
		cout << "                                1 - Simple threshold based variable transfer cost" << endl;
		cout << "                                2 - General threshold based variable transfer cost" << endl;
		cout << "      --thr <whole number>    threshold value for use with types 1 and 2 above (default 10)" << endl;	
		cout << "      --add <whole number>    additional transfer cost for use with types 1 and 2 above" << endl;
		cout << "                              (default value is transfer cost divided by 2, rounded down)" << endl;
		cout << "      --seed <whole number>   set a user defined random number generator seed." << endl;
		cout << "  -q, --quiet                 no process output" << endl;
		cout << "  -s, --summary               only output summary statistics" << endl;
		cout << "  -v, --version               version number" << endl;
		cout << "  -h, --help                  brief help message" << endl;
		return 0;
	}

	
	

	// create the reconciliation
	DTLreconciliation *reconciliation = new DTLreconciliation();

	// read input trees
	reconciliation->readTrees(input);

	
	// run timed
	time_t startTime = time(NULL);
	reconciliation->run(*out, ALL, summary);
	time_t endTime = time(NULL);

	delete reconciliation;

	// timing
	{
		long unsigned int duration = endTime - startTime;
		const int days = duration / (60 * 60 * 24); duration -= days * (60 * 60 * 24);
		const int hours = duration / (60 * 60); duration -= hours * (60 * 60);
		const int minutes = duration / 60; duration -= minutes * 60;
		const int seconds = duration;
		ostringstream os;
		bool on = false;
		if (days!=0) {on = true; os << days << "d ";}
		if (on || (hours!=0)) {on = true; os << hours << "h ";}
		if (on || (minutes!=0)) {on = true; os << minutes << "m ";}
		os << seconds << 's';
		msgout << "\n"<<"time: " << os.str() << endl;
	}

	#ifdef DEBUG
	for (map<string,int>::iterator itr=NodeID::namelist.begin();itr!=NodeID::namelist.end();itr++)
		cout << "Species namelist entry:" << itr->first << ' ' << itr->second << endl;
	msgout << "\nAll done!" << endl;
	#endif

	return 0;
}


