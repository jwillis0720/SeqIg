/*
 * OutputHandler.h
 *
 *  Created on: Mar 18, 2015
 *      Author: jordanwillis
 */

#ifndef SRC_OUTPUTHANDLER_H_
#define SRC_OUTPUTHANDLER_H_
#include "minicsv.h"
#include <iostream>
#include <vector>

class Output{
	public:
		Output();
		~Output();
		static const void WriteOutHeaders(const char * outputfilename)
		{
				std::vector<std::string> Headerlines =  {
						"Name", "RawQuerys","NucSequence","AASequence",
						"Top_V_Gene","Top_D_Gene","Top_J_Gene", "VScore",
						"DScore","JScore","FW1","FW2","FW3","CDR1","CDR2",
						"FW1_AA","FW2_AA","FW3_AA", "CDR1_AA","CDR2_AA"};
				csv::ofstream os(outputfilename, std::ios_base::out);
				os.set_delimiter(',');
				if(os.is_open())
				{
					for(std::vector<std::string>::iterator it = Headerlines.begin(); it != Headerlines.end(); ++it)
					{
				    	os << *it;
				    }
				    os << NEWLINE;
				}
		}
};


#endif /* SRC_OUTPUTHANDLER_H_ */
