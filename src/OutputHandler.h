/*
 * OutputHandler.h
 *
 *  Created on: Mar 18, 2015
 *      Author: jordanwillis
> */

#ifndef SRC_OUTPUTHANDLER_H_
#define SRC_OUTPUTHANDLER_H_
#include <iostream>
#include <vector>

#include "utility.h"
#include "minicsv.h"


class OutputHandler{
	public:
		OutputHandler(csv::ofstream &);
		static const void WriteOutHeaders(csv::ofstream & Outputstream)
		{
			if(Outputstream.is_open())
				{
					for(std::vector<std::string>::iterator it = Utility::Headerlines.begin(); it != Utility::Headerlines.end(); ++it)
					{
				    	Outputstream << *it;
				    }
				    Outputstream << NEWLINE;
				}
		}
};
#endif /* SRC_OUTPUTHANDLER_H_ */
