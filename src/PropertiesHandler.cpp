/*
 * PropertiesHandler.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: jordanwillis
 */
#include "PropertiesHandler.h"
#include "SeqIgUtility.h"
#include "boost/algorithm/string.hpp"

PropertiesHandler::PropertiesHandler(std::string const & path):
	_ppath(path)
{
	    std::ifstream ifs(_ppath);
	    if(!ifs)throw PropertiesHandlerExceptions("Can't find properties file " + _ppath);
	    std::string line;
	    while(std::getline(ifs,line))
	    {
	    	TSVector split_lines;
	    	boost::split(split_lines,line,boost::is_any_of("\t"));
	    	if(split_lines[0] == "#Gene") continue;
	    	_VGenePropertiesContainer[split_lines[0]]["FR1s"] = atoi(split_lines[1].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["FR1e"] = atoi(split_lines[2].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["CDR1s"] = atoi(split_lines[3].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["CDR1e"] = atoi(split_lines[4].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["FR2s"] = atoi(split_lines[5].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["FR2e"] = atoi(split_lines[6].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["CDR2s"] = atoi(split_lines[7].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["CDR2e"] = atoi(split_lines[8].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["FR3s"] = atoi(split_lines[9].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["FR3e"] = atoi(split_lines[10].c_str());
	    	_VGenePropertiesContainer[split_lines[0]]["CDR3s"] = atoi(split_lines[10].c_str());
	    }


}

void PropertiesHandler::PrintPretty(){
	TProperties::iterator itr1;
	std::map<Tcs,int>::iterator itr2;
	for(itr1 = _VGenePropertiesContainer.begin(); itr1 != _VGenePropertiesContainer.end(); itr1++){
		std::cout << "\n" << itr1->first << "\t";
		for(itr2 = itr1->second.begin() ; itr2 != itr1->second.end(); itr2++){
			std::cout << itr2->first << "\t" << itr2->second << "\t";
		}
	}
}
