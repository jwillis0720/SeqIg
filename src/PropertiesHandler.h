/*
 * PropertiesHandler.h
 *
 *  Created on: Mar 15, 2015
 *      Author: jordanwillis
 */

#ifndef SRC_PROPERTIESHANDLER_H_
#define SRC_PROPERTIESHANDLER_H_

#include <sstream>
#include <stdio.h>
#include <seqan/sequence.h>
#include "SeqIg.h"
#include "StructDefs.h"
#include "boost/algorithm/string.hpp"



class PropertiesHandler
{
private :
    std::string _ppath;
    TProperties _VGenePropertiesContainer;

public:
    PropertiesHandler(std::string const &);
    ~PropertiesHandler(){};
    void Open();
    void PrintPretty();
};

class PropertiesHandlerExceptions : public std::exception
{
private:
    std::string err_msg;
public:
    PropertiesHandlerExceptions(const std::string msg) : err_msg(msg) {};
    ~PropertiesHandlerExceptions() throw() {};
    const char *what() const throw() { return this->err_msg.c_str(); };
};

#endif /* SRC_PROPERTIESHANDLER_H_ */
