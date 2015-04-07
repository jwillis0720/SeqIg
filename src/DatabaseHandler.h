//
//  DatabaseIO.h
//  seqan_sandbox
//
//  Created by Jordan Willis on 12/6/14.
//
//

#include <stdio.h>
#include <seqan/sequence.h>
#include "SeqIg.h"
#include "StructDefs.h"


#ifndef __DatabaseHandler__
#define __DatabaseHandler__



class DatabaseHandler
{
private :
    const char * _dbname;
    seqan::StringSet<seqan::CharString> _ids;
    seqan::StringSet<seqan::Dna5String> _seqs;
    Tdbcontainer _dbcontainer;
    int _rresult;
public:
    DatabaseHandler(std::string const &);
    ~DatabaseHandler(){};
    void Open();
    void PrintPretty();
    seqan::StringSet<seqan::Dna5String> GetAllSeqs();
    Tdbcontainer GetDbContainer();
};

class DatabaseHandlerExceptions : public std::exception
{
private:
    std::string err_msg;
public:
    DatabaseHandlerExceptions(const std::string msg) : err_msg(msg) {};
    ~DatabaseHandlerExceptions() throw() {};
    const char *what() const throw() { return this->err_msg.c_str(); };
};

#endif /* defined(__seqan_sandbox__DatabaseHandler__) */
