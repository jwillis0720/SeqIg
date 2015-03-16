/*
 * Utility.h
 *
 *  Created on: Mar 15, 2015
 *      Author: jordanwillis
 */

#ifndef SRC_UTILITY_H_
#define SRC_UTILITY_H_

#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "DatabaseHandler.h"

namespace Utility {
	typedef std::map<seqan::CharString,char> Tcodonmap;
	Tcodonmap CreateCodonTable()
	{

		Tcodonmap m;
		m["TTT"] = 'F';
		m["TTC"] = 'F';
		m["TTA"] = 'L';
		m["TTG"] = 'L';
		m["CTT"] = 'L';
		m["CTC"] = 'L';
		m["CTA"] = 'L';
		m["CTG"] = 'L';
		m["ATT"] = 'I';
		m["ATC"] = 'I';
		m["ATA"] = 'I';
		m["ATG"] = 'M';
		m["GTT"] = 'V';
		m["GTC"] = 'V';
		m["GTA"] = 'V';
		m["GTG"] = 'V';
		m["TCT"] = 'S';
		m["TCC"] = 'S';
		m["TCA"] = 'S';
		m["TCG"] = 'S';
		m["CCT"] = 'P';
		m["CCC"] = 'P';
		m["CCA"] = 'P';
		m["CCG"] = 'P';
		m["ACT"] = 'T';
		m["ACC"] = 'T';
		m["ACA"] = 'T';
		m["ACG"] = 'T';
		m["GCT"] = 'A';
		m["GCC"] = 'A';
		m["GCA"] = 'A';
		m["GCG"] = 'A';
		m["TAT"] = 'Y';
		m["TAC"] = 'Y';
		m["TAA"] = '*';
		m["TAG"] = '*';
		m["CAT"] = 'H';
		m["CAC"] = 'H';
		m["CAA"] = 'Q';
		m["CAG"] = 'Q';
		m["AAT"] = 'N';
		m["AAC"] = 'N';
		m["AAA"] = 'K';
		m["AAG"] = 'K';
		m["GAT"] = 'D';
		m["GAC"] = 'D';
		m["GAA"] = 'E';
		m["GAG"] = 'E';
		m["TGT"] = 'C';
		m["TGC"] = 'C';
		m["TGA"] = '*';
		m["TGG"] = 'W';
		m["CGT"] = 'R';
		m["CGC"] = 'R';
		m["CGA"] = 'R';
		m["CGG"] = 'R';
		m["AGT"] = 'S';
		m["AGC"] = 'S';
		m["AGA"] = 'R';
		m["AGG"] = 'R';
		m["GGT"] = 'G';
		m["GGC"] = 'G';
		m["GGA"] = 'G';
		m["GGG"] = 'G';
		return m;
	}

	TCMap GetDatabaseFiles(const std::string & db_path, const bool & verbose){
	    //You could put any file path in here, e.g. "/home/me/mwah" to list that directory
	    boost::filesystem::path idbpath(db_path);
	    boost::filesystem::directory_iterator end_itr;
	    TCMap map_of_files;

	    // cycle through the directory
	    //This used to not segfault!!!
	    //TODO -- Segfaults here
	    for (boost::filesystem::directory_iterator itr(idbpath); itr != end_itr; ++itr)
	    {
	        // If it's not a directory, list it. If you want to list directories too, just remove this check.
	        if (boost::filesystem::is_regular_file(itr->path()))
	        {
	            // assign current file name to current_file and echo it out to the console.
	            std::string current_file = itr->path().string();

	            //Split by directory path and return a vector of the split
	            //TSVector split_lines = Split(current_file,'/');
	            TSVector split_lines;
	            boost::split(split_lines,current_file,boost::is_any_of("/"));

	            //Then we get the last thing in there to get the filename
	            std::string file_name = split_lines.back();

	            //Gene DB class will handle parsing this into memory
	            DatabaseHandler GeneDB(current_file);
	            //Try to open
	            try {
	                GeneDB.Open();
	                if(verbose)
	                {
	                    std::cout << "Loading DB at -> " << current_file << std::endl;
	                    GeneDB.PrintPretty();
	                }
	            }catch(DatabaseHandlerExceptions &msg){
	                std::cerr << "Couldn't open Database file" << std::endl;
	                std::cerr << msg.what();
	                exit(1);
	            }
	            seqan::CharString sub_file = file_name.substr(0, file_name.size()-6);
	            map_of_files[sub_file] = GeneDB.GetDbContainer();
	        }
	    }
	    return map_of_files;
	}

	//An function to quickly check if a file or directory exists
	bool check_if_dir_exists (const std::string &name) {
	    struct stat buffer;
	    return (stat (name.c_str(), &buffer) == 0);
	}

	//Codon Table when we need it
	Tcodonmap CodonTable = CreateCodonTable();
}
#endif /* SRC_UTILITY_H_ */
