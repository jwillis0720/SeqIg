//
//  SeqIg.h
//  seqan
//
//  Created by Jordan Willis on 12/7/14.
//
//
#ifndef seqan_SeqIg_h
#define seqan_SeqIg_h

//Seqan Headers
#include "seqan/arg_parse.h"
#include "seqan/seq_io.h"

//Project Headers
#include "DatabaseHandler.h"
#include "StructDefs.h"
#include "PropertiesHandler.h"


void SetUpArgumentParser(seqan::ArgumentParser &);
void SetDatabaseFastas(SeqIgOptions const &options, DatabasePaths &);
seqan::ArgumentParser::ParseResult ExtractOptions(seqan::ArgumentParser const &, SeqIgOptions &);
#endif
