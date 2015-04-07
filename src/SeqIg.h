//
//  SeqIg.h
//  seqan
//
//  Created by Jordan Willis on 12/7/14.
//
//
#ifndef seqan_SeqIg_h
#define seqan_SeqIg_h
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/arg_parse.h>
#include <sstream>
#include "StructDefs.h"
#include "DatabaseHandler.h"
#include "AlignAntibody.h"
#include "AntibodyJunction.h"
#include "PropertiesHandler.h"
#include "OutputHandler.h"
#include "minicsv.h"
#include "Utility.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"
void SetUpArgumentParser(seqan::ArgumentParser &);
void SetDatabaseFastas(SeqIgOptions const &options, DatabasePaths &);
seqan::ArgumentParser::ParseResult ExtractOptions(seqan::ArgumentParser const &, SeqIgOptions &);
#endif
