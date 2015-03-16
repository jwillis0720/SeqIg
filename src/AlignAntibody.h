//
//  AlignAntibody.h
//  seqan_sandbox
//
//  Created by Jordan Willis on 1/11/15.
//
//

#ifndef __AlignAntibody__
#define __AlignAntibody__

#include <stdio.h>
#include <seqan/sequence.h>
#include "DatabaseHandler.h"
#include <seqan/align.h>

class AlignAntibody
{

private:
    //primatives
    bool _verbose;
    int _top_score;
    
    //best alignment holders
    seqan::Dna5String _current_seq;
    seqan::CharString _current_id;
    seqan::CharString _best_gene;
    TAlignmnet _best_align;
    

    //private methods
    void DoPairWiseLocalAgainstDb(seqan::CharString const &,
                                  seqan::Dna5String const &);

public:
    AlignAntibody(seqan::CharString const &,
                      seqan::Dna5String const &,
                      Tdbcontainer const &,
                      bool const &);
    
    //Getters
    int GetTopScore() const;
    TAlignmnet GetTopAlignment() const;
    seqan::CharString GetTopGene() const;
    int GetBeginQueryMatch() const;
    int GetEndQueryMatch() const;
    int GetBeginGeneMatch() const;
    int GetEndGeneMatch() const;
    
    //For verbose
    void PrintBestAlignment() const;
};

#endif /* defined(__seqan_sandbox__AlignAntibody__) */
