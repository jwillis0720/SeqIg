//
//  AntibodyJunction.cpp
//  seqan_sandbox
//
//  Created by Jordan Willis on 1/14/15.
//
//
#include "AntibodyJunction.h"
#include <seqan/translation.h>

//Constructor One if we have just V and J Gene
AntibodyJunction::AntibodyJunction(AlignAntibody const & VGene,
		AlignAntibody const & JGene,
		Tds const & raw_sequence,
		bool const & verbose):

    //VGene Members Assignment From AlignAntibody Class
    _VGeneAlignment(VGene.GetTopAlignment()),
    _VGeneScore(VGene.GetTopScore()),
    _VGeneGene(VGene.GetTopGene()),
    _VGeneGeneStart(VGene.GetBeginGeneMatch()),
    _VGeneGeneEnd(VGene.GetEndGeneMatch()),
    _VGeneQueryStart(VGene.GetBeginGeneMatch()),
    _VGeneQueryEnd(VGene.GetEndQueryMatch()),

    //JGene Members Assignment From AlignAntibody Class
    _JGeneAlignment(JGene.GetTopAlignment()),
    _JGeneScore(JGene.GetTopScore()),
    _JGeneGene(JGene.GetTopGene()),
    _JGeneGeneStart(JGene.GetBeginGeneMatch()),
    _JGeneGeneEnd(JGene.GetEndGeneMatch()),
    _JGeneQueryStart(JGene.GetBeginGeneMatch()),
    _JGeneQueryEnd(JGene.GetEndQueryMatch()),

    _raw_sequence(raw_sequence),
    _verbose(verbose)
{

    _setVGeneQueryStartTranslation();

};

AntibodyJunction::AntibodyJunction(AlignAntibody const & VGene, AlignAntibody const & JGene,AlignAntibody const & DGene, Tds const & raw_sequence, bool const & verbose):

    //VGene Members Assignment From AlignAntibody Class
    _VGeneAlignment(VGene.GetTopAlignment()),
    _VGeneScore(VGene.GetTopScore()),
    _VGeneGene(VGene.GetTopGene()),
    _VGeneGeneStart(VGene.GetBeginGeneMatch()),
    _VGeneGeneEnd(VGene.GetEndGeneMatch()),
    _VGeneQueryStart(VGene.GetBeginGeneMatch()),
    _VGeneQueryEnd(VGene.GetEndQueryMatch()),

    //DGene Members Assignment From AlignAntibody Class
    _DGeneAlignment(DGene.GetTopAlignment()),
    _DGeneScore(DGene.GetTopScore()),
    _DGeneGene(DGene.GetTopGene()),
    _DGeneGeneStart(DGene.GetBeginGeneMatch()),
    _DGeneGeneEnd(DGene.GetEndGeneMatch()),
    _DGeneQueryStart(DGene.GetBeginGeneMatch()),
    _DGeneQueryEnd(DGene.GetEndQueryMatch()),
    
    //JGene Members Assignment From AlignAntibody Class
    _JGeneAlignment(JGene.GetTopAlignment()),
    _JGeneScore(JGene.GetTopScore()),
    _JGeneGene(JGene.GetTopGene()),
    _JGeneGeneStart(JGene.GetBeginGeneMatch()),
    _JGeneGeneEnd(JGene.GetEndGeneMatch()),
    _JGeneQueryStart(JGene.GetBeginGeneMatch()),
    _JGeneQueryEnd(JGene.GetEndQueryMatch()),

    _raw_sequence(raw_sequence),
    _verbose(verbose)

{
    _setVGeneQueryStartTranslation();
    _setJunctions();

};

void AntibodyJunction::_setVGeneQueryStartTranslation() {

	for(int i = _VGeneGeneStart, j = _VGeneQueryStart; i < _VGeneGeneEnd ; i++, j++)
    {
        if(i % 3 == 0)
        {
            _VGeneQueryStartTranslation = j;
            break;
        }
    }

    
    if(_VGeneQueryStartTranslation > _JGeneQueryEnd){
        std::cerr << "J gene is becoming before V Gene \n \n";
        throw AntibodyJunctionException("J Gene coming before V Gene for Gene \n");

    }
    
    _EntireAntibodySeq = infix(_raw_sequence, _VGeneQueryStartTranslation, _JGeneQueryEnd);

    if(_verbose)
    {
    	std::cout << "V Gene Start " << _VGeneQueryStartTranslation << std::endl;;
    	std::cout << "Entire Antibody In Frame : " << _EntireAntibodySeq << std::endl;
    }

    if (seqan::translate(_AbAASeq,_EntireAntibodySeq) != 0)
        throw AntibodyJunctionException("Can't Translate Sequence");
    
    if(_verbose)std::cout << "\n\n" << _AbAASeq.concat << std::endl;
};

void AntibodyJunction::_setJunctions() {
std::cout << "blah" << std::endl;


}




