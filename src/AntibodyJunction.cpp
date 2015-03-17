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
		TProperties const & vproperties,
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

	_vproperties(vproperties),
    _raw_sequence(raw_sequence),
    _verbose(verbose)

{

    _setVGeneQueryStartTranslation();

};

AntibodyJunction::AntibodyJunction(AlignAntibody const & VGene,
		AlignAntibody const & JGene,
		AlignAntibody const & DGene,
		Tds const & raw_sequence,
		TProperties const & vproperties,
		bool const & verbose):

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

	_vproperties(vproperties),
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
            _VGeneGeneStart = i;
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

	std::string start_block = "";

	int fw1_s = _vproperties["IGHV3-33*01"]["FR1s"];
	int fw1_e = _vproperties[_VGeneGene]["FR1e"];
	int cdr1_s = _vproperties[_VGeneGene]["CDR1s"];
	int cdr1_e = _vproperties[_VGeneGene]["CDR1e"];
	int fw2_s = _vproperties[_VGeneGene]["FR2s"];
	int fw2_e = _vproperties[_VGeneGene]["FR2e"];
	int cdr2_s = _vproperties[_VGeneGene]["CDR2s"];
	int cdr2_e = _vproperties[_VGeneGene]["CDR2e"];
	int fw3_s = _vproperties[_VGeneGene]["FR3s"];
	int fw3_e = _vproperties[_VGeneGene]["FR3e"];
	int cdr3_s = _vproperties[_VGeneGene]["CDR3s"];

	//Do all the tests to figure out where the Query starts
	if(_VGeneGeneStart >= fw1_s && _VGeneGeneStart <= fw1_e) start_block = "FW1";
	else if(_VGeneGeneStart >= cdr1_s && _VGeneGeneStart <= cdr1_e) start_block = "CDR1";
	else if(_VGeneGeneStart >= fw2_s && _VGeneGeneStart <= fw2_e) start_block = "FW2";
	else if(_VGeneGeneStart >= cdr2_s && _VGeneGeneStart <= cdr2_e) start_block = "CDR2";
	else if(_VGeneGeneStart >= fw3_s && _VGeneGeneStart <= fw3_e) start_block = "FW3";
	else if(_VGeneGeneStart >= cdr3_s) start_block = "CDR3";
	else start_block = "Unresolved";

	std::cout << "Starting V Match in " << start_block << std::endl;


	/*
	 * Todo the FW1, CDR1 should look like this
	 * FW1 = Raw_sequence[VGeneTranslationStart:(VGenequerystart+(fw1e - Vgenengenestart)]
	 * See board to validate that, its easy
	 */
}




