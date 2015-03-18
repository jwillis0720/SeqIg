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
    _VGeneQueryStart(VGene.GetBeginQueryMatch()),
    _VGeneQueryEnd(VGene.GetEndQueryMatch()),

    //JGene Members Assignment From AlignAntibody Class
    _JGeneAlignment(JGene.GetTopAlignment()),
    _JGeneScore(JGene.GetTopScore()),
    _JGeneGene(JGene.GetTopGene()),
    _JGeneGeneStart(JGene.GetBeginGeneMatch()),
    _JGeneGeneEnd(JGene.GetEndGeneMatch()),
    _JGeneQueryStart(JGene.GetBeginQueryMatch()),
    _JGeneQueryEnd(JGene.GetEndQueryMatch()),

    _raw_sequence(raw_sequence),
	_vproperties(vproperties),
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
    _VGeneQueryStart(VGene.GetBeginQueryMatch()),
    _VGeneQueryEnd(VGene.GetEndQueryMatch()),

    //DGene Members Assignment From AlignAntibody Class
    _DGeneAlignment(DGene.GetTopAlignment()),
    _DGeneScore(DGene.GetTopScore()),
    _DGeneGene(DGene.GetTopGene()),
    _DGeneGeneStart(DGene.GetBeginGeneMatch()),
    _DGeneGeneEnd(DGene.GetEndGeneMatch()),
    _DGeneQueryStart(DGene.GetBeginQueryMatch()),
    _DGeneQueryEnd(DGene.GetEndQueryMatch()),
    
    //JGene Members Assignment From AlignAntibody Class
    _JGeneAlignment(JGene.GetTopAlignment()),
    _JGeneScore(JGene.GetTopScore()),
    _JGeneGene(JGene.GetTopGene()),
    _JGeneGeneStart(JGene.GetBeginGeneMatch()),
    _JGeneGeneEnd(JGene.GetEndGeneMatch()),
    _JGeneQueryStart(JGene.GetBeginQueryMatch()),
    _JGeneQueryEnd(JGene.GetEndQueryMatch()),

	_raw_sequence(raw_sequence),
	_vproperties(vproperties),
    _verbose(verbose)


{
    _setVGeneQueryStartTranslation();
    _setJunctions();

};

void AntibodyJunction::_setVGeneQueryStartTranslation() {

	for(int i = _VGeneGeneStart, j = 0; i < _VGeneGeneEnd ; i++, j++)
    {
        if(i % 3 == 0)
        {
        	_VGeneQueryStart += j;
            _VGeneGeneStart += j;
            break;
        }
    }

    if(_VGeneQueryStart > _JGeneQueryEnd){
        std::cerr << "J gene is becoming before V Gene \n \n";
        throw AntibodyJunctionException("J Gene coming before V Gene for Gene \n");

    }
    
    _EntireAntibodySeq = infix(_raw_sequence, _VGeneQueryStart, _JGeneQueryEnd);

    if(_verbose)
    {
    	std::cout << "V Gene Query Start " << _VGeneQueryStart << std::endl;
    	std::cout << "V Gene Start " << _VGeneGeneStart;
    	std::cout << "To J gene end " << _JGeneQueryEnd << std::endl;
    	std::cout << "Entire Antibody In Frame : " << _EntireAntibodySeq << std::endl;
    }

    if (seqan::translate(_AbAASeq,_EntireAntibodySeq) != 0)
        throw AntibodyJunctionException("Can't Translate Sequence");
    
    if(_verbose)std::cout << "\n\n" << _AbAASeq.concat << std::endl;
};

void AntibodyJunction::_setJunctions() {

	std::string start_block = "";

	int fw1_s = _vproperties[_VGeneGene]["FR1s"];
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

	int cdr1_length = cdr1_e - cdr1_s;
	int fw2_length = fw2_e - fw2_s;
	int cdr2_length = cdr2_e - cdr2_s;
	int fw3_length = fw3_e - fw3_s;


	//Do all the tests to figure out where the Query starts and assign junctions
	if(_VGeneGeneStart >= fw1_s && _VGeneGeneStart <= fw1_e){
		start_block = "FW1";
		int q_fw1_start = _VGeneQueryStart;
		int q_fw1_end = _VGeneQueryStart + (fw1_e - _VGeneGeneStart) + 1;
		int q_cdr1_start = q_fw1_end;
		int q_cdr1_end = q_fw1_end + cdr1_length + 1;
		int q_fw2_start = q_cdr1_end;
		int q_fw2_end = q_fw2_start + fw2_length + 1;
		int q_cdr2_start = q_fw2_end;
		int q_cdr2_end = q_cdr2_start + cdr2_length + 1;
		int q_fw3_start = q_cdr2_end;
		int q_fw3_end = q_fw3_start + fw3_length + 1;
		_setJunctionLenghts("FW1", q_fw1_start, q_fw1_end);
		_setJunctionLenghts("CDR1", q_cdr1_start, q_cdr1_end);
		_setJunctionLenghts("FW2", q_fw2_start, q_fw2_end);
		_setJunctionLenghts("CDR2", q_cdr2_start, q_cdr2_end);
		_setJunctionLenghts("FW3", q_fw3_start, q_fw3_end);
	}

	else if(_VGeneGeneStart >= cdr1_s && _VGeneGeneStart <= cdr1_e){
		start_block = "CDR1";
		int q_cdr1_start = _VGeneQueryStart;
		int q_cdr1_end = _VGeneQueryStart + (cdr1_e - _VGeneGeneStart);
		int q_fw2_start = q_cdr1_end;
		int q_fw2_end = q_fw2_start + fw2_length + 1;
		int q_cdr2_start = q_fw2_end;
		int q_cdr2_end = q_cdr2_start + cdr2_length + 1;
		int q_fw3_start = q_cdr2_end;
		int q_fw3_end = q_fw3_start + fw3_length + 1;
		_setJunctionLenghts("CDR1", q_cdr1_start, q_cdr1_end);
		_setJunctionLenghts("FW2", q_fw2_start, q_fw2_end);
		_setJunctionLenghts("CDR2", q_cdr2_start, q_cdr2_end);
		_setJunctionLenghts("FW3", q_fw3_start, q_fw3_end);
	}

	else if(_VGeneGeneStart >= fw2_s && _VGeneGeneStart <= fw2_e){
		start_block = "FW2";
		int q_fw2_start = _VGeneQueryStart;
		int q_fw2_end = _VGeneQueryStart + (fw2_e - _VGeneGeneStart);
		int q_cdr2_start = q_fw2_end;
		int q_cdr2_end = q_cdr2_start + cdr2_length + 1;
		int q_fw3_start = q_cdr2_end;
		int q_fw3_end = q_fw3_start + fw3_length + 1;
		_setJunctionLenghts("FW2", q_fw2_start, q_fw2_end);
		_setJunctionLenghts("CDR2", q_cdr2_start, q_cdr2_end);
		_setJunctionLenghts("FW3", q_fw3_start, q_fw3_end);
	}


	else if(_VGeneGeneStart >= cdr2_s && _VGeneGeneStart <= cdr2_e){
		start_block = "CDR2";
		int q_cdr2_start = _VGeneQueryStart;
		int q_cdr2_end = _VGeneQueryStart + (cdr2_e - _VGeneGeneStart);
		int q_fw3_start = q_cdr2_end;
		int q_fw3_end = q_fw3_start + fw3_length + 1;
		_setJunctionLenghts("CDR2", q_cdr2_start, q_cdr2_end);
		_setJunctionLenghts("FW3", q_fw3_start, q_fw3_end);
	}

	else if(_VGeneGeneStart >= fw3_s && _VGeneGeneStart <= fw3_e){
		start_block = "FW3";
		int q_fw3_start = _VGeneQueryStart;
		int q_fw3_end = _VGeneQueryStart + (fw3_e - _VGeneGeneStart);
		_setJunctionLenghts("FW3", q_fw3_start, q_fw3_end);
	}

	else if(_VGeneGeneStart >= cdr3_s){
		start_block = "CDR3";
		int q_cdr3_start = _VGeneQueryStart;
		int q_vh_end = _VGeneQueryEnd;
		_setJunctionLenghts("CDR3_v", q_cdr3_start, q_vh_end);
	}
	else start_block = "Unresolved";

	std::cout << start_block << std::endl;
	TJunctionsSE::iterator itrSE;
	TJunctionsNuc::iterator itrNuc;
	TJunctionsAA::iterator itrAA;
	for(itrSE = _JunctionStart_Ends.begin(), itrNuc = _JunctionsNuc.begin(), itrAA = _JunctionsAA.begin();
			itrSE != _JunctionStart_Ends.end(); itrSE++, itrNuc++, itrAA++)
	{
		std::cout << "\n" << itrSE->first << "\tstart: " <<
				itrSE->second.first << "\tend: " << itrSE->second.second;

		std::cout << "\nNuc Seq:\t" << itrNuc->second;
		std::cout << "\nAA Seq:\t" << itrAA->second << std::endl;
	}

}

void AntibodyJunction::_setJunctionLenghts(std::string const & region, int const & start, int const & end){
	_JunctionStart_Ends[region] = std::make_pair(start,end);
	Tds Nuc = infix(_raw_sequence, start, end);
	TAASeq AA;
	seqan::translate(AA, Nuc);
	_JunctionsNuc[region] = Nuc;
	_JunctionsAA[region] = AA.concat;
}


