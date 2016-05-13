/***********
	GA.cpp

	‹K–ñi2016/05/09j
		- std::vector< std::vector<int> >‚Ì•Ï”–¼‚Í~Population
************/

#include "../include/Ga.h"
#include "../include/GaCommon.h"


GA::GA(int numVariable)
{
	std::cout << "Calls constructor." << std::endl;

	int geneLength		= 20;	// ˆâ“`Žq’·
	int population		= 120;	// ŒÂ‘Ì”

	this->_population	= population;
	this->_geneLength	= geneLength;
	this->_numVariable	= numVariable;
}

GA::~GA()
{
}

/*
	’Tõ•êW’c‚ð‰Šú‰»‚·‚éD
	@param &searchPopulation
*/
void GA::_initSearchPopulation(std::vector<std::vector<int> > &searchPopulation)
{
	while (searchPopulation.size() < this->_population)
	{
		std::vector<int> tmpGene(this->_geneLength*this->_numVariable);
		this->_createRandomlyIndividual(tmpGene);

		if (!this->_isDuplicatedGene(searchPopulation, tmpGene))
			searchPopulation.push_back(tmpGene);
	}
}

/*
	ƒ‰ƒ“ƒ_ƒ€‚ÉŒÂ‘Ì‚ð1‚Â¶¬‚·‚é
	@param &individual ¶¬‚³‚ê‚éŒÂ‘Ìi—v‘f‚Ì—Ìˆæ‚ÍŠm•Û‚µ‚Ä‚¨‚­j
*/
void GA::_createRandomlyIndividual(std::vector<int> &individual)
{
	// [0.0, 1.0]‚Ìƒ‰ƒ“ƒ_ƒ€’l‚ðì¬
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<> randomValue(0.0, 1.0);

	for (int numBinary = 0; numBinary < this->_geneLength*this->_numVariable; ++numBinary)
		individual[numBinary]	= (randomValue(mt) > 0.5) ? 1 : 0;
}

/*
	1ŒÂ‘Ì‚ð•\Ž¦‚·‚é
	@param &individual ŒÂ‘Ì
*/
void GA::_outputGene(const std::vector<int> &individual)
{
	int numBinary;

	for (numBinary = 0; numBinary < this->_geneLength*this->_numVariable; ++numBinary)
		std::cout << individual[numBinary];

	std::cout << std::endl;
}

/*
	ŒÂ‘ÌW’c‚ð•\Ž¦‚·‚éD
	@param &targetPopulation ŒÂ‘ÌŒQ
*/
void GA::_outputPopulation(const std::vector<std::vector<int>> &targetPopulation)
{
	int numIndividual;

	for (numIndividual = 0; numIndividual < targetPopulation.size(); ++numIndividual)
		this->_outputGene(targetPopulation[numIndividual]);
}

/*
	ŒÂ‘ÌŒQ‚Ì’†‚É“¯ŒÂ‘Ì‚ª‘¶Ý‚µ‚Ä‚¢‚é‚©ƒ`ƒFƒbƒN
	@param &searchPopulation ŒÂ‘ÌW’c‚Ì2ŽŸŒ³”z—ñ
	@param &targetGene ŒŸõŒÂ‘Ì
*/
bool GA::_isDuplicatedGene(std::vector<std::vector<int>> &searchPopulation, const std::vector<int> &targetGene)
{
	// ’TõŒÂ‘Ì‚ª‘¶Ý‚µ‚È‚¢ê‡
	if (searchPopulation.size() == 0)
		return false;

	int numGene;
	for (numGene = 0; numGene < searchPopulation.size(); ++numGene)
		if (searchPopulation[numGene] == targetGene)
			return true;

	return false;
}

/*
	1ŒÂ‘Ì‚Ì2i”ƒf[ƒ^‚ð•\Œ»Œ^‚É•ÏŠ·
	@param &binary 1ŒÂ‘Ì‚Ìˆâ“`Žq
	@param &phenotype ˆâ“`ŽqŒ^‚ð•\Œ»Œ^‚É•ÏŠ·‚µ‚½‚à‚Ì
*/
void GA::_convertPhenotype(const std::vector<int> &binary, std::vector<double> &phenotype)
{
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	for (int numVar = 0; numVar < this->_numVariable; ++numVar)
	{
		for (int numBinary = 0; numBinary < this->_geneLength; ++numBinary)
			tmpGene.push_back(binary[numVar*this->_geneLength + numBinary]);

		phenotype.push_back(this->_binary2Phenotype(tmpGene));
		tmpGene.clear();
	}
}

/*
	1ˆâ“`Žq‚Ì2i”ƒf[ƒ^‚ð•\Œ»Œ^‚É•ÏŠ·
	@param &binary 1ˆâ“`Žq’·‚Ì’·‚³‚ðŽ‚Â2i”ƒf[ƒ^
*/
double GA::_binary2Phenotype(const std::vector<int> &binary)
{
	int numGene, place	= 0;
	double decimal	= 0.;

	for (numGene = this->_geneLength - 1; numGene >= 0 ; --numGene)
	{
		if (binary[numGene] == 1)
		{
			decimal	+= pow(2.0, (double)place);
			place	+= 1;
		}
	}

	return decimal / (pow(2.0, (double)this->_geneLength) - 1.0);
}

/*
	–Ú“IŠÖ”‚ðŒvŽZ
	@param &variable •Ï”
	@param &objectiveValue –Ú“IŠÖ”‚Ì’li–Ú“IŠÖ”‚Ì”‚¾‚¯—Ìˆæ‚ðŠm•Û‚µ‚Ä‚¨‚­j
*/
void GA::_getObjectiveFunc(const std::vector<double> &variable, std::vector<double> &objectiveValue)
{
	objectiveValue.at(0)	= this->_f1(variable);
	objectiveValue.at(1)	= this->_f2(variable);
}

double GA::_f1(const std::vector<double> &variable)
{
	return 1. - exp(-4.*variable[0])*pow(sin(6.*M_PI*variable[0]),6);
}

double GA::_f2(const std::vector<double> &var)
{
	double g, h;

	g = 1. + 9.*pow((var[1]+var[2]+var[3]+var[4]+var[5]+var[6]+var[7]+var[8]+var[9])/(9.),0.25);
	h = 1. - pow(this->_f1(var)/g,2);

	return g*h;
}

/*
	1ŒÂ‘Ì‚©‚ç–Ú“IŠÖ”’l‚ðŒvŽZ
	@param &binary 1ŒÂ‘Ì‚Ìˆâ“`Žq
	@param &obj –Ú“IŠÖ”‚Ì’li–Ú“IŠÖ”‚Ì”‚¾‚¯—Ìˆæ‚ðŠm•Û‚µ‚Ä‚¨‚­j
*/
void GA::_binary2ObjectiveFunc(const std::vector<int> &binary, std::vector<double> &obj)
{
	// •\Œ»Œ^‚ðˆêŽž“I‚É•Û‘¶
	std::vector<double> tmpPhenotype(this->_numVariable);

	this->_convertPhenotype(binary, tmpPhenotype);
	this->_getObjectiveFunc(tmpPhenotype, obj);
}

/*
	NSGA2ŽÀs—pƒƒ\ƒbƒh
*/
void GA::nsga2Run()
{
	int generation = 0;
	int maxGeneration	= 500;
	std::vector<std::vector<int> > margedPopulation, nextRankPopulation, crowdingSortedPopulation;
	std::vector<std::vector<std::vector<int> > > archivePopulation(maxGeneration), searchPopulation(maxGeneration), classifiedByRankGene;

	/*** Step1 ***/
	this->_initSearchPopulation(searchPopulation[generation]);

	while(1)
	{
		/*** Step2 ***/
		// TODO:•]‰¿•û–@‚ÌŠm—§
//		this->_outputObjectiveValue(searchPopulation[generation], generation);

		/*** Step3 ***/
		GaCommon::joinPopulation(archivePopulation[generation], searchPopulation[generation], margedPopulation);
		this->_nonSuperioritySort(margedPopulation, classifiedByRankGene);
		margedPopulation.clear();

		std::cout << classifiedByRankGene[0].size() << std::endl;
		return;

		/*** Step4 ***/
		this->_updateArchivePopulation(classifiedByRankGene, archivePopulation[generation+1], nextRankPopulation);

		/*** Step5 ***/
		this->_crowdingSort(nextRankPopulation, crowdingSortedPopulation);
		this->_insertIndividuals(archivePopulation[generation+1], nextRankPopulation);
		nextRankPopulation.clear();
		crowdingSortedPopulation.clear();

		/*** Step6 ***/
		if (generation < maxGeneration)
			break;

		/*** Step7 ***/
		// ‘I‘ð‚ÆŒð³‚ð“¯Žž‚És‚Á‚Ä‚¢‚é
		this->_crowdedTournamentSelection(archivePopulation[generation+1], searchPopulation[generation+1], classifiedByRankGene);
		classifiedByRankGene.clear();

		/*** Step8 ***/
		this->_mutationGene(searchPopulation[generation+1], 1/this->_geneLength*this->_numVariable);

		generation	+= 1;
	}
}

/*
	•êW’c‚É‘Î‚µ‚Ä”ñ—D‰zƒ\[ƒg‚ðs‚¤
	@param gene ƒ\[ƒg‚ðs‚¤•êW’c
	@param &classifiedByRankGene ƒ‰ƒ“ƒN‚²‚Æ‚ÉƒNƒ‰ƒX•ª‚¯‚µ‚½ŒÂ‘ÌW’c 
*/
void GA::_nonSuperioritySort(
	const std::vector <std::vector<int> > &targetPopulation,
	std::vector<std::vector<std::vector<int> > > &classifiedByRankGene)
{
	int tmp, numGene;
	std::vector<std::vector<int> > sortingPopulation, tmpRankedGene;

	// ŒÂ‘Ì‚Éƒ‰ƒ“ƒN‚ð•t‚¯‚Ä‚¢‚«Cƒ‰ƒ“ƒN•t‚¯‚³‚ê‚½ŒÂ‘Ì‚Íœ‚­
	std::copy(targetPopulation.begin(), targetPopulation.end(), std::back_inserter(sortingPopulation));
	while(sortingPopulation.size() > 0)
	{
		for (numGene = 0; numGene < sortingPopulation.size(); ++numGene)
			if (this->_isSuperior(sortingPopulation[numGene], sortingPopulation))
				tmpRankedGene.push_back(sortingPopulation[numGene]);

		// ƒ‰ƒ“ƒN‚²‚Æ‚ÉŒÂ‘Ì‚ð•Û‘¶‚µCŒÂ‘ÌŒQ‚ðXV
		classifiedByRankGene.push_back(tmpRankedGene);
		for (tmp = 0; tmp < tmpRankedGene.size(); ++tmp)
			GaCommonTemp<int>::removeElement(sortingPopulation, tmpRankedGene[tmp]);
		tmpRankedGene.clear();
	}
}

/*
	Žw’è‚µ‚½ŒÂ‘Ì‚ª—D‰z‚µ‚Ä‚¢‚é‚©”»’è
	@param &targetGene —D‰z‚µ‚Ä‚¢‚é‚©”»’è‚µ‚½‚¢ŒÂ‘Ì
	@param &comparedPopulation targetGene‚ª‘®‚µ‚Ä‚¢‚éŒÂ‘ÌŒQ
*/
bool GA::_isSuperior(
	const std::vector<int> &targetGene,
	const std::vector<std::vector<int> > &comparedPopulation)
{
	int numGene, numObj;
	bool isSuperior;
	std::vector<double> targetObjeciveValue(2), comparedObjectiveValue(2);	// –Ú“IŠÖ”‚Ì”

	isSuperior	= true;
	this->_binary2ObjectiveFunc(targetGene, targetObjeciveValue);

	for (numGene = 0; numGene < comparedPopulation.size(); ++numGene)
	{
		if (comparedPopulation[numGene] == targetGene)
			continue;

		this->_binary2ObjectiveFunc(comparedPopulation[numGene], comparedObjectiveValue);
		for (numObj = 0; numObj < 2; ++numObj)	// –Ú“IŠÖ”‚Ì”
		{
			if (targetObjeciveValue[numObj] < comparedObjectiveValue[numObj])
				isSuperior	= false;
		}
	}

	return isSuperior;
}

/*
	V‚½‚ÈƒA[ƒJƒCƒuW’c‚ð¶¬
	@param &classifiedByRankGene ƒ‰ƒ“ƒN‚²‚Æ‚ÌŒÂ‘Ì
	@param &newArchivePopulation ƒ‰ƒ“ƒNãˆÊ‚©‚çŽæ“¾‚µ‚½XV—pƒA[ƒJƒCƒu•êW’c
	@param &nextRankPopulation ƒA[ƒJƒCƒu•êW’c‚É“ü‚è‚«‚ç‚È‚©‚Á‚½Å‚ƒ‰ƒ“ƒN‚ÌŒÂ‘ÌW’c
*/
void GA::_updateArchivePopulation(
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene,
	std::vector<std::vector<int> > &newArchivePopulation,
	std::vector<std::vector<int> > &nextRankPopulation)
{
	int rank;

	std::cout << classifiedByRankGene.size() << std::endl;
	for (rank = 0; rank < classifiedByRankGene.size(); ++rank)
	{
		if (this->_population - newArchivePopulation.size() > classifiedByRankGene[rank].size())
			GaCommon::pushBackAll(newArchivePopulation, classifiedByRankGene[rank]);
		else
			break;
	}

	// ƒA[ƒJƒCƒu•êW’c‚É“ü‚è‚«‚ç‚È‚©‚Á‚½Å‚ƒ‰ƒ“ƒN‚ÌŒÂ‘ÌW’c‚ðŠi”[
	std::copy(classifiedByRankGene[rank].begin(), classifiedByRankGene[rank].end(), std::back_inserter(nextRankPopulation));
}

/*
	¬ŽG“xƒ\[ƒg
	@param &certainRankPopulation ‚ ‚éƒ‰ƒ“ƒN‚ÌŒÂ‘ÌŒQ
	@param &crowdingSortedPopulation ¬ŽG“x‚²‚Æ‚Éƒ\[ƒg‚³‚ê‚½ŒÂ‘ÌŒQ
*/
void GA::_crowdingSort(
	const std::vector<std::vector<int > > &certainRankPopulation,
	std::vector<std::vector<int> > &crowdingSortedPopulation)
{
	// ƒJƒEƒ“ƒg•Ï”
	int tmp1, tmp2;

	int numObj;
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);
	std::vector<double> tmpObjectFunc(2);	// –Ú“IŠÖ”‚Ì”
	std::vector<std::vector<std::vector<int> > > objectiveSortedGene(2);	// –Ú“IŠÖ”‚Ì”	

	// –Ú“IŠÖ”‚²‚Æ‚É–Ú“IŠÖ”’l‚ªˆ«‚¢‡‚É•À‚×‚é
	this->_putObjectiveSortedGeneEveryObjectiveFunc(certainRankPopulation, objectiveSortedGene);

	// ¬ŽG“x‚ª‘å‚«‚¢‡‚ÉŒÂ‘Ì‚ðƒ\[ƒg
	std::copy(certainRankPopulation.begin(), certainRankPopulation.end(), std::back_inserter(crowdingSortedPopulation));
	double distance1 = 0., distance2 = 0.;
	for (tmp1 = 1; tmp1 < certainRankPopulation.size()-1; ++tmp1)
	{
		for (tmp2 = certainRankPopulation.size()-2;  tmp2 > tmp1; --tmp2)
		{
			distance1	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, crowdingSortedPopulation[tmp2]);
			distance2	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, crowdingSortedPopulation[tmp2-1]);
			if (distance1 > distance2)
			{
				tmpGene	= crowdingSortedPopulation[tmp2];
				crowdingSortedPopulation[tmp2]	= crowdingSortedPopulation[tmp2-1];
				crowdingSortedPopulation[tmp2-1]	= tmpGene;
			}
		}
	}
}

/*
	–Ú“IŠÖ”‚²‚Æ‚É–Ú“IŠÖ”’l‚ªˆ«‚¢‡‚É•À‚×‚é
	@param &targetPopulation ƒ\[ƒg‚µ‚½‚¢ŒÂ‘ÌŒQ
	@param &objectiveSortedGene –Ú“IŠÖ”‚²‚Æ‚Éƒ\[ƒg‚³‚ê‚½ŒÂ‘Ìi–Ú“IŠÖ”‚Ì”‚¾‚¯—Ìˆæ‚ðŠm•Û‚µ‚Ä‚¨‚­j
*/
void GA::_putObjectiveSortedGeneEveryObjectiveFunc(
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<std::vector<std::vector<int> > > &objectiveSortedGene)
{
	int numObj;
	for (numObj = 0; numObj < 2; ++numObj)	// –Ú“IŠÖ”‚Ì”
	{
		this->_sortByObjectiveValue(targetPopulation, objectiveSortedGene[numObj], numObj);
		std::reverse(objectiveSortedGene[numObj].begin(), objectiveSortedGene[numObj].end());
	}
}

/*
	Žw’è‚µ‚½–Ú“IŠÖ”’l‚ª¬‚³‚¢‡‚ÉŒÂ‘Ì‚ðƒoƒuƒ‹ƒ\[ƒg‚·‚é
	@param &targetGene ƒ\[ƒg‚µ‚½‚¢ŒÂ‘ÌŒQ
	@param &objectiveSortedGene ƒ\[ƒgŒã‚ÌŒÂ‘ÌŒQ
	@param num ‘ÎÛ‚Æ‚·‚é–Ú“IŠÖ”‚Ì”Ô†
*/
void GA::_sortByObjectiveValue(
	const std::vector<std::vector<int> > &targetGene,
	std::vector<std::vector<int> > &objectiveSortedPopulation,
	int num)
{
	// ƒJƒEƒ“ƒg•Ï”
	int tmp1, tmp2;

	std::vector<double> tmpObject1(2);	// –Ú“IŠÖ”‚Ì”
	std::vector<double> tmpObject2(2);	// –Ú“IŠÖ”‚Ì”
	std::vector<int> tmpGene(this->_geneLength*this->_numVariable);

	std::copy(targetGene.begin(), targetGene.end(), std::back_inserter(objectiveSortedPopulation));

	for (tmp1 = 0; tmp1 < targetGene.size(); ++tmp1)
	{
		for (tmp2 = targetGene.size() - 1; tmp2 > tmp1; --tmp2)
		{
			this->_binary2ObjectiveFunc(targetGene[tmp2], tmpObject1);
			this->_binary2ObjectiveFunc(targetGene[tmp2-1], tmpObject2);
			if (tmpObject1[num] < tmpObject2[num])
			{
				tmpGene	= objectiveSortedPopulation[tmp2];
				objectiveSortedPopulation[tmp2]	= objectiveSortedPopulation[tmp2-1];
				objectiveSortedPopulation[tmp2-1]	= tmpGene;
			}
		}
	}
}

/*
	Žw’è‚µ‚½–Ú“IŠÖ”‚É‘Î‚µ‚Ä¬ŽG“x‚ðŒvŽZ‚·‚é
	@param &objectiveSortedGene –Ú“IŠÖ”‚²‚Æ‚Éƒ\[ƒg‚³‚ê‚½ŒÂ‘ÌŒQ
	@param numGene ŒÂ‘Ì‚Ì”Ô†i‹«ŠEŒÂ‘Ì‚Í‘I‘ð‚Å‚«‚È‚¢j
	@param numObj –Ú“IŠÖ”‚Ì”Ô†
*/
double GA::_culcCrowdingDistanse(
	const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene,
	int numGene,
	int numObj)
{
	std::vector<double> tmpObjLeft(2);		// –Ú“IŠÖ”‚Ì”
	std::vector<double> tmpObjRight(2);	// –Ú“IŠÖ”‚Ì”
	std::vector<double> tmpObjMax(2);		// –Ú“IŠÖ”‚Ì”
	std::vector<double> tmpObjMin(2);		// –Ú“IŠÖ”‚Ì”

	// ‹«ŠEŒÂ‘Ì‚ðœ‚¢‚ÄŒÂ‘Ì‚ð‘I‘ð
	double distance	= 0.;
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][numGene-1], tmpObjLeft);
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][numGene+1], tmpObjRight);
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][0], tmpObjMax);
	this->_binary2ObjectiveFunc(objectiveSortedGene[numObj][objectiveSortedGene.front().size()-1], tmpObjMin);
	distance	= (tmpObjLeft[numObj] - tmpObjRight[numObj])/(tmpObjMax[numObj] - tmpObjMin[numObj]);

	return distance;
}

/*
	Žw’è‚µ‚½ŒÂ‘Ì‚Ì‘¬ŽG“x‚ðŒvŽZ‚·‚é
	@param &objectiveSortedGene –Ú“IŠÖ”‚²‚Æ‚É–Ú“IŠÖ”’n‚ªˆ«‚¢‡‚Éƒ\[ƒg‚³‚ê‚½ŒÂ‘ÌŒQ
	@param &individual ŒÂ‘Ì‚Ìˆâ“`Žq
*/
double GA::_culcCrowdingDistanseForIndividual(
	const std::vector<std::vector<std::vector<int> > > &objectiveSortedGene,
	const std::vector<int> &individual)
{
	int numObj, numGene;
	std::vector<double> distance(2);	// –Ú“IŠÖ”‚Ì”

	for (numObj = 0; numObj < 2; ++numObj)	// –Ú“IŠÖ”‚Ì”
	{
		for (numGene = 1; numGene < objectiveSortedGene.front().size()-1; ++numGene)
		{
			if (objectiveSortedGene[numObj][numGene] == individual)
			{
				distance[numObj]	= this->_culcCrowdingDistanse(objectiveSortedGene, numGene, numObj);
				break;
			}
		}
	}

	double resultDistanse	= 0.;
	for (numObj = 0; numObj < distance.size(); ++numObj)
	{
		resultDistanse	+= distance[numObj];
	}

	return  resultDistanse;
}

/*
	ŒÂ‘Ì”‚ªN‚É‚È‚é‚Ü‚ÅŒÂ‘Ì‚ð’Ç‰Á‚·‚é
	@param &insertedPopulation ŒÂ‘Ì‚ð’Ç‰Á‚³‚ê‚éŒÂ‘ÌW’ciŒÂ‘Ì”‚ÍNˆÈ‰ºj
	@param &insertPopulation ’Ç‰Á‚·‚éŒÂ‘Ì‚ðŠÜ‚ÞŒÂ‘ÌW’c
*/
void GA::_insertIndividuals(
	std::vector<std::vector<int> > &insertedPopulation,
	const std::vector<std::vector<int> > &insertPopulation)
{
	if (insertedPopulation.size() >= this->_population)
	{
		std::cout << "ERROR:Excess individual!!" << std::endl;
		return;
	}

	for (int tmp = 0; insertedPopulation.size() < this->_population; ++tmp)
		insertedPopulation.push_back(insertPopulation[tmp]);
}

/*
	¬ŽG“xƒg[ƒiƒƒ“ƒg‘I‘ð
	@param &selectedPopulation ‘I‘ð‚³‚ê‚éŒÂ‘ÌW’c
	@param &newSearchPopulation V‚½‚È’Tõ•êW’c
	@param &classifiedByRankGene ƒ‰ƒ“ƒN‚²‚Æ‚ÉƒNƒ‰ƒX•ª‚¯‚³‚ê‚½ŒÂ‘Ì
*/
void GA::_crowdedTournamentSelection(
	const std::vector<std::vector<int> > &selectedPopulation,
	std::vector<std::vector<int> > &newSearchPopulation,
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene)
{
	std::vector<int> parentGene1, parentGene2, childGene1, childGene2;
	std::vector<std::vector<int> > tmpSelectionPopulation, highRankPopulation;

	// ŒÂ‘Ì”‚ªN‚É‚È‚é‚Ü‚Å‘I‘ð‚ðŽÀs
	// eŒÂ‘Ì‚ðƒ‰ƒ“ƒ_ƒ€‚É‘I‘ð‚µCˆê—lŒð³‚ðŽÀs
	// e~2‚ÆŽq~2‚Ì‚¤‚¿ƒ‰ƒ“ƒN‚ªãˆÊ‚ÌŒÂ‘Ì‚ð1‚Â‘I‘ð‚µCV‚½‚È’Tõ•êW’c‚É’Ç‰Á‚·‚é
	for (int numGene = 0; newSearchPopulation.size() < this->_population; ++numGene)
	{
		this->_select2GenesFromPopulation(selectedPopulation, parentGene1, parentGene2);
		this->_uniformCrossover(parentGene1, parentGene2, childGene1, childGene2);

		tmpSelectionPopulation.push_back(parentGene1);
		tmpSelectionPopulation.push_back(parentGene2);
		tmpSelectionPopulation.push_back(childGene1);
		tmpSelectionPopulation.push_back(childGene2);

		this->_highRankGeneSelection(classifiedByRankGene, tmpSelectionPopulation, highRankPopulation, 1);
		newSearchPopulation.push_back(highRankPopulation[1]);
	}
}

/*
	ŒÂ‘Ì‚Ìƒ‰ƒ“ƒN‚ð•Ô‚·
	@param &classifiedByRankGene ƒ‰ƒ“ƒN‚²‚Æ‚ÉƒNƒ‰ƒX•ª‚¯‚³‚ê‚½ŒÂ‘Ì
	@param &targetGene ƒ‰ƒ“ƒN‚ð’m‚è‚½‚¢ŒÂ‘Ì
*/
int GA::_returnGeneRank(
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene,
	const std::vector<int> &targetGene)
{
	int rank, numGene;

	for (rank = 0; rank < classifiedByRankGene.size(); ++rank)
	{
		for (numGene = 0; numGene < classifiedByRankGene[rank].size(); ++numGene)
		{
			if (classifiedByRankGene[rank][numGene] == targetGene)
				return rank;
		}
	}
}

/*
	ŒÂ‘ÌŒQ‚©‚çƒ‰ƒ“ƒ_ƒ€‚É2ŒÂ‘Ì‚ð‘I‘ð‚·‚é
	@param &targetPopulation ‘I‘ð‚·‚éŒÂ‘ÌŒQ
	@param &gene1 ‘I‘ð‚³‚ê‚½ŒÂ‘Ì1
	@param &gene2 ‘I‘ð‚³‚ê‚½ŒÂ‘Ì2
*/
void GA::_select2GenesFromPopulation(
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<int> &gene1,
	std::vector<int> &gene2)
{
	int geneNum1, geneNum2;

	// ŒÂ‘Ì‚ð‘I‘ð‚·‚é‚½‚ß‚Ìƒ‰ƒ“ƒ_ƒ€’l‚ð¶¬
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_int_distribution<int> randomValue(0, targetPopulation.size()-1);

	// ƒ‰ƒ“ƒ_ƒ€‚ÉŒÂ‘Ì‚ð‘I‘ð
	do
	{
		geneNum1	= randomValue(mt);
		geneNum2	= randomValue(mt);
	} while (geneNum1 == geneNum2);

	gene1	= targetPopulation[geneNum1];
	gene2	= targetPopulation[geneNum2];
}

/*
	ˆê—lŒð³‚ð‚¨‚±‚È‚¤
	@param &parentGene1 eŒÂ‘Ì1
	@param &parentGene2 eŒÂ‘Ì2
	@param &childGene1 ŽqŒÂ‘Ì1
	@param &childGene2 ŽqŒÂ‘Ì2
*/
void GA::_uniformCrossover(
	const std::vector<int> &parentGene1,
	const std::vector<int> &parentGene2,
	std::vector<int> &childGene1,
	std::vector<int> &childGene2)
{
	// ƒJƒEƒ“ƒg•Ï”
	int tmp;

	// ƒ}ƒXƒNƒpƒ^[ƒ“—p‚Ìƒ‰ƒ“ƒ_ƒ€’l¶¬Ší
	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_int_distribution<int> randomValue(0, 1);

	// ƒ}ƒXƒNƒpƒ^[ƒ“‚ð¶¬
	std::vector<int> maskPattern(this->_geneLength*this->_numVariable);
	for (tmp = 0; tmp < this->_geneLength; tmp++)
		maskPattern.push_back(randomValue(mt));

	// Œð³
	for (tmp = 0; tmp < this->_geneLength*this->_numVariable; tmp++)
	{
		if (maskPattern[tmp] == 0)
		{
			childGene1.push_back(parentGene1[tmp]);
			childGene2.push_back(parentGene2[tmp]);
		}
		else if (maskPattern[tmp] == 1)
		{
			childGene1.push_back(parentGene2[tmp]);
			childGene2.push_back(parentGene1[tmp]);
		}
	}
}

/*
	Žw’è‚µ‚½ŒÂ‘Ì‚©‚çãˆÊƒ‰ƒ“ƒNŒÂ‘Ì‚ð‘I‘ð
	@param &classifiedByRankGene ƒ‰ƒ“ƒN‚²‚Æ‚ÉƒNƒ‰ƒX•ª‚¯‚³‚ê‚½ŒÂ‘Ì
	@param &targetPopulation ‘I‘ð‚·‚éŒÂ‘ÌŒQ
	@param &highRankPopulation ãˆÊƒ‰ƒ“ƒN‚ÌŒÂ‘ÌŒQ
	@param num ‘I‘ð‚·‚éŒÂ‘Ì”
*/
void GA::_highRankGeneSelection(
	const std::vector<std::vector<std::vector<int> > > &classifiedByRankGene,
	const std::vector<std::vector<int> > &targetPopulation,
	std::vector<std::vector<int> > &highRankPopulation,
	int num)
{
	int geneRank;
	double geneDistance, longestDistance = 0.;
	std::vector<int> tmpHighRankGene(this->_geneLength*this->_numVariable);
	std::vector<std::vector<std::vector<int> > > objectiveSortedGene(2);	// –Ú“IŠÖ”‚Ì”	

	// ‘ÎÛ‚ÌŒÂ‘ÌŒQ‚ðƒ‰ƒ“ƒN‚²‚Æ‚É•ª‚¯‚é
	std::vector<std::vector<std::vector<int> > > tmpClassifiedByRankGene;
	this->_nonSuperioritySort(targetPopulation, tmpClassifiedByRankGene);

	// ãˆÊƒ‰ƒ“ƒN‚ÌŒÂ‘Ì‚ðnumŒÂ‘I‘ð‚·‚é
	for (int rank = 0; highRankPopulation.size() < num; ++rank)
	{
		if (tmpClassifiedByRankGene[rank].size() == 1)
			highRankPopulation.push_back(tmpClassifiedByRankGene[rank][0]);
		else if (tmpClassifiedByRankGene.size() > 1)
		{
			// “¯ƒ‰ƒ“ƒN‚É•¡”ŒÂ‘Ì‚ª‘¶Ý‚µ‚½ê‡
			// Å‚à¬ŽG‹——£‚ª’·‚¢ŒÂ‘Ì‚ð‘I‘ð‚µ‚Ä’Ç‰Á‚·‚é
			for (int tmp = 0; tmp < tmpClassifiedByRankGene[rank].size(); ++tmp)
			{
				// ‘SŒÂ‘Ì’†‚Ìƒ‰ƒ“ƒN‚©‚ç¬ŽG“x‹——£‚ð‹‚ß‚é
				geneRank	= this->_returnGeneRank(classifiedByRankGene, tmpClassifiedByRankGene[rank][tmp]);
				this->_putObjectiveSortedGeneEveryObjectiveFunc(classifiedByRankGene[geneRank], objectiveSortedGene);
				geneDistance	= this->_culcCrowdingDistanseForIndividual(objectiveSortedGene, tmpClassifiedByRankGene[rank][tmp]);
				if (longestDistance < geneDistance)
				{
					longestDistance = geneDistance;
					tmpHighRankGene	= tmpClassifiedByRankGene[rank][tmp];
				}
			}

			highRankPopulation.push_back(tmpHighRankGene);
		}
	}
}

/*
	5”Ô–Ú‚ÌŒÂ‘Ì‚É“Ë‘R•ÏˆÙ‚ðs‚¤D
	“K‰ž“xŒvŽZ‘O‚ÌŒÂ‘Ì‚Ìˆâ“`Žq(0 or 1)‚ðƒ‰ƒ“ƒ_ƒ€‚É“ü‚ê‘Ö‚¦‚é
	@param &targetPopulation ‘ÎÛ‚ÌŒÂ‘ÌŒQ
	@param mutationRate “Ë‘R•ÏˆÙ—¦
*/
void GA::_mutationGene(
	std::vector<std::vector<int> > &targetPopulation,
	double mutationRate)
{
	int geneNum = 5;

	std::random_device seedGen;
	std::mt19937 mt(seedGen());
	std::uniform_real_distribution<double> randomValue(0.0, 1.0);

	for (int tmp = 0; tmp < this->_geneLength*this->_numVariable; tmp++)
	{
		if (mutationRate > randomValue(mt))
		{
			switch (targetPopulation[geneNum][tmp])
			{
			case 0:
				targetPopulation[geneNum][tmp] = 1;
			case 1:
				targetPopulation[geneNum][tmp] = 0;
			default:
				break;
			}
		}
	}
}

/*
	Žw’è‚µ‚½¢‘ã‚ÌŒÂ‘ÌW’c‚Ìx‚Æ“K‰ž“x‚ð•\Ž¦‚·‚é
	@param &targetPopulation ‘ÎÛ‚ÌŒÂ‘ÌŒS
	@param generation Œ»Žž“_‚Å‚Ì¢‘ã
*/
void GA::_outputObjectiveValue(
	std::vector<std::vector<int> > targetPopulation,
	int generation)
{
	double x;
	std::vector<std::vector<double> > objectiveParameter(targetPopulation.size(), std::vector<double>(2));
	std::vector<double> objectiveValue(2);

	std::cout << generation << "-generation" << std::endl;
	for (int numGene = 0; numGene < targetPopulation.size(); ++numGene)
	{
		this->_binary2ObjectiveFunc(targetPopulation[numGene], objectiveParameter[numGene]);
		this->_getObjectiveFunc(objectiveParameter[numGene], objectiveValue);

		for (int numVar = 0; numVar < this->_numVariable; ++numVar)
		{
			std::cout << objectiveParameter[numGene][numVar] << ",";
		}
		std::cout << "\t";
		for (int numObj = 0; numObj < 2; ++numObj)	// –Ú“IŠÖ”‚Ì”
		{
			std::cout << objectiveValue[numObj] << ",";
		}
	}
	std::cout << std::endl;
}