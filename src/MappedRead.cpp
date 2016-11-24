#include "MappedRead.h"

#include <stdlib.h>

#undef module_name
#define module_name "READ"

#ifdef INSTANCE_COUNTING
volatile int MappedRead::sInstanceCount = 0;
volatile int LocationScore::sInstanceCount = 0;
#endif

MappedRead::MappedRead(int const readid, int const qrymaxlen) :
		ReadId(readid), Calculated(-1), qryMaxLen(qrymaxlen), numTopScores(
				1), Scores(0), Alignments(0), iScores(0), Paired(0), Status(
				0), mappingQlty(255), s(0), length(0), polyATrimmed(0), RevSeq(0), Seq(0), qlty(0), name(0), AdditionalInfo(0) {
#ifdef INSTANCE_COUNTING
	AtomicInc(&sInstanceCount);
	maxSeqCount = std::max(sInstanceCount, maxSeqCount);
#endif
	//Name
	static size_t const MAX_READNAME_LENGTH = 100;
	name = new char[MAX_READNAME_LENGTH];

	Seq = new char[qryMaxLen];
	memset(Seq, '\0', qryMaxLen);

	qlty = new char[qryMaxLen];
	memset(qlty, '\0', qryMaxLen);
}

static inline char cpl(char c) {
	if (c == 'A')
		return 'T';
	else if (c == 'T')
		return 'A';
	else if (c == 'C')
		return 'G';
	else if (c == 'G')
		return 'C';
	else
		return c;
}

//// swaps two bases and complements them
//static inline void rc(char & c1, char & c2)
//{
//	char x = c1;
//	c1 = cpl(c2);
//	c2 = cpl(x);
//}

char const * MappedRead::computeReverseSeq() {
	if (RevSeq == 0) {
		//static int const qryMaxLen = Config.GetInt("qry_max_len");
		RevSeq = new char[qryMaxLen];
		memset(RevSeq, 0, qryMaxLen);

		char * fwd = Seq;
		char * rev = RevSeq + length - 1;

		for (int i = 0; i < length; ++i) {
			*rev-- = cpl(*fwd++);
			//Seq[i] = cpl(RevSeq[length - i - 1]);
		}
	}
	return RevSeq;
}

void MappedRead::DeleteReadSeq() {
	if (Seq != 0)
		delete[] Seq;
	Seq = 0;
	if (RevSeq != 0)
		delete[] RevSeq;
	RevSeq = 0;
}

MappedRead::~MappedRead() {

	if (Scores != 0) {
		delete[] Scores;
		Scores = 0;
		iScores = 0;
	}
	if(AdditionalInfo != 0) {
		delete[] AdditionalInfo; AdditionalInfo = 0;
	}
	if (Alignments != 0) {
		for (int i = 0; i < Calculated; ++i) {
			if (Alignments[i].pBuffer1 != 0)
				delete[] Alignments[i].pBuffer1;
			Alignments[i].pBuffer1 = 0;
			if (Alignments[i].pBuffer2 != 0)
				delete[] Alignments[i].pBuffer2;
			Alignments[i].pBuffer2 = 0;
			if(Alignments[i].ExtendedData != 0)
				delete[] ((AlignmentPosition *)Alignments[i].ExtendedData);
			Alignments[i].ExtendedData = 0;
		}
		delete[] Alignments;
		Alignments = 0;
	}
	DeleteReadSeq();
	if (qlty != 0) {
		delete[] qlty;
		qlty = 0;
	}
	if (name != 0)
		delete[] name;
	name = 0;

#ifdef INSTANCE_COUNTING
	AtomicDec(&sInstanceCount);
#endif
}

void MappedRead::AllocScores(LocationScore * tmp, int const n) {
	iScores = n;
	if (iScores > 0) {
		Scores = new LocationScore[iScores];
		memcpy(Scores, tmp, n * sizeof(LocationScore));
	}
}

//void MappedRead::reallocScores(int const n) {
//	if (Scores != 0) {
//		Scores = (LocationScore *) realloc(Scores, n * sizeof(LocationScore));
//	}
//}

void MappedRead::clearScores(int const TopScore) {
	if (Scores != 0) {
		if (TopScore != -1) {
			LocationScore tmp = Scores[TopScore];
			delete[] Scores;
			Scores = new LocationScore[1];
			//Scores = (LocationScore *)realloc(Scores, 1 * sizeof(LocationScore));
			Scores[0] = tmp;
			iScores = 1;
		} else {
			delete[] Scores;
			iScores = 0;
			Scores = 0;
		}
	}

}

int MappedRead::numScores() const {
	return iScores;
}

bool MappedRead::hasCandidates() const {
	return iScores > 0;
}
