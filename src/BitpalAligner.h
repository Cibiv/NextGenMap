#ifndef __BITPALALIGNER_H__
#define __BITPALALIGNER_H__

#include "IAlignment.h"
#include "IConfig.h"

class BitpalAligner: public IAlignment {

protected:

	int Score(char const * ref, char const * read);

public:
	BitpalAligner();
	virtual ~BitpalAligner();

//	void setCorridorSize(int corr);

	virtual int GetScoreBatchSize() const {
		return 8192;
	}
	virtual int GetAlignBatchSize() const {
		return 8192;
	}

	virtual int BatchScore(int const mode, int const batchSize,
				char const * const * const refSeqList,
				char const * const * const qrySeqList,
				char const * const * const qalSeqList, float * const results,
				void * extData);

	virtual int BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList,
			char const * const * const qalSeqList, Align * const results,
			void * extData);

	virtual int SingleAlign(int const mode, int const corridor,
			char const * const refSeq, char const * const qrySeq,
			Align & result, void * extData);

};

#endif//__BITPALALIGNER_H__
