#ifndef __CIGARWRITER_H__
#define __CIGARWRITER_H__

#include "GenericReadWriter.h"

class SAMWriter: public GenericReadWriter {
public:
//		SAMWriter(char const * const filename) :
//			GenericReadWriter(filename) {
	SAMWriter(FileWriter * writer) :
			GenericReadWriter(), slamSeq(Config.GetInt(SLAM_SEQ)), trimPolyA(Config.GetInt(MAX_POLYA) >= 0) {

		if(Config.Exists(RG_ID)) {
			RG = Config.GetString(RG_ID);
		} else {
			RG = 0;
		}
		m_Writer = writer;
	}

	virtual ~SAMWriter() {
		m_Writer->Flush(bufferPosition, BUFFER_LIMIT, writeBuffer, true);
	}

protected:
	virtual void DoWriteProlog();
	virtual void DoWriteRead(MappedRead const * const read, int const scoreID);
	virtual void DoWriteRead(MappedRead const * const read, int const * scoreIDs, int const scoreIdLength);
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1, MappedRead const * const read2, int const scoreId2);
	virtual void DoWriteReadGeneric(MappedRead const * const read, int const scoreID, char const * pRefName, int const pLoc, int const pDist, int const mappingQlty, int flags =
			0);
	virtual void DoWriteUnmappedReadGeneric(MappedRead const * const read, int const refId, char const pRefName, int const loc, int const pLoc, int const pDist, int const mappingQlty, int flags);
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags = 0x4);
	virtual void DoWriteEpilog();

private:

	FileWriter * m_Writer;

	char const * RG;

	int const slamSeq;

	bool const trimPolyA;
};

#endif
