#ifndef __GENERICREADWRITER_H__
#define __GENERICREADWRITER_H__

#include <map>

#include <stdio.h>
#include <stdarg.h>

#include "ILog.h"
#include "Config.h"
#include "NGM.h"

#include "MappedRead.h"
#include "FileWriter.h"

#undef module_name
#define module_name "OUTPUT"

class GenericReadWriter {

public:

	GenericReadWriter() :
			baseNumber(5) {
		writeBuffer = new char[BUFFER_SIZE];
		bufferPosition = 0;

		identity = 0;
		if (Config.Exists("identity_tresh")) {
			identity = Config.GetFloat("identity_tresh", 0, 100);
			identity = identity / 100.0f;
		}
		writeUnmapped = !Config.GetInt("no_unal");
	}
	virtual ~GenericReadWriter() {
		if(writeBuffer != 0) {
			delete[] writeBuffer;
			writeBuffer = 0;
		}
	}
protected:

	virtual void DoWriteProlog() = 0;
	virtual void DoWriteRead(MappedRead const * const read,
			int const scoreID) = 0;
	virtual void DoWriteRead(MappedRead const * const read,
				int const * scoreIDs, int const scoreIdLength) = 0;
	virtual void DoWritePair(MappedRead const * const read1, int const scoreId1,
			MappedRead const * const read2, int const scoreId2) = 0;
	virtual void DoWriteUnmappedRead(MappedRead const * const read, int flags =
			0x4) = 0;
	virtual void DoWriteEpilog() = 0;

	float identity;

	bool writeUnmapped;

	int const baseNumber;

	static int const BUFFER_SIZE = 20000000;
	static int const BUFFER_LIMIT = 16000000;

	char * writeBuffer;
	int bufferPosition;

	int Print(const char *format, ...) {
		int done;
		va_list arg;

		va_start(arg, format);
		done = vsprintf(writeBuffer + bufferPosition, format, arg);
		bufferPosition += done;
		va_end(arg);
		return done;
	}

private:

	static int const MAX_PASSED = 1000;
	int passed[MAX_PASSED];

public:
	void WriteProlog() {
		DoWriteProlog();
	}

	int computeSlaSeqTags(MappedRead const * const read, int const scoreID,
			char * & mp, size_t & mpIndex, char * & ra, size_t & raIndex) {

		// Rates contains the number and type of matches/mismatches in the alignment
		//               Read
		//			 A	 C	 G	 T	 N
		//  	A 	 0	 1	 2	 3	 4
		// R	C 	 5	 6	 7	 8	 9
		// e	G 	10	11	12	13	14
		// f	T 	15	16	17	18	19
		//  	N 	20	21	22	23	24
		int * rates = new int[baseNumber * baseNumber + 1];
		memset(rates, 0, sizeof(int) * baseNumber * baseNumber);

		AlignmentPosition * alignmentPositons =
				(AlignmentPosition *) read->Alignments[scoreID].ExtendedData;

		// MismatchPositions
		// <type1>:<readPos1>:<refPos1>,<type2>:<readPos2>:<refPos2>,...
		// DEPRECATED: For reverse reads order should be reversed -> Last position of the read is reported as first
		mp = new char[read->length * 100];
		mpIndex = 0;

//		if (read->Scores[scoreID].Location.isReverse()) {
//			//Get number of alignment positions
//			//TODO: pass length
//			size_t alignIndex = 0;
//			while (alignmentPositons[alignIndex].type > -1) {
//				alignIndex += 1;
//			}
//
//			for (int i = alignIndex - 1; i >= 0; i--) {
//				if (alignmentPositons[i].type >= 0
//						&& alignmentPositons[i].type
//								< (baseNumber * baseNumber)) {
//					//Increase count in rates array
//					rates[alignmentPositons[i].type] += 1;
//					if (!alignmentPositons[i].match) {
//						mpIndex += sprintf(mp + mpIndex, "%d:%d,",
//								alignmentPositons[i].type,
//								read->length
//										- alignmentPositons[i].readPosition);
//					}
//
//				} else {
//					Log.Error("Error while printing T->C rates for SlamSeq output.");
//					Log.Message("%d - %d: %d", i, alignmentPositons[i].type, alignmentPositons[i].readPosition);
//					exit(1);
//				}
//
//			}
//		} else {
		// Same but for forward reads
		size_t alignIndex = 0;
		while (alignmentPositons[alignIndex].type > -1) {
			if (alignmentPositons[alignIndex].type >= 0
					&& alignmentPositons[alignIndex].type
							< (baseNumber * baseNumber)) {
				//Increase count in rates array
				rates[alignmentPositons[alignIndex].type] += 1;
				if (!alignmentPositons[alignIndex].match) {
					mpIndex += sprintf(mp + mpIndex, "%d:%d:%d,",
							alignmentPositons[alignIndex].type,
							alignmentPositons[alignIndex].readPosition + 1,
							alignmentPositons[alignIndex].refPosition + 1);
				}

			} else {
				Log.Error("Error while printing T->C rates for SlamSeq output.");
				Log.Message("%d: %d", alignmentPositons[alignIndex].type, alignmentPositons[alignIndex].readPosition);
				exit(1);
			}

			alignIndex += 1;
		}
//		}
		if (mpIndex > 0) {
			mp[mpIndex - 1] = '\0';
		}

		ra = new char[baseNumber * baseNumber * 4];
		raIndex = 0;
		for (int i = 0; i < baseNumber; ++i) {
			for (int j = 0; j < baseNumber; ++j) {
				int number = rates[i * baseNumber + j];
				raIndex += sprintf(ra + raIndex, "%d,", number);
			}
		}
		ra[raIndex] = '\0';

		int tcCount = 0;
		if (read->Scores[scoreID].Location.isReverse()) {
			tcCount = rates[baseNumber * 0 + 2];
		} else {
			tcCount = rates[baseNumber * 3 + 1];
		}
		delete[] rates;
		rates = 0;

		return tcCount;

	}

	void WriteRead(MappedRead * read, bool mapped = true) {
		if (Config.Exists(ARGOS)) {
			if (mapped) {
				DoWriteRead(read, 0);
				NGM.AddMappedRead(read->ReadId);
			} else {
				DoWriteUnmappedRead(read);
			}
		} else {
			//TODO: remove max passed limit!
			if (mapped && read->Calculated < MAX_PASSED) {
				int indexPassed = 0;
				std::map<SequenceLocation, bool> iTable;
				bool mappedOnce = false;
				for (int i = 0; i < read->Calculated; ++i) {

					float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
					float minResidues = Config.GetFloat("min_residues", 0, 1000);

					if (minResidues <= 1.0f) {
						minResidues = read->length * minResidues;
					}

					mapped = mapped && (read->Alignments[i].Identity >= minIdentity);
					mapped = mapped && ((float)(read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd) >= minResidues);

					Log.Debug(4, "READ_%d\tOUTPUT\tChecking alignment CRM_%d (read length: %d - %d - %d)\t%f >= %f\t%f >= %f", read->ReadId, i, read->length, read->Alignments[i].QStart, read->Alignments[i].QEnd, read->Alignments[i].Identity, minIdentity, (float)(read->length - read->Alignments[i].QStart - read->Alignments[i].QEnd), minResidues);

					if (mapped) {
						mappedOnce = true;
						if (iTable.find(read->Scores[i].Location) == iTable.end()) {
							iTable[read->Scores[i].Location] = true;
							Log.Debug(4, "READ_%d\tOUTPUT\tWriting alignment CRM_%d", read->ReadId, i);
							//DoWriteRead(read, i);
							passed[indexPassed++] = i;
						} else {
							Log.Debug(4, "READ_%d\tOUTPUT\tIgnoring duplicated alignment CRM_%d", read->ReadId, i);
						}
					}
				}
				//read->Calculated = indexPassed;
				//for(int i = 0; i < indexPassed; ++i) {
					//DoWriteRead(read, passed[i]);
				//}
				DoWriteRead(read, passed, indexPassed);

				if (mappedOnce) {
					Log.Debug(4, "READ_%d\tOUTPUT\tRead was mapped", read->ReadId);
					NGM.AddMappedRead(read->ReadId);
				} else {
					if(read->HasFlag(NGMNames::Empty)) {
						Log.Debug(4, "READ_%d\tOUTPUT\tRead empty (discard read)", read->ReadId);
					} else {
						Log.Debug(4, "READ_%d\tOUTPUT\tRead unmapped", read->ReadId);
						DoWriteUnmappedRead(read);
					}
				}
			} else {
				if(read->HasFlag(NGMNames::Empty)) {
					Log.Debug(4, "READ_%d\tOUTPUT\tRead empty (discard read)", read->ReadId);
				} else {
					Log.Debug(4, "READ_%d\tOUTPUT\tRead unmapped", read->ReadId);
					DoWriteUnmappedRead(read);
				}
			}
		}
	}

	void WritePair(MappedRead * const read1, int const scoreId1, MappedRead * const read2, int const scoreId2) {
		if(read1->HasFlag(NGMNames::Empty) || read2->HasFlag(NGMNames::Empty)) {
			Log.Debug(LOG_OUTPUT_DETAILS, "Empty read found in pair: %s/%s. Both reads will be discarded and not written to output.", read1->name, read2->name);
		} else {

			//TODO: fix paired end!!! MULTI MAP!

			static float const minIdentity = Config.GetFloat("min_identity", 0.0f, 1.0f);
			static float const minResiduesConfig = Config.GetFloat("min_residues", 0, 1000);
			float minResidues = minResiduesConfig;
			static int const min_mq = Config.GetInt(MIN_MQ);

			float minResidues1 = 0.0f;
			float minResidues2 = 0.0f;

			if (minResidues <= 1.0f) {
				minResidues1 = read1->length * minResidues;
				minResidues2 = read2->length * minResidues;
			} else {
				minResidues1 = minResidues2 = minResidues;
			}

			bool mapped1 = read1->hasCandidates() && read1->mappingQlty >= min_mq && (read1->Alignments[scoreId1].Identity >= minIdentity)
			&& ((float)(read1->length - read1->Alignments[scoreId1].QStart - read1->Alignments[scoreId1].QEnd) >= minResidues1);
			bool mapped2 = read2->hasCandidates() && read2->mappingQlty >= min_mq && (read2->Alignments[scoreId2].Identity >= minIdentity)
			&& ((float)(read2->length - read2->Alignments[scoreId2].QStart - read2->Alignments[scoreId2].QEnd) >= minResidues2);

			if (!mapped1) {
				read1->clearScores();
			} else {
				NGM.AddMappedRead(read1->ReadId);
			}
			if (!mapped2) {
				read2->clearScores();
			} else {
				NGM.AddMappedRead(read2->ReadId);

			}

			//Log.Message("Output paired 1: hC %d, R: %d %d, I: %f %f", read1->hasCandidates(), ((read1->length - read1->Alignments[scoreId2].QStart - read1->Alignments[scoreId2].QEnd)), minResidues, read1->Alignments[scoreId2].Identity >= minIdentity, minIdentity);
			//Log.Message("%d %d %d", read1->length, read1->Alignments[scoreId2].QStart, read1->Alignments[scoreId2].QEnd);
			//Log.Message("Output paired 2: hC %d, R: %d %d, I: %f %f", read2->hasCandidates(), ((read2->length - read2->Alignments[scoreId2].QStart - read2->Alignments[scoreId2].QEnd)), minResidues, read2->Alignments[scoreId2].Identity >= minIdentity, minIdentity);
			//Log.Message("%d %d %d", read2->length, read2->Alignments[scoreId2].QStart, read2->Alignments[scoreId2].QEnd);

			DoWritePair(read1, scoreId1, read2, scoreId2);
		}
	}

	void WriteEpilog() {
		DoWriteEpilog();
	}
}
;

#endif
