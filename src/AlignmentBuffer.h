#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include "GenericReadWriter.h"

#include "SAMWriter.h"
#include "BAMWriter.h"
#include "ScoreWriter.h"

#undef module_name
#define module_name "OUTPUT"

class AlignmentBuffer {

private:

	struct Alignment {
		MappedRead * read;
		int scoreId;
	};

	int const outputformat;
	int const alignmode;
	int const batchSize;
	int const corridor;
	uloc const refMaxLen;
	int const min_mq;

	Alignment * reads;
	int nReads;
	char const * * qryBuffer;
	char const * * refBuffer;
	char const * * qalBuffer;
	char * m_DirBuffer;
	int dbLen;
	Align * alignBuffer;
	char * dBuffer;
	char * dummy;

	long pairInsertCount;
	long pairInsertSum;
	long brokenPairs;

	float alignTime;

	GenericReadWriter* m_Writer;

	static bool first;

	IAlignment * aligner;

	bool const argos;
	bool const m_EnableBS;
	int const slamSeq;

	void debugAlgnFinished(MappedRead * read);

public:

	static ulong alignmentCount;

	AlignmentBuffer(const char* const filename, IAlignment * mAligner) :
			batchSize(mAligner->GetAlignBatchSize() / 2), outputformat(
					NGM.GetOutputFormat()),
					alignmode(Config.GetInt(MODE, 0, 1)),
					corridor(Config.GetInt("corridor")),
					refMaxLen((Config.GetInt("qry_max_len") + corridor) | 1 + 1),
					min_mq(Config.GetInt(MIN_MQ)),
					aligner(mAligner), argos(Config.Exists(ARGOS)), m_EnableBS(Config.GetInt("bs_mapping", 0, 1) == 1), slamSeq(Config.GetInt(SLAM_SEQ)) {
						pairInsertSum = 0;
						pairInsertCount = 0;
						brokenPairs = 0;
						m_Writer = 0;
						nReads = 0;

						int const outputformat = NGM.GetOutputFormat();

						if(Config.Exists(ARGOS)) {
							m_Writer = (GenericReadWriter*) new ScoreWriter((FileWriter*)NGM.getWriter());
						} else {
							switch (outputformat) {
								case 0:
								Log.Error("This output format is not supported any more.");
								Fatal();
								break;
								case 1:
								m_Writer = (GenericReadWriter*) new SAMWriter((FileWriter*)NGM.getWriter());
								break;
								case 2:
								m_Writer = (GenericReadWriter*) new BAMWriter((FileWriterBam*)NGM.getWriter(), filename);
								break;
								default:
								break;
							}
						}

						if(first) {
							m_Writer->WriteProlog();
							first = false;
						}

						Log.Verbose("Alignment batchsize = %i", batchSize);

						reads = new Alignment[batchSize];

						qryBuffer = new char const *[batchSize];
						refBuffer = new char const *[batchSize];
						qalBuffer = new char const *[batchSize];

						for (int i = 0; i < batchSize; ++i) {
							refBuffer[i] = new char[refMaxLen];
						}

						m_DirBuffer = new char[batchSize];

						alignBuffer = new Align[batchSize];
						dbLen = std::max(1, Config.GetInt("qry_max_len")) * 8;
						dBuffer = new char[dbLen];

						dummy = new char[refMaxLen];
						memset(dummy, '\0', refMaxLen);
						//dummy[Config.GetInt("qry_max_len") - 1] = '\0';

						alignTime = 0.0f;

					}

					virtual ~AlignmentBuffer() {
						delete m_Writer;
						delete[] m_DirBuffer;
						m_DirBuffer = 0;

						delete[] dummy;
						dummy = 0;

						for (int i = 0; i < batchSize; ++i) {
							delete[] refBuffer[i];
							refBuffer[i] = 0;
						}
						delete[] qryBuffer;
						delete[] refBuffer;
						delete[] qalBuffer;
						delete[] alignBuffer;

						delete[] reads;
						delete[] dBuffer;

						//m_Writer->WriteEpilog();

						//delete m_Writer;
					}

					void DoRun();

					int GetStage() const {
						return 4;
					}

					inline const char* GetName() const {
						return "Output";
					}

					void addRead(MappedRead * read, int scoreID);
					void flush();

					float getTime() {
						float tmp = alignTime;
						alignTime = 0;
						return tmp;
					}

					void SaveRead(MappedRead* read, bool mapped = true);
					void WriteRead(MappedRead* read, bool mapped);
				};

#endif
