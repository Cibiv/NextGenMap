/*
 * Parser.h
 *
 *  Created on: Aug 22, 2012
 *      Author: fritz
 */

#ifndef PARSER_H_
#define PARSER_H_

#include "MappedRead.h"

#include <zlib.h>
#include "kseq.h"
#include "../SAMRecord.h"
#include <string.h>

#include "IConfig.h"

#include <iostream>

KSEQ_INIT(gzFile, gzread)

class IParser {

public:

	static size_t const MAX_READNAME_LENGTH = 100;

	IParser(int const qrymaxlen) : qryMaxLen(qrymaxlen), trim5p(Config.GetInt(TRIM5)), trim3p(0) {

	}

	virtual ~IParser() {

	}

	virtual void init(char const * fileName, bool const keepTags) = 0;

	int parseRead(MappedRead * pRead) {
		assert(pRead != 0);
		return doParseRead(pRead);
	}
	int parseSAMRecord(SAMRecord * pRead) {
		assert(pRead != 0);
		return doParseRead(pRead);
	}


protected:

	int const qryMaxLen;
	int const trim5p;
	int const trim3p;

	virtual int doParseRead(MappedRead * pRead) = 0;
	virtual int doParseRead(SAMRecord * pRead) = 0;

	int copyToRead(MappedRead * read, kseq_t * kseq, int const l) {
		int nameLength = 0;
		if (l >= 0) {
			if (kseq->seq.l == kseq->qual.l || kseq->qual.l == 0) {

				nameLength = std::min(MAX_READNAME_LENGTH - 1, kseq->name.l);
				memcpy(read->name, kseq->name.s, nameLength);
				read->name[nameLength] = '\0';

				//Sequence
				memset(read->Seq, '\0', qryMaxLen);
				if (kseq->seq.l != 0 && kseq->seq.l > trim5p) {
					read->length = std::min(kseq->seq.l - trim5p,
							(size_t) qryMaxLen - 1);
					int nCount = 0;
					for (int i = trim5p; i < read->length + trim5p; ++i) {
						char c = toupper(kseq->seq.s[i]);
						if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
							read->Seq[i - trim5p] = c;
						} else {
							read->Seq[i - trim5p] = 'N';
							nCount += 1;
						}

					}
					read->Seq[read->length] = '\0';
				} else {
					read->length = 1;
					read->Seq[0] = 'N';
//					read->length = qryMaxLen - 2;
//					memset(read->Seq, 'N', read->length);
					read->SetFlag(NGMNames::Empty);
				}

				//Quality
				if (kseq->qual.l > 0 && kseq->qual.l > trim5p) {
					memcpy(read->qlty, kseq->qual.s + trim5p, read->length);
					read->qlty[read->length] = '\0';
				} else {
					read->qlty[0] = '*';
					read->qlty[1] = '\0';
				}
			} else {
				throw "Error while parsing. Read length not equal to length of quality values!";
				//Log.Error("Discarding read %s. Length of read not equal length of quality values.", parser->read->name.s);
			}
		} else {
			switch (l) {
			case -1:				//End of file
				break;
			case -2:
				//Length of read not equal to length of quality values
				nameLength = std::min(MAX_READNAME_LENGTH - 1, kseq->name.l);
				memcpy(read->name, kseq->name.s, nameLength);
				read->name[nameLength] = '\0';
				break;
			default:
				//Unknown error. Should not happen.
				throw "Unknown error while parsing. Please check whether the input file is corrupted!";
			}
		}
		return l;
	}

	int copyToRead(SAMRecord * read, kseq_t * kseq, int const l) {

		if (l >= 0) {
			if (kseq->seq.l == kseq->qual.l || kseq->qual.l == 0) {
				read->set_read_name(string(kseq->name.s));
				read->set_sequence(string(kseq->seq.s));
				read->set_qualities(string(kseq->qual.s));
			} else {
				throw "Error while parsing. Read length not equal to length of quality values!";
			}
		} else {
			switch (l) {
			case -1:				//End of file
				break;
			case -2:
				//Length of read not equal to length of quality values
				read->set_read_name(string(kseq->name.s));
				break;
			default:
				//Unknown error. Should not happen.
				throw "Unknown error while parsing. Please check whether the input file is corrupted!";
			}
		}
		return l;
	}


}
;

#endif /* PARSER_H_ */
