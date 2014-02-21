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

#include <iostream>

KSEQ_INIT(gzFile, gzread)

class IParser {

public:

	static size_t const MAX_READNAME_LENGTH = 100;

	int qryMaxLen;

	virtual ~IParser() {

	}

	virtual void init(char const * fileName, bool const keepTags) = 0;
	int parseRead(MappedRead * pRead) {
		assert(pRead != 0);
		return doParseRead(pRead);
	}

protected:

	virtual int doParseRead(MappedRead * pRead) = 0;

	int copyToRead(MappedRead * read, kseq_t * kseq, int const l) {

		int nameLength = 0;
		if (l >= 0) {
			if (kseq->seq.l == kseq->qual.l || kseq->qual.l == 0) {

				nameLength = std::min(MAX_READNAME_LENGTH - 1, kseq->name.l);
				memcpy(read->name, kseq->name.s, nameLength);
				read->name[nameLength] = '\0';

				//Sequence
				memset(read->Seq, '\0', qryMaxLen);
				if (kseq->seq.l != 0) {
					read->length = std::min(kseq->seq.l,
							(size_t) qryMaxLen - 1);
					int nCount = 0;
					for (int i = 0; i < read->length; ++i) {
						char c = toupper(kseq->seq.s[i]);
						if (c == 'A' || c == 'T' || c == 'C' || c == 'G') {
							read->Seq[i] = c;
						} else {
							read->Seq[i] = 'N';
							nCount += 1;
						}

					}
				} else {
					read->length = qryMaxLen - 2;
					memset(read->Seq, 'N', read->length);
					read->SetFlag(NGMNames::Empty);
				}

				//Quality
				if (kseq->qual.l > 0) {
					memcpy(read->qlty, kseq->qual.s, read->length);
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

}
;

#endif /* PARSER_H_ */
