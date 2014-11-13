#include "PrefixTable.h"
//#include "NGM.h"
#include "CS.h"

//#include "SequenceLocation.h"

#include "Timing.h"

#include <stdexcept>

#undef module_name
#define module_name "PREPROCESS"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <stdio.h>

extern int lastSeqTotal;

static uint const refTabCookie = 0x1701D;

int lastSeqTotal = 0;

int CompactPrefixTable::maxPrefixFreq = 1000;

ulong CompactPrefixTable::lastPrefix;
int CompactPrefixTable::lastBin;
int CompactPrefixTable::lastPos;

uint CompactPrefixTable::skipCount;
uint CompactPrefixTable::skipBuild;


TableUnit* CompactPrefixTable::CurrentUnit;

//Used to control which kmers should be counted for index building, only locations
//that will be in the unit should also be in the index
uloc CompactPrefixTable::kmerCountMinLocation;
uloc CompactPrefixTable::kmerCountMaxLocation;

static const unsigned char ReverseTable16[] = { 0x00, 0x04, 0x08, 0x0C, 0x01, 0x05, 0x09, 0x0D, 0x02, 0x06, 0x0A, 0x0E, 0x03, 0x07, 0x0B,
		0x0F };

//static const unsigned char ReverseTable256[] = { 0, 64, 128, 192, 16, 80, 144,
//		208, 32, 96, 160, 224, 48, 112, 176, 240, 4, 68, 132, 196, 20, 84, 148,
//		212, 36, 100, 164, 228, 52, 116, 180, 244, 8, 72, 136, 200, 24, 88, 152,
//		216, 40, 104, 168, 232, 56, 120, 184, 248, 12, 76, 140, 204, 28, 92,
//		156, 220, 44, 108, 172, 236, 60, 124, 188, 252, 1, 65, 129, 193, 17, 81,
//		145, 209, 33, 97, 161, 225, 49, 113, 177, 241, 5, 69, 133, 197, 21, 85,
//		149, 213, 37, 101, 165, 229, 53, 117, 181, 245, 9, 73, 137, 201, 25, 89,
//		153, 217, 41, 105, 169, 233, 57, 121, 185, 249, 13, 77, 141, 205, 29,
//		93, 157, 221, 45, 109, 173, 237, 61, 125, 189, 253, 2, 66, 130, 194, 18,
//		82, 146, 210, 34, 98, 162, 226, 50, 114, 178, 242, 6, 70, 134, 198, 22,
//		86, 150, 214, 38, 102, 166, 230, 54, 118, 182, 246, 10, 74, 138, 202,
//		26, 90, 154, 218, 42, 106, 170, 234, 58, 122, 186, 250, 14, 78, 142,
//		206, 30, 94, 158, 222, 46, 110, 174, 238, 62, 126, 190, 254, 3, 67, 131,
//		195, 19, 83, 147, 211, 35, 99, 163, 227, 51, 115, 179, 243, 7, 71, 135,
//		199, 23, 87, 151, 215, 39, 103, 167, 231, 55, 119, 183, 247, 11, 75,
//		139, 203, 27, 91, 155, 219, 43, 107, 171, 235, 59, 123, 187, 251, 15,
//		79, 143, 207, 31, 95, 159, 223, 47, 111, 175, 239, 63, 127, 191, 255 };

///////////////////////////////////////////////////////////////////////////////
// HELPERS
///////////////////////////////////////////////////////////////////////////////

//Works only for 4 byte
inline ulong revComp(ulong prefix) {
	static const int shift = 32 - CS::prefixBits;

	//Compute complement
	ulong compPrefix = (prefix ^ 0xAAAAAAAA) & CS::prefixMask;
	//Reverse
	compPrefix = compPrefix << shift;
	ulong compRevPrefix = (ReverseTable16[compPrefix & 0x0f] << 28) | (ReverseTable16[(compPrefix >> 4) & 0x0f] << 24)
			| (ReverseTable16[(compPrefix >> 8) & 0x0f] << 20) | (ReverseTable16[(compPrefix >> 12) & 0x0f] << 16)
			| (ReverseTable16[(compPrefix >> 16) & 0x0f] << 12) | (ReverseTable16[(compPrefix >> 20) & 0x0f] << 8)
			| (ReverseTable16[(compPrefix >> 24) & 0x0f] << 4) | (ReverseTable16[(compPrefix >> 28) & 0x0f]);

//	ulong compRevPrefix = (ReverseTable256[compPrefix & 0xff] << 24)
//			| (ReverseTable256[(compPrefix >> 8) & 0xff] << 16)
//			| (ReverseTable256[(compPrefix >> 16) & 0xff] << 8)
//			| (ReverseTable256[(compPrefix >> 24) & 0xff]);

	return compRevPrefix;
}

static inline int calc_binshift(int corridor) {
	corridor >>= 1;
	int l = 0;
	while ((corridor >>= 1) > 0)
		++l;
	return l;
}

inline int GetBin(uint pos) {
	static int shift = calc_binshift(12);
	return pos >> shift;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

CompactPrefixTable::CompactPrefixTable(bool const dualStrand, bool const skip) :
		m_RECount(0), m_RRCount(0), m_BCalls(0), m_TotalLocs(0), DualStrand(dualStrand), skipRep(skip) {

	bool m_EnableBS = false;
	m_EnableBS = (Config.GetInt("bs_mapping", 0, 1) == 1);

	maxPrefixFreq = (Config.Exists("kmer_maxFreq") ? Config.GetInt("kmer_maxFreq") : maxPrefixFreq);
	m_RefSkip = (Config.Exists("kmer_skip") ? Config.GetInt("kmer_skip", 0, -1) : 0);
	if (m_EnableBS) {
		m_RefSkip = 0;
		Log.Verbose("BS mapping enabled. Kmer skip on ref is set to 0");
	}
	m_PrefixLength = CS::prefixBasecount;
	uint indexLength = (int) pow(4.0, (double) m_PrefixLength);

	std::stringstream refFileName;
	refFileName << std::string(Config.GetString("ref")) << "-ht-" << m_PrefixLength << "-" << m_RefSkip << ".ngm";

	char * cacheFile = new char[refFileName.str().size() + 1];
	strcpy(cacheFile, refFileName.str().c_str());

	if (!FileExists(cacheFile)) {
		kmerCountMinLocation = 0;
		kmerCountMaxLocation = c_tableLocMax;

		uint genomeSize = SequenceProvider.GetTotalLen(); //TODO: Does this function work as expected?
		m_UnitCount = 1 + genomeSize / c_tableLocMax;
		m_Units = new TableUnit[ m_UnitCount ];

		Log.Message("Allocated %d hashtable units",m_UnitCount);

		CreateTable(indexLength);
		
		saveToFile(cacheFile, indexLength);
	} else {
		readFromFile(cacheFile);
	}
}

CompactPrefixTable::~CompactPrefixTable() {
	Log.Verbose("Clearing prefix table");
	Clear();
	Log.Verbose("Cleanup done");
}

void CompactPrefixTable::Clear() {
	delete[] m_Units;
}

int * CompactPrefixTable::CountKmerFreq(uint length) {

	Log.Message("Number of k-mers: %d", length);
	int * freq = new int[length];
	memset(freq, 0, length);

	for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
		lastPrefix = 111111;
		lastBin = -1;
		m_CurGenSeq = i;

		if (!DualStrand || !(m_CurGenSeq % 2)) {

			uint offset = SequenceProvider.GetRefStart(m_CurGenSeq);
			uint len = SequenceProvider.GetRefLen(m_CurGenSeq);
			char * seq = new char[len + 2];
			SequenceProvider.DecodeRefSequence(seq, m_CurGenSeq, offset, len);

			if(skipRep) {
				CS::PrefixIteration(seq, len, &CompactPrefixTable::CountKmer, 0, 0, freq, m_RefSkip, offset);
			} else {
				CS::PrefixIteration(seq, len, &CompactPrefixTable::CountKmerwoSkip, 0, 0, freq, m_RefSkip, offset);
			}
			delete[] seq;
			seq = 0;
		}
	}
	return freq;
}

void CompactPrefixTable::Generate() {

	int i = 0;
	for (int i = 0; i < SequenceProvider.GetRefCount(); ++i) {
		lastPrefix = 111111;
		lastBin = -1;

		m_CurGenSeq = i;

		if (!DualStrand || !(m_CurGenSeq % 2)) {
			Timer t;
			t.ST();

			uint offset = SequenceProvider.GetRefStart(m_CurGenSeq);
			uint len = SequenceProvider.GetRefLen(m_CurGenSeq);
			char * seq = new char[len + 2];
			SequenceProvider.DecodeRefSequence(seq, m_CurGenSeq, offset, len);

			if(skipRep) {
				CS::PrefixIteration(seq, len, &CompactPrefixTable::BuildPrefixTable, 0, 0, this, m_RefSkip, offset);
			} else {
				CS::PrefixIteration(seq, len, &CompactPrefixTable::BuildPrefixTablewoSkip, 0, 0, this, m_RefSkip, offset);
			}
			Log.Verbose("Create table for chr %d. Start: %d, Length: %u (%.2fs)", m_CurGenSeq, 0, len, t.ET());
			delete[] seq;
			seq = 0;
		}
	}

	if (skipBuild == skipCount) {
		Log.Message("Number of repetitive k-mers ignored: %d", skipBuild);
	} else {
		Log.Error("SkipBuild (%d) != SkipCount (%d)", skipCount, skipBuild);
	}
}

uint CompactPrefixTable::createRefTableIndex(uint const length) {

	Timer freqT;
	freqT.ST();
	int * freqs = CountKmerFreq(length);
	Log.Message("Counting kmers took %.2fs", freqT.ET());

	Timer t;
	t.ST();


	CurrentUnit->RefTableIndex = new Index[length];

	uint next = 0;
	int ignoredPrefixes = 0;
	int usedPrefixes = 0;

	int maxFreq = 0;
	long sum = 0;

	uint i = 0;
	for (i = 0; i < length; i++) {
		//Add for each kmer ref to the reverse complement kmer
		ulong compRevPrefix = revComp(i);

		//Create index based on kmer frequencies
		int freq = freqs[i];
		int total_freq = freq + freqs[compRevPrefix];
		maxFreq = std::max(maxFreq, total_freq);

		if (freq > 0 && total_freq < maxPrefixFreq) {
			CurrentUnit->RefTableIndex[i].m_TabIndex = next + 1;
			CurrentUnit->RefTableIndex[i].m_RevCompIndex = (maxPrefixFreq - total_freq) * 100.0f / maxPrefixFreq;
			next += freq;
			sum += freq;
			usedPrefixes += 1;
		} else {
			CurrentUnit->RefTableIndex[i].m_TabIndex = next + 1;
			if (freq > 0) {
				ignoredPrefixes += 1;
			}
		}
	}
	CurrentUnit->RefTableIndex[i].m_TabIndex = next + 1;
	float avg = sum * 1.0 / (usedPrefixes + ignoredPrefixes) * 1.0;
	delete[] freqs;
	freqs = 0;

	Log.Message("Average number of positions per prefix: %f", avg);
	Log.Message("%d prefixes are ignored due to the frequency cutoff (%d)", ignoredPrefixes, maxPrefixFreq);
	Log.Message("Index size: %d byte (%d x %d)", length * sizeof(Index), length, sizeof(Index));
	Log.Message("Generating index took %.2fs", t.ET());
	return next;
}

void CompactPrefixTable::CreateTable(uint const length) {
	for( int i = 0; i < m_UnitCount; ++ i )
	{
		CurrentUnit = &m_Units[ i ];

		Log.Message("Building RefTable (kmer length: %d, reference skip: %d)", m_PrefixLength, m_RefSkip);
		Timer gtmr;
		gtmr.ST();
		CurrentUnit->cRefTableLen = createRefTableIndex(length);

		Timer tmr;
		tmr.ST();

		CurrentUnit->RefTable = new Location[CurrentUnit->cRefTableLen + 1];

		for (uint i = 0; i < CurrentUnit->cRefTableLen + 1; ++i) {
			CurrentUnit->RefTable[i].m_Location = 0;
		}
		Log.Message("Allocating and initializing prefix Table took %.2fs", tmr.ET());
		Log.Message("Number of prefix positions is %d (%d)", CurrentUnit->cRefTableLen, sizeof(Location));
		Log.Message("Size of RefTable is %ld", (ulong)CurrentUnit->cRefTableLen * (ulong)sizeof(Location));

		Generate();
		Log.Message("Overall time for creating RefTable: %.2fs", gtmr.ET());

		kmerCountMinLocation += c_tableLocMax;
		kmerCountMaxLocation += c_tableLocMax;
	}

//	long count = 0;
//	RefEntry * dummy = new RefEntry(0);
//	for (uint i = 0; i < length; i++) {
//		RefEntry const * entry = GetRefEntry(i, dummy);
//		int lastBin = -1;
//		for (int j = 0; j < entry->refCount; ++j) {
//			int currentBin = GetBin(entry->ref[j].m_Location);
//			if (currentBin == lastBin) {
////				entry->ref[j].m_RefId = -1;
//				count += 1;
//				Log.Message("Prefix %d:\t%d\t%d (%d)\t(%ld)", i, entry->ref[j].m_RefId, entry->ref[j].m_Location, currentBin, count);
//			}
//			lastBin = currentBin;
//		}
////		Log.Message("-----------------------------------------------");
//	}

}

void CompactPrefixTable::CountKmer(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo, void* data) {
	if( pos < kmerCountMinLocation || pos > kmerCountMaxLocation )
		return;

	int * freq = (int *) data;
	if (prefix == lastPrefix) {
		int currentBin = GetBin(pos);
		if (currentBin != lastBin || lastBin == -1) {
			freq[prefix] += 1;
		} else {
			skipCount += 1;
//			Log.Message("Prefix %d (skip):\t%d (%d)\t%d (%d)\t(%ld)", prefix, lastPos, lastBin, pos, currentBin, skipCount);
		}
		lastBin = currentBin;
		lastPos = pos;
	} else {
		lastBin = -1;
		lastPos = -1;
		freq[prefix] += 1;
	}
	lastPrefix = prefix;
}

void CompactPrefixTable::CountKmerwoSkip(ulong prefix, uloc pos, ulong mutateFrom, ulong mutateTo, void* data) {
	if( pos < kmerCountMinLocation || pos > kmerCountMaxLocation )
		return;

	int * freq = (int *) data;
	freq[prefix] += 1;

}

void CompactPrefixTable::BuildPrefixTable(ulong prefix, uloc real_pos, ulong mutateFrom, ulong mutateTo, void* data) {
	if( real_pos < kmerCountMinLocation || real_pos > kmerCountMaxLocation )
		return;

	//Rebase position using current hashtable unit offset
	uloc temp_pos = real_pos - CurrentUnit->Offset;
	uint reduced_pos( temp_pos );

	CompactPrefixTable * _this = (CompactPrefixTable*) data;
	_this->m_BCalls++;

	if (prefix == lastPrefix) {
		int currentBin = GetBin(real_pos);
		if (currentBin != lastBin || lastBin == -1) {
			if (CurrentUnit->RefTableIndex[prefix].used()) {
				Location tmp = {reduced_pos};
				_this->SaveToRefTable(prefix, tmp);
			}
		} else {
			skipBuild += 1;
//			Log.Message("Prefix %d (skip):\t%d (%d)\t%d (%d)\t(%ld)", prefix, lastPos, lastBin, pos, currentBin, skipCount);
		}
		lastBin = currentBin;
		lastPos = real_pos;
	} else {
		lastBin = -1;
		lastPos = -1;
		if (CurrentUnit->RefTableIndex[prefix].used()) {
			Location tmp = {reduced_pos};
			_this->SaveToRefTable(prefix, tmp);
		}
	}
	lastPrefix = prefix;
}

void CompactPrefixTable::BuildPrefixTablewoSkip(ulong prefix, uloc real_pos, ulong mutateFrom, ulong mutateTo, void* data) {
	if( real_pos < kmerCountMinLocation || real_pos > kmerCountMaxLocation )
		return;

	//Rebase position using current hashtable unit offset
	uloc temp_pos = real_pos - CurrentUnit->Offset;
	uint reduced_pos( temp_pos );

	CompactPrefixTable * _this = (CompactPrefixTable*) data;
	_this->m_BCalls++;

	if (CurrentUnit->RefTableIndex[prefix].used()) {
		Location tmp = {reduced_pos};
		_this->SaveToRefTable(prefix, tmp);
	}
}

void CompactPrefixTable::SaveToRefTable(ulong prefix, Location loc) {

	uint start = CurrentUnit->RefTableIndex[prefix].m_TabIndex - 1;
	uint maxLength = CurrentUnit->RefTableIndex[prefix + 1].m_TabIndex - 1 - start;

	uint i = 0;

	while (CurrentUnit->RefTable[start + i].used() && i < maxLength) {
		i += 1;
	}
 if (CurrentUnit->RefTable[start + i].used()) {
		Log.Message(
				"Tried to insert kmer %d starting at position %d, number of slots %d. Position: %d",
				prefix, start, maxLength, i);
		throw;
	} else {
		CurrentUnit->RefTable[start + i] = loc;
	}

	++m_TotalLocs;
}

RefEntry const * CompactPrefixTable::GetRefEntry(ulong prefix, RefEntry * initial_entry) const {

	RefEntry* entry = initial_entry;

	for( int i = 0; i < m_UnitCount; ++ i )
	{
		uint      cRefTableLen = m_Units[i].cRefTableLen;
		Location* RefTable = m_Units[i].RefTable;
		Index*    RefTableIndex = m_Units[i].RefTableIndex;

		uint start = 0;
		uint maxLength = 0;
		if (RefTableIndex[prefix].used()) {
			start = RefTableIndex[prefix].m_TabIndex - 1;
			//TODO: Fix Invalid read of size 4
			maxLength = RefTableIndex[prefix + 1].m_TabIndex - 1 - start;

			entry->ref = RefTable + start;
			entry->reverse = false;
			entry->weight = RefTableIndex[prefix].m_RevCompIndex;
	//		entry->weight = 1.0f;
			entry->refCount = maxLength;
			entry->refTotal = maxLength;
			entry->offset = m_Units[i].Offset;
		} else {
			entry->weight = 0.0f;
			entry->refCount = 0;
			entry->refTotal = 0;
			entry->offset = m_Units[i].Offset;
		}

		ulong compRevPrefix = revComp(prefix);
		RefEntry * revEntry = entry->nextEntry;

		if (RefTableIndex[compRevPrefix].used()) {
			start = RefTableIndex[compRevPrefix].m_TabIndex - 1;
			//TODO: Fix Invalid read of size 4
			maxLength = RefTableIndex[compRevPrefix + 1].m_TabIndex - 1 - start;

			revEntry->ref = RefTable + start;
			revEntry->reverse = true;
			revEntry->weight = RefTableIndex[compRevPrefix].m_RevCompIndex;
	//		revEntry->weight = 1.0f;
			revEntry->refCount = maxLength;
			entry->refTotal = revEntry->refTotal = entry->refTotal + maxLength;
			revEntry->offset = m_Units[i].Offset;
		} else {
			revEntry->weight = 0.0f;
			revEntry->refCount = 0;
			revEntry->refTotal = 0;
			revEntry->offset = m_Units[i].Offset;
		}
		//}

		if( m_UnitCount > 1 )
		{
			RefEntry* new_sense = new RefEntry(0);
			RefEntry* new_antisense = new RefEntry(0);

			revEntry->nextEntry = new_sense;
			new_sense->nextEntry = new_antisense;

			entry = new_sense;
		}
	}

	return initial_entry;
}

void CompactPrefixTable::saveToFile(char const * fileName, uint const refIndexSize) {
	if (!Config.GetInt("skip_save")) {
		Timer wtmr;
		wtmr.ST();
		Log.Message("Writing RefTable to %s", fileName);
		FILE *fp;
		fp = fopen(fileName, "wb");
		if (fp != 0) {
			fwrite(&refTabCookie, sizeof(uint), 1, fp);
			fwrite(&m_PrefixLength, sizeof(uint), 1, fp);
			fwrite(&m_RefSkip, sizeof(uint), 1, fp);
			fwrite(&m_UnitCount, sizeof(int), 1, fp);
			fwrite(&refIndexSize, sizeof(uint), 1, fp);

			for( int i = 0; i < m_UnitCount; ++ i )
			{
				TableUnit& curr = m_Units[ i ];

				fwrite(&curr.cRefTableLen, sizeof(uint), 1, fp);
				fwrite(curr.RefTableIndex, sizeof(Index), refIndexSize, fp);
				fwrite(curr.RefTable, sizeof(Location), curr.cRefTableLen, fp);
				fwrite(&curr.Offset, sizeof(uint), 1, fp);
			}

			fclose(fp);
		} else {
			Log.Error("Error while opening file %s for writing.", fileName);
			Fatal();
		}
		Log.Message("Writing to disk took %.2fs", wtmr.ET());
	} else {
		Log.Warning("RefTable is not saved to disk! (--skip-save)");
	}
}

void CompactPrefixTable::readFromFile(char const * fileName) {
	Log.Message("Reading RefTable from %s", fileName);
	Timer wtmr;
	wtmr.ST();
	uint refIndexSize = 0;
	uint refTableSize = 0;
	uint prefixBasecount = 0;
	uint refskip = 0;
	uint cookie = 0;
	FILE *fp;
	fp = fopen(fileName, "rb");
	if (!fp) {
		Log.Error("Couldn't open file %s for reading.", fileName);
		Fatal();
	}
	fread(&cookie, sizeof(uint), 1, fp);
	fread(&prefixBasecount, sizeof(uint), 1, fp);
	fread(&refskip, sizeof(uint), 1, fp);
	if (cookie != refTabCookie || prefixBasecount != m_PrefixLength || refskip != m_RefSkip) {
		fclose(fp);
		Log.Error("Invalid reference table found: %s.", fileName);
		Log.Error("Please delete it and run NGM again.");
		Fatal();
	}

	fread(&m_UnitCount, sizeof(int), 1, fp );
	fread(&refIndexSize, sizeof(uint), 1, fp);

	m_Units = new TableUnit[ m_UnitCount ];

	for( int i = 0; i < m_UnitCount; ++ i )
	{
		TableUnit& curr = m_Units[ i ];

		fread(&curr.cRefTableLen, sizeof(uint), 1, fp);
		curr.RefTableIndex = new Index[refIndexSize];
		fread(curr.RefTableIndex, sizeof(Index), refIndexSize, fp);
		curr.RefTable = new Location[curr.cRefTableLen + 1];
		fread(curr.RefTable, sizeof(Location), curr.cRefTableLen, fp);
		fread(&curr.Offset, sizeof(uint), 1, fp);
	}

	fclose(fp);
	Log.Message("Reading from disk took %.2fs", wtmr.ET());
}
