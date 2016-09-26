/*
 * SWOclCigar.cpp
 *
 *  Created on: May 25, 2011
 *      Author: philipp_
 */

#include "SWOclCigar.h"

#include <stdio.h>
#include <string.h>
#include <sstream>

#include "Timing.h"
#include "OclHost.h"
#include <pthread.h>

pthread_mutex_t mutext_batch_align;

using std::stringstream;

#ifndef ECLIPSE
#include "oclSwCigar.h"
#else
char const * const oclSwCigar = "dummy";
#endif

#define pCigar pBuffer1
#define pMD pBuffer2

int const result_number = 4;

SWOclCigar::SWOclCigar(OclHost * host) :
		SWOcl(oclSwCigar,
				" -D result_number=4 -D CIGAR_M=0 -D CIGAR_I=1 -D CIGAR_D=2 -D CIGAR_N=3 -D CIGAR_S=4 -D CIGAR_H=5 -D CIGAR_P=6 -D CIGAR_EQ=7 -D CIGAR_X=8 ",
				host) {
	batch_size_align = computeAlignmentBatchSize();
	Log.Debug(LOG_INFO, "Batchsize (alignment): %d", batch_size_align);
	swAlignScoreKernel = host->setupKernel(clProgram, "oclSW_Score");
	swAlignScoreKernelGlobal = host->setupKernel(clProgram,
			"oclSW_ScoreGlobal");
	swAlignBacktrackingKernel = host->setupKernel(clProgram,
			"oclSW_Backtracking");
}

SWOclCigar::~SWOclCigar() {
	clReleaseKernel(swAlignBacktrackingKernel);
	clReleaseKernel(swAlignScoreKernel);
	clReleaseKernel(swAlignScoreKernelGlobal);
}

void SWOclCigar::runSwBatchKernel(cl_kernel swScoreAlign, const int batchSize,
		const char * const * const qrySeqList,
		const char * const * const refSeqList, char * bsDirection,
		cl_mem & results_gpu, cl_mem & alignments, short * const result,
		short * calignments, cl_mem & matrix_gpu, cl_mem & bsdirection_gpu) {
	const size_t cnDim = batch_size_align / alignments_per_thread;
	const size_t cBlockSize = threads_per_block;

	int runbatchSize = std::min(batch_size_align, batchSize);
	copySeqDataToDevice(cpu_read_data, cpu_ref_data, qrySeqList, refSeqList,
			runbatchSize, batch_size_align);
	for (int i = 0; i < batchSize; i += batch_size_align) {

		if (host->isGPU()) {
			host->writeToDevice(scaffold_gpu, CL_FALSE, 0,
					ref_data_size * sizeof(cl_char), cpu_ref_data);
			host->writeToDevice(reads_gpu, CL_FALSE, 0,
					read_data_size * sizeof(cl_char), cpu_read_data);

			if (bsdirection_gpu != 0) {
				host->writeToDevice(bsdirection_gpu, CL_FALSE, 0,
						std::min(runbatchSize, batchSize - i) * sizeof(cl_char),
						bsDirection + i);
			}

			host->waitForDevice();
			host->executeKernel(interleaveKernel, cnDim, cBlockSize);
		}

		host->executeKernel(swScoreAlign, cnDim, cBlockSize);
		host->executeKernel(swAlignBacktrackingKernel, cnDim, cBlockSize);

		int nextBatch = (i + batch_size_align);
		int nextRunbatchSize = std::min(batch_size_align,
				batchSize - nextBatch);
		if (nextRunbatchSize > 0) {
			copySeqDataToDevice(cpu_read_data, cpu_ref_data,
					qrySeqList + nextBatch, refSeqList + nextBatch,
					nextRunbatchSize, batch_size_align);
		}
		host->readFromDevice(results_gpu, CL_FALSE, 0,
				result_number * runbatchSize, result + i * result_number,
				sizeof(cl_short));
		host->readFromDevice(alignments, CL_FALSE, 0,
				runbatchSize * alignment_length * 2,
				calignments + alignment_length * 2 * i, sizeof(cl_short));
		runbatchSize = nextRunbatchSize;

	}
	host->waitForDevice();
}

int SWOclCigar::BatchAlign(int const mode, int const batchSize_,
		char const * const * const refSeqList_,
		char const * const * const qrySeqList_,
		char const * const * const qalSeqList, Align * const results,
		void * extData) {
	if (batchSize_ <= 0) {
		Log.Warning("Align for batchSize <= 0");
		return 0;
	}

	bool batchSizeDif = !host->isGPU() && (batchSize_ % 4 != 0);
//	int batchSize = (batchSizeDif) ? batchSize_ + 4 : batchSize_;
	int batchSize =
			(batchSizeDif) ? ((int) (batchSize_ / 4) + 1) * 4 : batchSize_;
	char const * * const tmpRefSeqList = new char const *[batchSize];
	char const * * const tmpQrySeqList = new char const *[batchSize];

	if (batchSizeDif) {
		for (int i = 0; i < batchSize_; ++i) {
			tmpRefSeqList[i] = refSeqList_[i];
			tmpQrySeqList[i] = qrySeqList_[i];
		}
		for (int i = batchSize_; i < batchSize; ++i) {
			tmpRefSeqList[i] = refSeqList_[0];
			tmpQrySeqList[i] = qrySeqList_[0];
		}
	}
	char const * const * const refSeqList =
			batchSizeDif ? tmpRefSeqList : refSeqList_;
	char const * const * const qrySeqList =
			batchSizeDif ? tmpQrySeqList : qrySeqList_;

//	//	char * const * const qrySeqList = (char * const * const)qrySeqList_;
//	char * * const qrySeqList = new char *[batchSize * 2];
//
//	long corridor = (Config.GetInt("corridor"));
//	long qry_len = (Config.GetInt("qry_max_len"));
//
//	for (int i = 0; i < batchSize; ++i) {
//		std::cout << "Read: ";
//		for (int j = 0; j < (qry_len); ++j) {
//			if (qrySeqList[i][j] != '\0') {
//				std::cout << qrySeqList[i][j];
//			} else {
//				std::cout << "X";
//			}
//		}
//		std::cout << std::endl;
//		std::cout << "Ref:  ";
//		for (int j = 0; j < (corridor + qry_len); ++j) {
//			if (refSeqList[i][j] != '\0') {
//				std::cout << refSeqList[i][j];
//			} else {
//				std::cout << "X";
//			}
//		}
//		std::cout << std::endl;
//		//		std::cout << std::endl;
//		//		for (int j = 0; j < (corridor + qry_len); ++j) {
//		//			std::cout << refSeqList[i][j];
//		//		}
//		//		std::cout << std::endl;
//	}
//	std::cout << "===========================================================" << std::endl;
//
//	bool found = false;
//	for (int i = 0; i < batchSize; ++i) {
//		qrySeqList[i] = new char[qry_len + 1];
//		found = false;
//		for (int j = 0; j < (qry_len); ++j) {
//			found = found || (qrySeqList_[i][j] == '\0');
//			if (found) {
//				qrySeqList[i][j] = '\0';
//			} else {
//				qrySeqList[i][j] = qrySeqList_[i][j];
//			}
//			//			std::cout << qrySeqList[i][j];
//		}
//		//		std::cout << std::endl;
//		//		for (int j = 0; j < (corridor + qry_len); ++j) {
//		//			std::cout << refSeqList[i][j];
//		//		}
//		//		std::cout << std::endl;
//	}

//	Log.Error("Batch align: %d %d", batchSize, Config.GetInt("qry_max_len"));

	if (host->isGPU()) {
		pthread_mutex_lock(&mutext_batch_align);
	}

	cl_kernel scoreKernel;
	switch ((mode & 0xFF)) {
	case 0: {

		Log.Debug(LOG_INFO, "Alignment mode: local");

		scoreKernel = swAlignScoreKernel;
	}
		break;
	case 1:
		Log.Debug(LOG_INFO, "Alignment mode: end-free");
		scoreKernel = swAlignScoreKernelGlobal;
		break;
	default:
		Log.Error("Unsupported alignment mode %i", mode & 0xFF);
		exit(-1);
	}

	Timer timer;
	timer.ST();

	cl_int ciErrNum = 0;
	cl_mem c_scaff_gpu = scaffold_gpu;
	if (host->isGPU()) {
		c_scaff_gpu = host->allocate(CL_MEM_READ_WRITE,
				ref_data_size * sizeof(cl_char));
//#pragma omp critical
		{
			ciErrNum |= clSetKernelArg(interleaveKernel, 0, sizeof(cl_mem),
					(void *) (&scaffold_gpu));
			ciErrNum |= clSetKernelArg(interleaveKernel, 1, sizeof(cl_mem),
					&c_scaff_gpu);
		}
	} else {
		c_scaff_gpu = scaffold_gpu;
	}

	cl_mem results_gpu = host->allocate(CL_MEM_READ_WRITE,
			result_number * batch_size_align * sizeof(cl_short));
	cl_mem matrix_gpu = host->allocate(CL_MEM_READ_WRITE,
			batch_size_align * (Config.GetInt("corridor") + 2)
					* (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char));

	cl_mem alignments_gpu = host->allocate(CL_MEM_READ_WRITE,
			batch_size_align * alignment_length * 2 * sizeof(cl_short));

//#pragma omp critical
	{
		//Set parameter
		ciErrNum |= clSetKernelArg(scoreKernel, 0, sizeof(cl_mem),
				(void *) (&c_scaff_gpu));
		ciErrNum |= clSetKernelArg(scoreKernel, 1, sizeof(cl_mem),
				(void *) (&reads_gpu));
		ciErrNum |= clSetKernelArg(scoreKernel, 2, sizeof(cl_mem),
				&results_gpu);
		ciErrNum |= clSetKernelArg(scoreKernel, 3, sizeof(cl_mem), &matrix_gpu);
	}
	cl_mem bsdirection_gpu = 0;
	static bool const bsMapping = Config.GetInt("bs_mapping") == 1;
	if (bsMapping || (slamSeq & 0x2)) {
		if (host->isGPU()) {
			bsdirection_gpu = host->allocate(CL_MEM_READ_ONLY,
					batch_size_align * sizeof(cl_char));
		} else {
			bsdirection_gpu = host->allocate(
			CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
					batch_size_align * sizeof(cl_char), (char *) extData);
		}
//#pragma omp critical
		{
			ciErrNum |= clSetKernelArg(scoreKernel, 4, sizeof(cl_mem),
					(void *) (&bsdirection_gpu));
		}
	}

//#pragma omp critical
	{
		ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 0, sizeof(cl_mem),
				(void *) (&c_scaff_gpu));
		ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 1, sizeof(cl_mem),
				(void*) (&reads_gpu));
		ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 2, sizeof(cl_mem),
				&results_gpu);
		ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 3, sizeof(cl_mem),
				&matrix_gpu);
		ciErrNum |= clSetKernelArg(swAlignBacktrackingKernel, 4, sizeof(cl_mem),
				&alignments_gpu);
	}

	host->checkClError("Unable to set kernel parameters", ciErrNum);

	short * calignments = new short[batchSize * alignment_length * 2];
	short * gpu_return_values = new short[result_number * batchSize];

	runSwBatchKernel(scoreKernel, batchSize, qrySeqList, refSeqList,
			(char *) extData, results_gpu, alignments_gpu, gpu_return_values,
			calignments, matrix_gpu, bsdirection_gpu);

	for (int i = 0; i < batchSize_; ++i) {
		short * gpuCigar = calignments + i * alignment_length * 2;
		//results[i].pCigar = new char[alignment_length];
		//results[i].pMD = new char[alignment_length];
		int offset = (i / alignments_per_thread) * result_number
				* alignments_per_thread;
		int k = (i) % alignments_per_thread;

		char bsFrom = '0';
		char bsTo = '0';
		if (bsMapping) {
			if (((char *) extData)[i] == 1) {
				bsFrom = 'A';
				bsTo = 'G';
			} else {
				bsFrom = 'T';
				bsTo = 'C';
			}
		}
		if (slamSeq) {
			if (((char *) extData)[i] == 1) {
				bsFrom = 'G';
				bsTo = 'A';
			} else {
				bsFrom = 'C';
				bsTo = 'T';
			}
		}

		if (!computeCigarMD(results[i],
				gpu_return_values[offset + alignments_per_thread * 3 + k],
				gpuCigar, refSeqList[i] + gpu_return_values[offset + k],
				qrySeqList[i], qalSeqList[i], bsFrom, bsTo)) {
			results[i].Score = -1.0f;
		}
		results[i].PositionOffset = gpu_return_values[offset + k];
	}

	delete[] calignments;
	calignments = 0;
	delete[] gpu_return_values;
	gpu_return_values = 0;

	Log.Verbose("Releasing results.");
	clReleaseMemObject(results_gpu);

	Log.Verbose("Releasing alignments.");
	clReleaseMemObject(alignments_gpu);

	Log.Verbose("Releasing matrix.");
	clReleaseMemObject(matrix_gpu);

	if (host->isGPU()) {
		Log.Verbose("Releasing scaff.");

		clReleaseMemObject(c_scaff_gpu);
		if (bsMapping || (slamSeq & 0x2)) {
			clReleaseMemObject(bsdirection_gpu);
		}
	}

	if (host->isGPU()) {
		pthread_mutex_unlock(&mutext_batch_align);
	}

//	//TODO: remove
//	for (int i = 0; i < batchSize; ++i) {
//		delete[] qrySeqList[i];
//	}
//	delete[] qrySeqList;
	Log.Verbose("SW finished computing alignments for %d sequences (elapsed: %.3fs)", batchSize_, timer.ET());

	delete[] tmpRefSeqList;
	delete[] tmpQrySeqList;

	return batchSize_;

}

int printCigarElement(char const op, short const length, char * cigar) {
	int offset = 0;
	offset = sprintf(cigar, "%d%c", length, op);
	return offset;
}

bool debugCigar(int op, int length) {
	switch (op) {
	case CIGAR_M:
		Log.Message("CIGAR: %d M", length);
		break;
	case CIGAR_I:
		Log.Message("CIGAR: %d I", length);
		break;
	case CIGAR_D:
		Log.Message("CIGAR: %d D", length);
		break;
	case CIGAR_N:
		Log.Message("CIGAR: %d N", length);
		break;
	case CIGAR_S:
		Log.Message("CIGAR: %d S", length);
		break;
	case CIGAR_H:
		Log.Message("CIGAR: %d H", length);
		break;
	case CIGAR_P:
		Log.Message("CIGAR: %d P", length);
		break;
	case CIGAR_EQ:
		Log.Message("CIGAR: %d EQ", length);
		break;
	case CIGAR_X:
		Log.Message("CIGAR: %d X", length);
		break;
	default:
		//Log.Error("Invalid cigar operator.");
		//exit(1);
		return false;
	}
	return true;
}

int const baseNumber = 5;
int trans[256] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1,
		4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

char back[baseNumber] = { 'A', 'C', 'G', 'T', 'N' };

bool SWOclCigar::computeCigarMD(Align & result, int const gpuCigarOffset,
		short const * const gpuCigar, char const * const refSeq,
		char const * const qrySeq, char const * const qalSeq, char const bsFrom,
		char const bsTo) {

	static bool const bsMapping = Config.GetInt("bs_mapping") == 1;
	int cigar_offset = 0;
	int md_offset = 0;

	result.QStart = 0;
	result.QEnd = 0;

	AlignmentPosition * alignmentPositions = 0;
	if (slamSeq) {
		//Read length would be sufficient
		alignmentPositions = new AlignmentPosition[alignment_length];
		result.ExtendedData = (void *) alignmentPositions;
	}

	if ((gpuCigar[gpuCigarOffset] >> 4) > 0) {
		if (Config.GetInt("hard_clip") == 1) { //soft clipping
			cigar_offset += printCigarElement('H',
					gpuCigar[gpuCigarOffset] >> 4,
					result.pCigar + cigar_offset);
		} else if (Config.GetInt("silent_clip") != 1) {
			cigar_offset += printCigarElement('S',
					gpuCigar[gpuCigarOffset] >> 4,
					result.pCigar + cigar_offset);
		}
		result.QStart = gpuCigar[gpuCigarOffset] >> 4;
	}

//	if ((gpuCigar[gpuCigarOffset] >> 4) > 0) {
//		cigar_offset += printCigarElement('S', gpuCigar[gpuCigarOffset] >> 4, result.pCigar + cigar_offset);
//		result.QStart = gpuCigar[gpuCigarOffset] >> 4;
//	} else {
//		result.QStart = 0;
//	}

	int match = 0;
	int mismatch = 0;
	int total = 0;
	int cigar_m_length = 0;
	int md_eq_length = 0;
	int ref_index = 0;
	int read_index = result.QStart;
	for (int j = gpuCigarOffset + 1; j < (alignment_length - 1); ++j) {
		int op = gpuCigar[j] & 15;
		int length = gpuCigar[j] >> 4;

		//debugCigar(op, length);
		total += length;
		switch (op) {
		case CIGAR_X:
			if (slamSeq) {
				for (int k = 0; k < length; ++k) {
					// Log.Message("%c - %c (%d) -> %c (%d)", qalSeq[read_index], refSeq[ref_index], trans[refSeq[ref_index]], qrySeq[read_index], trans[qrySeq[read_index]]);
					alignmentPositions->type = baseNumber
							* trans[refSeq[ref_index + k]]
							+ trans[qrySeq[read_index + k]];
					alignmentPositions->readPosition = read_index + k;
					alignmentPositions->refPosition = ref_index + k;
					alignmentPositions->match = false;
//					rates[baseNumber * trans[refSeq[ref_index + k]]
//							+ trans[qrySeq[read_index + k]]]++;
					alignmentPositions++;
				}
			}

			cigar_m_length += length;
			if (!bsMapping && !slamSeq)
				mismatch += length;

			//Produces: 	[0-9]+(([A-Z]+|\^[A-Z]+)[0-9]+)*
			//instead of: 	[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*
			md_offset += sprintf(result.pMD + md_offset, "%d", md_eq_length);
			for (int k = 0; k < length; ++k) {
				if (bsMapping || slamSeq) {
					if (qrySeq[read_index] == bsFrom
							&& refSeq[ref_index] == bsTo) {
						match += 1;
					} else {
						mismatch += 1;
					}
				}
				md_offset += sprintf(result.pMD + md_offset, "%c",
						refSeq[ref_index++]);
				read_index += 1;
			}
			md_eq_length = 0;

			break;
		case CIGAR_EQ:
			if (slamSeq) {
				for (int k = 0; k < length; ++k) {
//					rates[baseNumber * trans[refSeq[ref_index + k]]
//							+ trans[qrySeq[read_index + k]]]++;
					alignmentPositions->type = baseNumber
							* trans[refSeq[ref_index + k]]
							+ trans[qrySeq[read_index + k]];
					alignmentPositions->readPosition = read_index + k;
					alignmentPositions->refPosition = ref_index + k;
					alignmentPositions->match = true;
					alignmentPositions++;
				}
			}
			match += length;
			cigar_m_length += length;
			md_eq_length += length;
			ref_index += length;
			read_index += length;
			break;
		case CIGAR_D:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pCigar + cigar_offset);
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('D', length,
					result.pCigar + cigar_offset);

			md_offset += sprintf(result.pMD + md_offset, "%d", md_eq_length);
			md_eq_length = 0;
			result.pMD[md_offset++] = '^';
			for (int k = 0; k < length; ++k) {
				result.pMD[md_offset++] = refSeq[ref_index++];
			}
			mismatch += length;

			break;
		case CIGAR_I:
			if (cigar_m_length > 0) {
				cigar_offset += printCigarElement('M', cigar_m_length,
						result.pCigar + cigar_offset);
				cigar_m_length = 0;
			}
			cigar_offset += printCigarElement('I', length,
					result.pCigar + cigar_offset);
			read_index += length;
			mismatch += length;
			break;
		default:
			Log.Message("Unable to compute alignment for:");
			Log.Message("Ref: %s", refSeq);
			Log.Message("Qry: %s", qrySeq);
			Log.Warning("This aligment will be discarded. No other alignments will be affected");
			return false;
		}
	}
	md_offset += sprintf(result.pMD + md_offset, "%d", md_eq_length);
	if (cigar_m_length > 0) {
		cigar_offset += printCigarElement('M', cigar_m_length,
				result.pCigar + cigar_offset);
		cigar_m_length = 0;
	}

	if ((gpuCigar[alignment_length - 1] >> 4) > 0) {
		if (Config.GetInt("hard_clip") == 1) { //soft clipping
			cigar_offset += printCigarElement('H',
					gpuCigar[alignment_length - 1] >> 4,
					result.pCigar + cigar_offset);
		} else if (Config.GetInt("silent_clip") != 1) {
			cigar_offset += printCigarElement('S',
					gpuCigar[alignment_length - 1] >> 4,
					result.pCigar + cigar_offset);
		}
		result.QEnd = gpuCigar[alignment_length - 1] >> 4;
	}

//	if ((gpuCigar[alignment_length - 1] >> 4) > 0) {
//		cigar_offset += printCigarElement('S', gpuCigar[alignment_length - 1] >> 4, result.pCigar + cigar_offset);
//		result.QEnd = gpuCigar[alignment_length - 1] >> 4;
//	} else {
//		result.QEnd = 0;
//	}

	result.pCigar[cigar_offset] = '\0';
	result.pMD[md_offset] = '\0';

	result.Identity = match * 1.0f / total;
	result.NM = mismatch;

	//Hack. Use read_index to check whether the alignment is correct. If read_index >= read length: not a valid cigar/md string
	result.Score = (float) read_index;
	return true;
}

#include <iostream>
//TODO: remove

long SWOclCigar::getMaxAllocSize(int const batch_size) {
	//	Log.Message("Batch size:\t %d", batch_size);
	//	Log.Message("Results:\t %d", result_number * batch_size * sizeof(cl_short));
	//
	//	Log.Message("%d %d %d %d", batch_size, (Config.GetInt("corridor") + 2), (Config.GetInt("qry_max_len") + 1), sizeof(cl_char));
	//	Log.Message("Matrix:\t %d", );
	//	Log.Message("Alignments:\t %d\n\n", );

	long corridor = (Config.GetInt("corridor") + 2);
	long qry_len = (Config.GetInt("qry_max_len") + 1);
	long s_char = sizeof(cl_char);

	long matrix = (long) batch_size * corridor * qry_len * s_char;
	//	std::cout << corridor << std::endl;
	//	std::cout << qry_len << std::endl;
	//	std::cout << s_char << std::endl;
	//	std::cout << batch_size << std::endl;
	//
	//	std::cout << std::endl << (corridor * s_char) << std::endl << (corridor * s_char * qry_len) <<
	//			std::endl << (corridor * s_char * qry_len * 92160) << std::endl;
	//	std::cout << "Matrix: " << matrix << std::endl;

	//std::cout << batch_size * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char) << std::endl;
//	batch_size_align * (Config.GetInt("corridor") + 2) * (Config.GetInt("qry_max_len") + 1) * sizeof(cl_char)
	long alignments = (long) batch_size * (long) alignment_length * (long) 2
			* (long) sizeof(cl_char);
	//std::cout << "alignments: " << alignments << std::endl;
	return std::max(matrix, alignments);
}

int SWOclCigar::GetAlignBatchSize() const {
	return batch_size_align * step_count;
}

using std::cout;

int SWOclCigar::computeAlignmentBatchSize() {
	cl_uint mpCount = host->getDeviceInfoInt(CL_DEVICE_MAX_COMPUTE_UNITS);

	if (host->isGPU()) {

//		cout << "MPCount: " << mpCount << " block_multiplier: " << block_multiplier << " Thread per Multi: " << host->getThreadPerMulti()
//				<< " thread_per_block: " << threads_per_block << std::endl;
		int block_count = mpCount * block_multiplier
				* (host->getThreadPerMulti() / threads_per_block);
		block_count = (block_count / mpCount) * mpCount;

		unsigned long largest_alloc = getMaxAllocSize(
				block_count * threads_per_block);

//		Log.Message("Largest alloc: %u %u", largest_alloc, mpCount);
		while (!host->testAllocate(largest_alloc)) {
			block_count -= mpCount;
			largest_alloc = getMaxAllocSize(block_count * threads_per_block);
			Log.Verbose("Reducing batch size to %d", block_count * threads_per_block);
		}

		Log.Verbose("Multi processor count: %d", mpCount);
		Log.Verbose("Max. threads per multi processor: %d", host->getThreadPerMulti());
		Log.Verbose("Threads per block used: %d", threads_per_block);
		Log.Verbose("Block number: %d", block_count);
		Log.Verbose("Batch size: %d", (block_count * threads_per_block));
		//TODO: Print debug info

		return block_count * threads_per_block;
	} else {
		//cout << "mp " << mpCount << std::endl;
		//if (!host->isGPU()) {
		//	block_count *= 4;
		//}
		return 2048;
	}
}
