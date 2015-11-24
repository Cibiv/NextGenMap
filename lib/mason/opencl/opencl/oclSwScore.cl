#ifdef __CPU__

#ifndef __ALT_SCORING__
__kernel void oclSW_Score(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix) {
#else
__kernel void oclSW_Score(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix, __global char const * direction) {

	int4 direction4 =  convert_int4(vload4(global_index, direction));
#endif

	matrix = matrix + (global_index * ((corridor_length + 1) * (read_length + 1)) * 4);
	read = read + (global_index * read_length * 4);
	scaff = scaff + (global_index * ref_length * 4);
	result = result + result_number * global_index * 4;

	float4 best_read_index4 = null4;
	float4 best_ref_index4 = null4;
	float4 qend4 = null4;

//	for(int k = 0; k < 4; ++k) {
//		printf("Ref  (GPU) %d:", k);
//		for(int i = 0; i < read_length + corridor_length - 1; ++i) {
//			printf("%c", scaff[i + k * ref_length]);
//		}
//		printf("\n");
//		printf("Read (GPU) %d:", k);
//		for(int i = 0; i < read_length; ++i) {
//			if(read[i + k * read_length] == '\0') {
//				printf("X");
//			} else {
//				printf("%c", read[i + k * read_length]);
//			}
//		}
//		printf("\n");
//	}

	if (read[0] != line_end) {
		float4 l_matrix_lines[MATRIX_LENGTH];
		float4 * local_matrix_line = l_matrix_lines + local_index * corridor_length;

		//Init matrix lines
		for (short i = 0; i < corridor_length; ++i) {
			local_matrix_line[i] = null4;
			vstore4(CIGAR_STOP, i, matrix);
		}
		vstore4(CIGAR_STOP, corridor_length, matrix);

		float4 curr_max = -1;
		float4 read_index = null4;

		for (short read_pos = 0; read_pos < read_length; ++read_pos) {
			int4 read_char_cache = (int4)(trans[read[read_pos + 0 * read_length]], trans[read[read_pos + 1 * read_length]], trans[read[read_pos + 2 * read_length]], trans[read[read_pos + 3 * read_length]]);
			//			printf("%d %d %d %d\n", test.s0, test.s1, test.s2, test.s3);
			matrix += (corridor_length + 1) * 4;
			float4 left_cell = null4;
			vstore4(CIGAR_STOP, 0, matrix);

			for (short ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {
				//init values
				left_cell += gap_ref4;
#ifndef __ALT_SCORING__
				float4 score = (float4)(scores[read_char_cache.s0][trans[scaff[ref_index]]], scores[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scores[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scores[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]);
				float4 diag_cell = local_matrix_line[ref_index] + score;
				float4 pointerX = select(CIGAR_X4, CIGAR_EQ4, (score == match));
#else
				float4 score = select((float4)(scoresREV[read_char_cache.s0][trans[scaff[ref_index]]], scoresREV[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scoresREV[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scoresREV[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]), (float4)(scoresFWD[read_char_cache.s0][trans[scaff[ref_index]]], scoresFWD[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scoresFWD[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scoresFWD[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]), direction4 == 0);

				float4 diag_cell = local_matrix_line[ref_index] + score;
				float4 pointerX = select(CIGAR_X4, CIGAR_EQ4, (read_char_cache == (int4)(trans[scaff[ref_index]], trans[scaff[ref_index + 1 * ref_length]], trans[scaff[ref_index + 2 * ref_length]], trans[scaff[ref_index + 3 * ref_length]])));
#endif

				float4 up_cell = local_matrix_line[(ref_index + 1)] + gap_read4;

				//find max
				float4 max_cell = null4;
				max_cell = max(left_cell, max_cell);
				max_cell = max(diag_cell, max_cell);
				max_cell = max(up_cell, max_cell);

				float4 pointer = 0;
				pointer = select(pointer, CIGAR_D, max_cell == left_cell);
				pointer = select(pointer, CIGAR_I, max_cell == up_cell);
				pointer = select(pointer, pointerX, max_cell == diag_cell || max_cell == (local_matrix_line[ref_index] + mismatch));
				pointer = select(pointer, CIGAR_STOP, pointer == 0 || max_cell <= 0);

				vstore4(convert_char4(pointer), (ref_index + 1), matrix);

				best_read_index4 = select(best_read_index4, read_index, (max_cell > curr_max));
				best_ref_index4 = select(best_ref_index4, ref_index, (max_cell > curr_max));

				curr_max = max(curr_max, max_cell);

				left_cell = max_cell;

				local_matrix_line[ref_index] = max_cell;
			}
			vstore4(CIGAR_STOP, corridor_length, matrix);
			scaff += 1;

			read_index = select(read_index + 1, read_index, read_char_cache == (int4) 6);
		}
		qend4 = read_index - best_read_index4 - 1;
	}
	vstore4(convert_short4(best_read_index4), param_best_read_index, result);
	vstore4(convert_short4(best_ref_index4), param_best_ref_index, result);
	vstore4(convert_short4(qend4), qend, result);
}


#ifndef __ALT_SCORING__
__kernel void oclSW(__global char const * scaff, __global char const * read, __global float * results) {
#else
__kernel void oclSW(__global char const * scaff, __global char const * read, __global float * results, __global char const * direction) {

	int4 direction4 =  convert_int4(vload4(global_index, direction));
#endif
	scaff = scaff + (global_index * ref_length * 4);
	read = read + (global_index * read_length * 4);

	float4 l_matrix_lines[MATRIX_LENGTH];
	float4 * matrix_lines = l_matrix_lines + local_index * corridor_length;

	float4 curr_max = (float4)(-1.0f);
	if (read[0] != line_end) {
		//printf("Ref:  %s\n", scaff);
		//printf("Read: %s\n", read);
		//Init matrix lines
		for (short i = 0; i < corridor_length; ++i) {
			matrix_lines[i] = null4;
		}
		for (int read_pos = 0; read_pos < read_length; ++read_pos) {

			int4 read_char_cache = (int4)(trans[read[read_pos]], trans[read[read_pos + 1 * read_length]], trans[read[read_pos + 2 * read_length]], trans[read[read_pos + 3 * read_length]]);
			float4 left_cell = 0;
			for (int ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {

				//Max can be set to left_cell because we are computing the score row-wise, so the actual max is the left_cell for the next cell
				left_cell = max(null4, left_cell + gap_ref);
				left_cell = max(matrix_lines[(ref_index + 1)] + gap_read, left_cell);
#ifndef __ALT_SCORING__
				float4 scores4 = (float4)(scores[read_char_cache.s0][trans[scaff[ref_index]]], scores[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scores[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scores[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]);
#else
				float4 scores4 = select((float4)(scoresREV[read_char_cache.s0][trans[scaff[ref_index]]], scoresREV[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scoresREV[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scoresREV[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]), (float4)(scoresFWD[read_char_cache.s0][trans[scaff[ref_index]]], scoresFWD[read_char_cache.s1][trans[scaff[ref_index + 1 * ref_length]]], scoresFWD[read_char_cache.s2][trans[scaff[ref_index + 2 * ref_length]]], scoresFWD[read_char_cache.s3][trans[scaff[ref_index + 3 * ref_length]]]), direction4 == 0);
#endif
				left_cell = max(matrix_lines[ref_index] + scores4, left_cell);
				curr_max = max(curr_max, left_cell);
				matrix_lines[ref_index] = left_cell;
			}
			scaff += 1;
		}
	}
	vstore4(curr_max, global_index, results);
	//printf("%f %f %f %f\n", curr_max.s0, curr_max.s1, curr_max.s2, curr_max.s3);
}

__kernel void oclSW_wSSE(__global char const * scaff, __global char const * read, __global float * results) {

	scaff = scaff + (global_index * ref_length);
	read = read + (global_index * read_length);

	float l_matrix_lines[MATRIX_LENGTH];
	float * matrix_lines = l_matrix_lines + local_index * corridor_length;

	float curr_max = -1.0f;
	if (read[0] != line_end) {
		//Init matrix lines
		for (short i = 0; i < corridor_length; ++i) {
			matrix_lines[i] = 0.0f;
		}
		for (int read_pos = 0; read_pos < read_length; ++read_pos) {
			int read_char_cache = trans[read[read_pos]];
			float left_cell = 0.0f;
			for (int ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {

				//Max can be set to left_cell because we are computing the score row-wise, so the actual max is the left_cell for the next cell
				left_cell = max(0.0f, left_cell + gap_ref);
				left_cell = max(matrix_lines[(ref_index + 1)] + gap_read, left_cell);

				float score = scores[read_char_cache][trans[scaff[ref_index]]];

				left_cell = max(matrix_lines[ref_index] + score, left_cell);
				curr_max = max(curr_max, left_cell);

				matrix_lines[ref_index] = left_cell;
			}
			scaff += 1;
		}
	}
	results[global_index] = curr_max;
}

#endif

__kernel
void interleaveSeq(__global char const * _src, __global char * dest) {

	__global char const * src = _src + ((global_index) * ref_length);
	dest = dest + ((global_index / interleave_number) * interleave_number * ref_length) + global_index % interleave_number;

	int i = 0;
	for (i = 0; i < ref_length / 4 * 4; i += 4) {
		char4 temp = vload4(0, src);
		dest[(i) * interleave_number] = temp.s0;
		dest[(i + 1) * interleave_number] = temp.s1;
		dest[(i + 2) * interleave_number] = temp.s2;
		dest[(i + 3) * interleave_number] = temp.s3;
		src += 4;
	}

	for (; i < ref_length; i += 1) {
		dest[(i) * interleave_number] = *src++;
	}
}


#ifdef __GPU__

#ifndef __ALT_SCORING__
__kernel void oclSW_Score(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix) {
#else
__kernel void oclSW_Score(__global char const * scaff, __global char const * read, __global short * result, __global char * matrix, __global char const * direction) {
#endif

	matrix = matrix + ((global_index / threads_per_block) * threads_per_block * ((corridor_length + 1) * (read_length + 1))) + global_index % threads_per_block;
	read = read + (global_index * read_length);
	scaff = scaff + ((global_index / interleave_number) * interleave_number * ref_length) + global_index % interleave_number;
	result = result + result_number * global_index;

	__local
	short l_matrix_lines[MATRIX_LENGTH];
	__local
	short * local_matrix_line = l_matrix_lines + local_index;

	//Init matrix lines
	for (short i = 0; i < corridor_length; ++i) {
		local_matrix_line[i * threads_per_block] = 0;
		matrix[i * threads_per_block] = CIGAR_STOP;
	}
	matrix[corridor_length * threads_per_block] = CIGAR_STOP;

	short curr_max = -1;
	short read_index = 0;

	//	if (read[0] != '\0') {
	//		//printf("GPU Ref : %s\n", scaff);
	//		//printf("GPU Read: %s\n", read);
	//		printf("S:\t");
	//		char c;
	//		int i = 0;
	//		while ((c = scaff[i * interleave_number]) != '\0') {
	//			printf("\t%c", c);
	//			i += 1;
	//		}
	//		printf("\n");
	//	}

	for (char read_char_cache; (read_char_cache = *read) != line_end; ++read) {
		//						printf("%c:\t", read_char_cache);
		//						for (short i = 0; i < read_index; ++i) {
		//							printf("\t");
		//						}
		matrix += (corridor_length + 1) * threads_per_block;
		short left_cell = 0;
		matrix[0] = CIGAR_STOP;
		for (short ref_index = 0; ref_index < corridor_length - 1; ++ref_index) {

			//init values
			left_cell += gap_ref;
			short diag_cell = local_matrix_line[ref_index * threads_per_block];
			//			printf("%c == %c\n", read_char_cache, scaff[ref_index * interleave_number]);
			int pointer = CIGAR_X;

#ifndef __ALT_SCORING__
			diag_cell += scores[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
			pointer = select(CIGAR_X, CIGAR_EQ, (read_char_cache == scaff[ref_index * interleave_number]));
#else
			short score = 0;
			if(direction[global_index] == 0) {
				score = scoresFWD[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
			} else {
				score = scoresREV[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
			}
			diag_cell += score;
			pointer = select(CIGAR_X, CIGAR_EQ, (read_char_cache == scaff[ref_index * interleave_number]));
#endif
			short up_cell = local_matrix_line[(ref_index + 1) * threads_per_block] + gap_read;

			//find max
			short max_cell = 0;
			max_cell = max(left_cell, max_cell);
			max_cell = max(diag_cell, max_cell);
			max_cell = max(up_cell, max_cell);

			if(max_cell <= 0) {
				pointer = CIGAR_STOP;
			} else if(max_cell == diag_cell ||
		 	   max_cell == (local_matrix_line[ref_index * threads_per_block] + mismatch)) {
			} else if(max_cell == up_cell) {
				pointer = CIGAR_I;
			} else {
				pointer = CIGAR_D;
			}

			//printf("\t%d", max_cell);
			matrix[(ref_index + 1) * threads_per_block] = pointer;

			if (max_cell > curr_max) {
				curr_max = max_cell;
				result[param_best_read_index] = read_index;
				result[param_best_ref_index] = ref_index;
			}

			left_cell = max_cell;
			local_matrix_line[ref_index * threads_per_block] = max_cell;
		}
		matrix[corridor_length * threads_per_block] = CIGAR_STOP;
		scaff += interleave_number;
		read_index += 1;

		//for (short i = 0; i < corridor_length + 1; ++i) {
		//	printf("%d ", matrix[i]);
		//}
		//printf("\n");
	}
	result[qend] = read_index - result[param_best_read_index] - 1;
	if (read_index == 0) {
		result[param_best_read_index] = result[param_best_read_index] = 0;
	}
}

#ifndef __ALT_SCORING__
__kernel void oclSW(__global char const * scaff, __global char const * read, __global float * results) {
#else
__kernel void oclSW(__global char const * scaff, __global char const * read, __global float * results, __global char const * direction) {
#endif

	read = read + (global_index * read_length);
	scaff = scaff + ((global_index / interleave_number) * interleave_number * ref_length) + (global_index % interleave_number);

	__local
	short l_matrix_lines[MATRIX_LENGTH];
	__local
	short * matrix_lines = l_matrix_lines + local_index;
	//Init matrix lines
	for (short i = 0; i < corridor_length; ++i) {
		matrix_lines[i * threads_per_block] = 0;
	}

	short curr_max = -1;
	for (int read_pos = 0; read_pos < read_length; ++read_pos) {

		char read_char_cache = *read++;
		short left_cell = 0;
		for (int ref_index = 0; ref_index < (corridor_length - 1); ++ref_index) {

			short diag_cell = matrix_lines[ref_index * threads_per_block];
#ifndef __ALT_SCORING__
			diag_cell += scores[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
#else
			if(direction[global_index] == 0) {
				diag_cell += scoresFWD[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
			} else {
				diag_cell += scoresREV[trans[read_char_cache]][trans[scaff[ref_index * interleave_number]]];
			}
#endif
			//Max can be set to left_cell because we are computing the score row-wise, so the actual max is the left_cell for the next cell
			left_cell = max(null, (short) (left_cell + gap_ref));
			left_cell = max((short) (matrix_lines[ref_index * threads_per_block + threads_per_block] + gap_read), left_cell);
			left_cell = max(diag_cell, left_cell);

			curr_max = max(curr_max, left_cell);
			matrix_lines[ref_index * threads_per_block] = left_cell;
		}
		scaff += interleave_number;
	}
	results[global_index] = curr_max;
}

#endif
