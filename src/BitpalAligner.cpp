#include "BitpalAligner.h"
#include "NGM.h"

/**
 * Copyright (c) 2013, Laboratory for Biocomputing and Informatics, Boston University
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *   + Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *   + Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *   + This source code may not be used in any program that is sold, any
 *    derivative work must allow free distribution and modification.
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE LABORATORY FOR BIOCOMPUTING AND INFORMATICS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//Uses 114 operations, counting |, &, ^, and +
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>

#define wordsize 64

int bitwise_alignment_0_1_1(char * s1, char* s2, int words){

	int debug = 1;
	int wordsizem1 = wordsize - 1;
	int i, j, k, n, m;
	unsigned long long int bitmask;
	unsigned long long int Matches;
	unsigned long long int **matchvec;
	unsigned long long int *matchvecmem;
	unsigned long long int *matchv;
	unsigned long long int NotMatches;
	unsigned long long int all_ones = ~0x0000000000000000;
	unsigned long long int one = 0x0000000000000001;
	unsigned long long int sixtythreebitmask = ~((unsigned long long int) 1 << wordsizem1);

	unsigned long long int DVpos1shift, DVzeroshift, DVneg1shift;
	unsigned long long int DHpos1temp, DHzerotemp, DHneg1temp;
	unsigned long long int  *DHpos1,  *DHzero,  *DHneg1;
	unsigned long long int INITpos1s, INITzeros;
	unsigned long long int INITpos1sprevbit, INITzerosprevbit;
	unsigned long long int OverFlow0, OverFlow1;

	unsigned long long int RemainDHneg1;
	unsigned long long int DHneg1tozero;
	unsigned long long int DVpos1shiftorMatch;
	unsigned long long int DVnot1to1shiftorMatch;
	unsigned long long int DHpos1orMatch;
	unsigned long long int add4;
	unsigned long long int add2;
	unsigned long long int add1;

	unsigned long long int sum;
	unsigned long long int highone = one << (wordsize - 1);
	unsigned long long int initval;
	char * iterate;
	int score;
	int maxscore;
	int counter = 0, w = 0;

	n = strlen (s1);
	m = strlen (s2);
	DHpos1 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHzero = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHneg1 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));



	//*************************encode match strings A C G T N for string1
	//loop through string1 and store bits in matchA, matchC, etc.
	//position zero corresponds to column one in the score matrix, i.e., first character
	//so we start with i = 0 and bitmask = 1
	//below is optimized
	bitmask = 0x0000000000000001;
	matchvec = (unsigned long long int **) calloc(256, sizeof(unsigned long long int *));
	matchvecmem = (unsigned long long int *) calloc(words*256, sizeof(unsigned long long int));
	for (i = 0 ; i < 256; ++i)
	matchvec[i] = &matchvecmem[i*words];
	for (iterate = s1, i = 0; i < n; ++i, ++iterate)
	{
		matchvec[(*iterate)][w] |= bitmask;
		bitmask <<= 1; ++counter;
		if (counter == 63)
		{
			counter = 0;
			w++;
			bitmask = one;
		}
	}



	//intialize top row (penalty for initial gap)
	for (i = 0; i < words; ++i){
		DHzero[i] = all_ones;
		DHneg1[i] = DHpos1[i] = 0;

	}
	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		//initialize OverFlow
		OverFlow0 = OverFlow1 = 0;
		INITpos1sprevbit = INITzerosprevbit = 0;

		matchv = matchvec[*iterate];
		for (j = 0; j < words; ++j){
			DHpos1temp = DHpos1[j];
			DHzerotemp = DHzero[j];
			DHneg1temp = DHneg1[j];

			Matches = *matchv;
			++matchv;
			//Complement Matches
			NotMatches = ~Matches;
			//Finding the vertical values
			//Find 1s
			INITpos1s = DHneg1temp & Matches;
			sum = (INITpos1s + DHneg1temp) + OverFlow0;
			DVpos1shift = ((sum  ^ DHneg1temp) ^ INITpos1s) & sixtythreebitmask;
			OverFlow0 = sum >> wordsizem1;

			//set RemainingDHneg1
			RemainDHneg1 = DHneg1temp ^ (INITpos1s) & sixtythreebitmask;
			//combine 1s and Matches
			DVpos1shiftorMatch = DVpos1shift | Matches;

			//set DVnot1to1shiftorMatch
			DVnot1to1shiftorMatch = ~(DVpos1shiftorMatch);
			//Find 0s
			INITzeros = ((DHzerotemp & DVpos1shiftorMatch) | (DHneg1temp & DVnot1to1shiftorMatch));
			DVzeroshift = (INITzeros << 1) | INITzerosprevbit;
			INITzerosprevbit = (INITzeros & sixtythreebitmask) >> (wordsizem1 - 1);
			//Find -1s
			DVneg1shift = all_ones ^ (DVpos1shift | DVzeroshift);

			//Finding the horizontal values
			//Remove matches from DH values except 1
			//combine 1s and Matches
			DHpos1orMatch = DHpos1temp| Matches;
			//group -1topos0
			DHneg1tozero = all_ones^(DHpos1orMatch);
			//Find 0s
			DHzerotemp = ((DVzeroshift & DHpos1orMatch) | (DVneg1shift & DHneg1tozero));
			//Find 1s
			DHpos1temp = ((DVneg1shift & DHpos1orMatch) );
			//Find -1s
			DHneg1temp = sixtythreebitmask & (all_ones^(DHzerotemp | DHpos1temp));DHpos1[j] = DHpos1temp;
			DHzero[j] = DHzerotemp;
			DHneg1[j] = DHneg1temp;

		}
	}
	//find scores in last row

	score = -1 * m;
	bitmask = one;
	for (j = 0; j < words; ++j){
		DHpos1temp = DHpos1[j];
		DHzerotemp = DHzero[j];
		DHneg1temp = DHneg1[j];
		add1 = DHzerotemp;
		add2 = DHpos1temp;
		add4 = 0;


		for (i = j*wordsizem1; i < (j + 1)*wordsizem1 && i < n; ++i)
		{
			score += (add1 & bitmask) * 1 + (add2 & bitmask) * 2 + (add4 & bitmask) * 4 - 1;
			add1>>= 1;
			add2>>= 1;
			add4>>= 1;

		}
		if (score > maxscore)
		maxscore = score;

	}
	free(DHpos1);
	free(DHzero);
	free(DHneg1);

	free(matchvecmem);
	free(matchvec);
	return score;
}

int bitwise_alignment_multiword_original(char * s1, char* s2, int words){

	int debug = 1;
	int wordsizem1 = wordsize - 1;
	int i, j, k, n, m;
	unsigned long long int bitmask;
	unsigned long long int Matches;
	unsigned long long int **matchvec;
	unsigned long long int *matchvecmem;
	unsigned long long int *matchv;
	unsigned long long int NotMatches;
	unsigned long long int all_ones = ~0x0000000000000000;
	unsigned long long int one = 0x0000000000000001;
	unsigned long long int sixtythreebitmask = ~((unsigned long long int) 1 << wordsizem1);

	unsigned long long int DVpos2shift, DVpos1shift, DVzeroshift, DVneg1shift;
	unsigned long long int DHpos2temp, DHpos1temp, DHzerotemp, DHneg1temp;
	unsigned long long int  *DHpos2,  *DHpos1,  *DHzero,  *DHneg1;
	unsigned long long int INITpos2s, INITpos1s, INITzeros;
	unsigned long long int INITpos2sprevbit, INITpos1sprevbit, INITzerosprevbit;
	unsigned long long int OverFlow0, OverFlow1, OverFlow2;

	unsigned long long int DVpos1shiftNotMatch;

	unsigned long long int RemainDHneg1;
	unsigned long long int DHneg1tozero;
	unsigned long long int DVpos2shiftorMatch;
	unsigned long long int DVnot2to1shiftorMatch;
	unsigned long long int DHpos2orMatch;
	unsigned long long int add4;
	unsigned long long int DHpos1temporDHpos2temp;
	unsigned long long int add2;
	unsigned long long int add1;

	unsigned long long int sum;
	unsigned long long int highone = one << (wordsize - 1);
	unsigned long long int initval;
	char * iterate;
	int score;
	int maxscore;
	int counter = 0, w = 0;

	n = strlen (s1);
	m = strlen (s2);
	DHpos2 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHpos1 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHzero = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHneg1 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));



	//*************************encode match strings A C G T N for string1
	//loop through string1 and store bits in matchA, matchC, etc.
	//position zero corresponds to column one in the score matrix, i.e., first character
	//so we start with i = 0 and bitmask = 1
	//below is optimized
	bitmask = 0x0000000000000001;
	matchvec = (unsigned long long int **) calloc(256, sizeof(unsigned long long int *));
	matchvecmem = (unsigned long long int *) calloc(words*256, sizeof(unsigned long long int));
	for (i = 0 ; i < 256; ++i)
	matchvec[i] = &matchvecmem[i*words];
	for (iterate = s1, i = 0; i < n; ++i, ++iterate)
	{
		matchvec[(*iterate)][w] |= bitmask;
		bitmask <<= 1; ++counter;
		if (counter == 63)
		{
			counter = 0;
			w++;
			bitmask = one;
		}
	}



	//intialize top row (penalty for initial gap)
	for (i = 0; i < words; ++i){
		DHzero[i] = all_ones;
		DHneg1[i] = DHpos1[i] = DHpos2[i] = 0;

	}
	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		//initialize OverFlow
		OverFlow0 = OverFlow1 = OverFlow2 = 0;
		INITpos2sprevbit = INITpos1sprevbit = INITzerosprevbit = 0;

		matchv = matchvec[*iterate];
		for (j = 0; j < words; ++j){
			DHpos2temp = DHpos2[j];
			DHpos1temp = DHpos1[j];
			DHzerotemp = DHzero[j];
			DHneg1temp = DHneg1[j];

			Matches = *matchv;
			++matchv;
			//Complement Matches
			NotMatches = ~Matches;
			//Finding the vertical values
			//Find 2s
			INITpos2s = DHneg1temp & Matches;
			sum = (INITpos2s + DHneg1temp) + OverFlow0;
			DVpos2shift = ((sum  ^ DHneg1temp) ^ INITpos2s) & sixtythreebitmask;
			OverFlow0 = sum >> wordsizem1;

			//set RemainingDHneg1
			RemainDHneg1 = DHneg1temp ^ (INITpos2s) & sixtythreebitmask;
			//combine 2s and Matches
			DVpos2shiftorMatch = DVpos2shift | Matches;

			//Find 1s
			INITpos1s = (DHzerotemp & DVpos2shiftorMatch) ;
			initval = ((INITpos1s << 1) | INITpos1sprevbit);
			INITpos1sprevbit = (initval >> wordsizem1);
			initval &= sixtythreebitmask;
			sum = initval + RemainDHneg1 + OverFlow1;
			DVpos1shift = (sum  ^ RemainDHneg1);
			DVpos1shiftNotMatch = DVpos1shift & NotMatches;
			OverFlow1= sum >> wordsizem1;
			//set DVnot2to1shiftorMatch
			DVnot2to1shiftorMatch = ~(DVpos2shiftorMatch | DVpos1shift);
			//Find 0s
			INITzeros = ((DHpos1temp & DVpos2shiftorMatch) | (DHzerotemp & DVpos1shiftNotMatch)| (DHneg1temp & DVnot2to1shiftorMatch));
			DVzeroshift = (INITzeros << 1) | INITzerosprevbit;
			INITzerosprevbit = (INITzeros & sixtythreebitmask) >> (wordsizem1 - 1);
			//Find -1s
			DVneg1shift = all_ones ^ (DVpos2shift | DVpos1shift | DVzeroshift);

			//Finding the horizontal values
			//Remove matches from DH values except 2
			DHpos1temp &= NotMatches;
			//combine 2s and Matches
			DHpos2orMatch = DHpos2temp| Matches;
			//group -1topos1
			DHneg1tozero = all_ones^(DHpos2orMatch | DHpos1temp);
			//Find 0s
			DHzerotemp = ((DVpos1shift & DHpos2orMatch) | (DVzeroshift & DHpos1temp)| (DVneg1shift & DHneg1tozero));
			//Find 1s
			DHpos1temp = ((DVzeroshift & DHpos2orMatch) | (DVneg1shift & DHpos1temp));
			//Find 2s
			DHpos2temp = ((DVneg1shift & DHpos2orMatch) );
			//Find -1s
			DHneg1temp = sixtythreebitmask & (all_ones^(DHzerotemp | DHpos1temp | DHpos2temp));DHpos2[j] = DHpos2temp;
			DHpos1[j] = DHpos1temp;
			DHzero[j] = DHzerotemp;
			DHneg1[j] = DHneg1temp;

		}
	}
	//find scores in last row

	score = -1 * m;
	bitmask = one;
	for (j = 0; j < words; ++j){
		DHpos2temp = DHpos2[j];
		DHpos1temp = DHpos1[j];
		DHzerotemp = DHzero[j];
		DHneg1temp = DHneg1[j];
		DHpos1temporDHpos2temp = (DHpos1temp|DHpos2temp);
		add1 = DHpos2temp | DHzerotemp;
		add2 = DHpos1temporDHpos2temp;
		add4 = 0;


		for (i = j*wordsizem1; i < (j + 1)*wordsizem1 && i < n; ++i)
		{
			score += (add1 & bitmask) * 1 + (add2 & bitmask) * 2 + (add4 & bitmask) * 4 - 1;
			add1>>= 1;
			add2>>= 1;
			add4>>= 1;

		}
		if (score > maxscore)
		maxscore = score;

	}
	free(DHpos2);
	free(DHpos1);
	free(DHzero);
	free(DHneg1);

	free(matchvecmem);
	free(matchvec);
	return score;
}

int bitwise_alignment_single_word_packed(char * s1, char * s2, int words)
/**
*
* s1 and s2 are the strings to be aligned.
* words is a dummy variable used to give the single word and
* multi word programs the same signature.
*
**/
{
	int debug = 1;
	int i, j, n, m;
	unsigned long long int bitmask;
	unsigned long long int Matches;
	unsigned long long int NotMatches;
	unsigned long long int all_ones = ~0x0000000000000000;
	unsigned long long int one = 0x0000000000000001;

	unsigned long long int DVpos2shift, DVpos1shift, DVzeroshift, DVneg1shift;
	unsigned long long int DHpos2, DHpos1, DHzero, DHneg1;
	unsigned long long int INITpos2s, INITpos1s, INITzeros;
	unsigned long long int DVpos1shiftNotMatch;
	unsigned long long int RemainDHneg1;
	unsigned long long int DHneg1tozero;
	unsigned long long int DVpos2shiftorMatch;
	unsigned long long int DVnot2to1shiftorMatch;
	unsigned long long int DHpos2orMatch;

	unsigned long long int notDVDHbits2;
	unsigned long long int notDVDHbits4;
	unsigned long long int notDVDHbits4DVDHbits2;
	unsigned long long int notDVDHbits4notDVDHbits2;
	unsigned long long int DVDHbits4DVDHbits2;
	unsigned long long int DVDHbits4notDVDHbits2;
	unsigned long long int notDVDHbits1;
	unsigned long long int DVbits4DVbits2;
	unsigned long long int notDVbits4;
	unsigned long long int notDVbits4DVbits2;
	unsigned long long int notDVbits2;
	unsigned long long int notDVbits1;
	unsigned long long int DVbits4notDVbits2;
	unsigned long long int notDVbits4notDVbits2;
	unsigned long long int DHbits1;
	unsigned long long int DVbits1;
	unsigned long long int DVDHbits1;
	unsigned long long int DHbits2;
	unsigned long long int DVbits2;
	unsigned long long int DVDHbits2;
	unsigned long long int DHbits4;
	unsigned long long int DVbits4;
	unsigned long long int DVDHbits4;
	unsigned long long int DVDHbits8;
	unsigned long long int carry0;
	unsigned long long int carry1;
	unsigned long long int carry2;
	unsigned long long int DVDHbithighcomp;
	unsigned long long int compDHneg1tozero;
	unsigned long long int DVpos1shiftorDVpos2shiftorMatch;
	unsigned long long int DHpos2orDHpos1orDHzero;
	unsigned long long int DHpos1orDHzero;

	int maxscore;//Holds the maximum score found in last row
	int score;//Holds the NW score calculated from the last row of DHbit values
	unsigned long long int  * matchvec = (unsigned long long int * )calloc(256, sizeof(unsigned long long int));


	char *iterate;
	n = strlen (s1);
	m = strlen (s2);

	//*************************encode match strings A C G T N for string1
	//loop through string1 and store bits in matchA, matchC, etc.
	//position zero corresponds to column one in the score matrix, i.e., first character
	//so we start with i = 0 and bitmask = 1

	bitmask = 0x0000000000000001;
	for (i = 0, iterate = s1; i < n; ++i, ++iterate)
	{
		matchvec[(int)((*iterate))] |= bitmask;
		bitmask <<= 1;
	}


	//intialize top row (penalty for initial gap, unless semi-global alignment is being done)
	DHzero = all_ones;
	DHneg1 = DHpos1 = DHpos2 = 0;
	DHpos1orDHzero = (DHpos1|DHzero);
	DHpos2orDHpos1orDHzero = (DHpos2|DHpos1orDHzero);
	DHbits1 = DHzero | DHpos2;
	DHbits2 = DHpos1orDHzero;
	DHbits4 = DHpos2orDHpos1orDHzero;

	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		Matches = matchvec[(int)(*iterate)];
		//Complement Matches
		NotMatches = ~Matches;

		//Finding the vertical values.         //Find 2s
		INITpos2s = DHneg1 & Matches;
		DVpos2shift = (((INITpos2s + DHneg1) ^DHneg1) ^ INITpos2s);

		//set RemainingDHneg1
		RemainDHneg1 = DHneg1 ^ (DVpos2shift >> 1);
		//combine 2s and Matches
		DVpos2shiftorMatch = DVpos2shift | Matches;

		//Find 1s
		INITpos1s = (DHzero & DVpos2shiftorMatch) ;
		DVpos1shift = (((INITpos1s << 1) + RemainDHneg1) ^ RemainDHneg1) & NotMatches;
		//set DVnot2to1shiftorMatch
		DVnot2to1shiftorMatch = ~(DVpos2shiftorMatch | DVpos1shift);
		DVpos1shiftorDVpos2shiftorMatch = (DVpos1shift|DVpos2shiftorMatch);
		DVbits1 = DVpos2shiftorMatch | DVnot2to1shiftorMatch;
		DVbits2 = DVpos1shiftorDVpos2shiftorMatch;
		DVbits4 = 0;
		carry0 = (DHbits1 & DVbits1);
		carry1 = ((DHbits2 & DVbits2) | ((DHbits2 ^ DVbits2) & carry0));
		carry2 = ((DHbits4 & DVbits4) | ((DHbits4 ^ DVbits4) & carry1));
		DVDHbits1 = (DHbits1 ^ DVbits1);
		DVDHbits2 = ((DHbits2 ^ DVbits2) ^ carry0);
		DVDHbits4 = ((DHbits4 ^ DVbits4) ^ carry1);
		DVDHbits8 = carry2;
		DVDHbithighcomp = ~DVDHbits4;
		DVDHbits1 &= DVDHbithighcomp;
		DVDHbits2 &= DVDHbithighcomp;
		DVDHbits4 &= DVDHbithighcomp;
		DVDHbits8 &= DVDHbithighcomp;
		DVDHbits1 <<= 1;
		DVDHbits2 <<= 1;
		DVDHbits4 <<= 1;
		DVDHbits8 <<= 1;
		DVbits1 = DVDHbits1;
		DVbits2 = DVDHbits2;
		DVbits4 = DVDHbits4;

		//Finding the horizontal values
		//Remove matches from DH values except 2
		DHneg1tozero = (DHneg1 | DHzero) & NotMatches;
		compDHneg1tozero = all_ones ^ DHneg1tozero;
		//Representing the range -1 to 0 as 0;
		DHbits1 |= DHneg1tozero;
		DHbits2 |= DHneg1tozero;
		DHbits4 |= DHneg1tozero;
		DHbits1 |= Matches;
		DHbits2 &= NotMatches;
		DHbits4 |= Matches;
		carry0 = (DHbits1 & DVDHbits1);
		carry1 = ((DHbits2 & DVDHbits2) | ((DHbits2 ^ DVDHbits2) & carry0));
		carry2 = ((DHbits4 & DVDHbits4) | ((DHbits4 ^ DVDHbits4) & carry1));
		DVDHbits1 = (DHbits1 ^ DVDHbits1);
		DVDHbits2 = ((DHbits2 ^ DVDHbits2) ^ carry0);
		DVDHbits4 = ((DHbits4 ^ DVDHbits4) ^ carry1);
		DVDHbits8 = carry2;
		//Convert everything that is positive to 0
		DVDHbits1 &= DVDHbits4;
		DVDHbits2 &= DVDHbits4;
		DVDHbits4 &= DVDHbits4;
		DVDHbits8 &= DVDHbits4;
		DHbits1 = DVDHbits1;
		DHbits2 = DVDHbits2;
		DHbits4 = DVDHbits4;
		notDVDHbits4 = (~DVDHbits4);
		notDVDHbits2 = ~DVDHbits2;
		DVDHbits4DVDHbits2 = DVDHbits4 & DVDHbits2;
		notDVDHbits4notDVDHbits2 = notDVDHbits4 & notDVDHbits2;
		notDVDHbits1 = ~DVDHbits1;
		DHneg1 = notDVDHbits4notDVDHbits2 & notDVDHbits1;
		DHzero = DVDHbits4DVDHbits2 & DVDHbits1;

	}

	//find scores in last row
	score = -1 * m;

	maxscore = score;

	bitmask = 0x0000000000000001;

	for (i = 0; i < n; i++)
	{
		score -= ((DVDHbits1 & bitmask) >> i) * 1 + ((DVDHbits2 & bitmask) >> i) * 2 - ((DVDHbits4 & bitmask) >> i) * 4 + 1;
		//printf ("%4d", score);

		if (score > maxscore)
		maxscore = score;

		bitmask <<= 1;
	}

	score = maxscore;

	free(matchvec);
	return score;
}


int bitwise_alignment(char * s1, char* s2, int words){

	int debug = 1;
	int wordsizem1 = wordsize - 1;
	int i, j, k, n, m;
	unsigned long long int bitmask;
	unsigned long long int Matches;
	unsigned long long int **matchvec;
	unsigned long long int *matchvecmem;
	unsigned long long int *matchv;
	unsigned long long int NotMatches;
	unsigned long long int all_ones = ~0x0000000000000000;
	unsigned long long int one = 0x0000000000000001;
	unsigned long long int sixtythreebitmask = ~((unsigned long long int) 1 << wordsizem1);

	unsigned long long int DVpos2shift, DVpos1shift, DVzeroshift, DVneg1shift;
	unsigned long long int DHpos2temp, DHpos1temp, DHzerotemp, DHneg1temp;
	unsigned long long int  *DHpos2,  *DHpos1,  *DHzero,  *DHneg1;
	unsigned long long int INITpos2s, INITpos1s, INITzeros;
	unsigned long long int INITpos2sprevbit, INITpos1sprevbit, INITzerosprevbit;
	unsigned long long int OverFlow0, OverFlow1, OverFlow2;

	unsigned long long int DVpos1shiftNotMatch;

	unsigned long long int RemainDHneg1;
	unsigned long long int DHneg1tozero;
	unsigned long long int DVpos2shiftorMatch;
	unsigned long long int DVnot2to1shiftorMatch;
	unsigned long long int DHpos2orMatch;
	unsigned long long int add4;
	unsigned long long int DHpos1temporDHpos2temp;
	unsigned long long int add2;
	unsigned long long int add1;

	unsigned long long int sum;
	unsigned long long int highone = one << (wordsize - 1);
	unsigned long long int initval;
	char * iterate;
	int score;
	int maxscore;
	int counter = 0, w = 0;

	n = strlen (s1);
	m = strlen (s2);
	DHpos2 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHpos1 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHzero = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));
	DHneg1 = (unsigned long long int*)calloc(words, sizeof(unsigned long long int));



	//*************************encode match strings A C G T N for string1
	//loop through string1 and store bits in matchA, matchC, etc.
	//position zero corresponds to column one in the score matrix, i.e., first character
	//so we start with i = 0 and bitmask = 1
	//below is optimized
	bitmask = 0x0000000000000001;
	matchvec = (unsigned long long int **) calloc(256, sizeof(unsigned long long int *));
	matchvecmem = (unsigned long long int *) calloc(words*256, sizeof(unsigned long long int));
	for (i = 0 ; i < 256; ++i)
	matchvec[i] = &matchvecmem[i*words];
	for (iterate = s1, i = 0; i < n; ++i, ++iterate)
	{
		matchvec[(*iterate)][w] |= bitmask;
		bitmask <<= 1; ++counter;
		if (counter == 63)
		{
			counter = 0;
			w++;
			bitmask = one;
		}
	}



	//intialize top row (penalty for initial gap)
	for (i = 0; i < words; ++i){
		DHzero[i] = all_ones;
		DHneg1[i] = DHpos1[i] = DHpos2[i] = 0;

	}
	//recursion
	for (i = 0, iterate = s2; i < m; ++i, ++iterate)
	{
		//initialize OverFlow
		OverFlow0 = OverFlow1 = OverFlow2 = 0;
		INITpos2sprevbit = INITpos1sprevbit = INITzerosprevbit = 0;

		matchv = matchvec[*iterate];
		for (j = 0; j < words; ++j){
			DHpos2temp = DHpos2[j];
			DHpos1temp = DHpos1[j];
			DHzerotemp = DHzero[j];
			DHneg1temp = DHneg1[j];

			Matches = *matchv;
			++matchv;
			//Complement Matches
			NotMatches = ~Matches;
			//Finding the vertical values
			//Find 2s
			INITpos2s = DHneg1temp & Matches;
			sum = (INITpos2s + DHneg1temp) + OverFlow0;
			DVpos2shift = ((sum  ^ DHneg1temp) ^ INITpos2s) & sixtythreebitmask;
			OverFlow0 = sum >> wordsizem1;

			//set RemainingDHneg1
			RemainDHneg1 = DHneg1temp ^ (INITpos2s) & sixtythreebitmask;
			//combine 2s and Matches
			DVpos2shiftorMatch = DVpos2shift | Matches;

			//Find 1s
			INITpos1s = (DHzerotemp & DVpos2shiftorMatch) ;
			initval = ((INITpos1s << 1) | INITpos1sprevbit);
			INITpos1sprevbit = (initval >> wordsizem1);
			initval &= sixtythreebitmask;
			sum = initval + RemainDHneg1 + OverFlow1;
			DVpos1shift = (sum  ^ RemainDHneg1);
			DVpos1shiftNotMatch = DVpos1shift & NotMatches;
			OverFlow1= sum >> wordsizem1;
			//set DVnot2to1shiftorMatch
			DVnot2to1shiftorMatch = ~(DVpos2shiftorMatch | DVpos1shift);
			//Find 0s
			INITzeros = ((DHpos1temp & DVpos2shiftorMatch) | (DHzerotemp & DVpos1shiftNotMatch)| (DHneg1temp & DVnot2to1shiftorMatch));
			DVzeroshift = (INITzeros << 1) | INITzerosprevbit;
			INITzerosprevbit = (INITzeros & sixtythreebitmask) >> (wordsizem1 - 1);
			//Find -1s
			DVneg1shift = all_ones ^ (DVpos2shift | DVpos1shift | DVzeroshift);

			//Finding the horizontal values
			//Remove matches from DH values except 2
			DHpos1temp &= NotMatches;
			//combine 2s and Matches
			DHpos2orMatch = DHpos2temp| Matches;
			//group -1topos1
			DHneg1tozero = all_ones^(DHpos2orMatch | DHpos1temp);
			//Find 0s
			DHzerotemp = ((DVpos1shift & DHpos2orMatch) | (DVzeroshift & DHpos1temp)| (DVneg1shift & DHneg1tozero));
			//Find 1s
			DHpos1temp = ((DVzeroshift & DHpos2orMatch) | (DVneg1shift & DHpos1temp));
			//Find 2s
			DHpos2temp = ((DVneg1shift & DHpos2orMatch) );
			//Find -1s
			DHneg1temp = sixtythreebitmask & (all_ones^(DHzerotemp | DHpos1temp | DHpos2temp));DHpos2[j] = DHpos2temp;
			DHpos1[j] = DHpos1temp;
			DHzero[j] = DHzerotemp;
			DHneg1[j] = DHneg1temp;

		}
	}
	//find scores in last row

	score = -1 * m;
	bitmask = one;
	for (j = 0; j < words; ++j){
		DHpos2temp = DHpos2[j];
		DHpos1temp = DHpos1[j];
		DHzerotemp = DHzero[j];
		DHneg1temp = DHneg1[j];
		DHpos1temporDHpos2temp = (DHpos1temp|DHpos2temp);
		add1 = DHpos2temp | DHzerotemp;
		add2 = DHpos1temporDHpos2temp;
		add4 = 0;


		for (i = j*wordsizem1; i < (j + 1)*wordsizem1 && i < n; ++i)
		{
			score += (add1 & bitmask) * 1 + (add2 & bitmask) * 2 + (add4 & bitmask) * 4 - 1;
			add1>>= 1;
			add2>>= 1;
			add4>>= 1;

		}
		if (score > maxscore)
		maxscore = score;

	}
	free(DHpos2);
	free(DHpos1);
	free(DHzero);
	free(DHneg1);

	free(matchvecmem);
	free(matchvec);
	return score;
}

BitpalAligner::BitpalAligner() {

	Log.Message("Using BitpalAligner");

}

BitpalAligner::~BitpalAligner() {

}

int BitpalAligner::Score(char const * ref, char const * read) {

//	Log.Message("Ref:   %s", ref);
//	Log.Message("Read:  %s", read);
	int result = bitwise_alignment_0_1_1(const_cast<char *>(read), const_cast<char *>(ref), (strlen(ref) / 63 + 1));
//	Log.Message("Score: %d", result);
	return result;

}

int BitpalAligner::BatchScore(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList,
			char const * const * const qalSeqList, float * const results,
			void * extData) {

	for (int i = 0; i < batchSize; ++i) {
		results[i] = Score(refSeqList[i], qrySeqList[i]);
	}

	return batchSize;
}


int BitpalAligner::BatchAlign(int const mode, int const batchSize,
			char const * const * const refSeqList,
			char const * const * const qrySeqList,
			char const * const * const qalSeqList, Align * const results,
			void * extData) {

//	throw "Not implemented yet";

	return batchSize;
}

int BitpalAligner::SingleAlign(int const mode, int const corridor,
		char const * const refSeq, char const * const qrySeq, Align & result,
		void * extData) {

//	throw "Not implemented yet.";

	return 1;
}

