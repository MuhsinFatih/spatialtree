#pragma once

#include <iostream>
#include <stdio.h>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <random>
#include <algorithm>
#include <chrono>
#include <functional>
using namespace std;


extern "C" {
	#include "spatial_tree.h"
}


#define REP(size) for(size_t i=0; i<size; ++i)
#define REPJ(size) for(size_t j=0; j<size; ++j)
#define REPW(size) size_t w=0; while(w<size)
#define vi vector<int>
#define vd vector<double_t>
#define vs vector<string>
#define st size_t
#define vul vector<unsigned long>
#define ul unsigned long

template <typename T>
void printArray(T arr, const st size, const char* heading = "") {
	cout << heading;
	st c = size;

	REP(c) printf(i<c-1 ? "%.3g, " : "%.3g\n", (double)(arr[i]));
}

template <typename T>
void printArray_2d(T arr, const st s1, const st s2, const char* heading = "") {
	cout << heading;

	REP(s1) {
		if(i) printf("\n");
		printf("[");
		REPJ(s2) printf(j<s2-1 ? "%.3g, " : "%.3g\n", (double)(arr[i][j]));
		printf("]");
	}
}

/**
 * fast 2^x
 */
template<typename T>
T _2_to_power(T x) {
	return 1 << x;
}
// returns the position of the first bit that is set to 1. Index starts from 0
#define firstsetbit(u) (size_t)(log2((u)&(-(u))))
#define firstNOTsetbit(u) firstsetbit(~u)
#define firstbitindex(u, b) (b ? firstsetbit(~(u)) : firstsetbit(u))
#define lastnbits(u, n) (u&((1ull<<n)-1))
#define reverse_lastnbits(u, n) (u^((1ull<<n)-1))
template<typename T>
string _dumpbits(T *n) {
	char* nc = (char*)n;
	string res = "";
	size_t s = sizeof(T);
	bool empty = true;
	for(size_t i=0; i<s; ++i){ // notice it's little endian so this is actually reverse. It will start from right
		char c = *(nc+i);
		if(c) empty = false;
		for(size_t j=0; j<sizeof(char)*8; ++j) {
			// if(!(c >> j)) break;
			res = (((c >> j) & 1) ? "1" : "0") + res;
		}
	}
	if(empty) return "";
	return res;
}

template<typename T>
string dumpbits(T *n, string delim) {
	string res = _dumpbits(n);
	if(res.empty()) res = "0";
	std::cout << res << delim;
	return res;
}

template<typename T>
string dumpbits(T *n, size_t min_size, string delim) {
	string res = _dumpbits(n);
	if(res.size() < min_size)
		res = string(min_size - res.size(), '0') + res;
	std::cout << res << delim;
	return res;
}

template<typename T>
string getbits(T n, size_t min_size) {
	T* num = &n;
	string res = _dumpbits(num);
	if(res.size() < min_size)
		res = string(min_size - res.size(), '0') + res;
	return res;
}
template<typename T>
string getbits(T n) {
	T* num = &n;
	string res = _dumpbits(num);
	if(res.empty()) res = "0";
	return res;
}

void dumpbytes(unsigned char* charPtr, size_t size) {
	for(int i=size-1; i>=0; --i) {
		printf("%02x ", charPtr[i]);
	}
	printf("\n");
}