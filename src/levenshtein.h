#ifndef LEVENSHTEIN_H
#define LEVENSHTEIN_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LEVLENGTH 50

unsigned int ndiff(const char *a, const char *b);

/* `levenshtein.h` - levenshtein
 * MIT licensed.
 * Copyright (c) 2015 Titus Wormer <tituswormer@gmail.com> */

/* Returns an unsigned integer, depicting
 * the difference between `a` and `b`.
 * See http://en.wikipedia.org/wiki/Levenshtein_distance
 * for more information.
 */

unsigned int levenshtein(const char *a, const char *b);

bool multiple_levenshtein(const char *a, const char *b, unsigned int len, unsigned int min_t);
bool multiple_ndiff(const char *a, const char *b, unsigned int len, unsigned int min_t);
bool multiple_comp(const char *a, const char *b, unsigned int len, unsigned int min_t,
				   unsigned int (*f)(const char *, const char *));

/* `levenshtein.c` - levenshtein
 * MIT licensed.
 * Copyright (c) 2015 Titus Wormer <tituswormer@gmail.com> */

/* Returns an unsigned integer, depicting
 * the difference between `a` and `b`.
 * See http://en.wikipedia.org/wiki/Levenshtein_distance
 * for more information. */

unsigned int ndiff(const char *a, const char *b){
    unsigned int length = strlen(a);
	unsigned int index = 0;
	unsigned int n = 0;
	
	while(index < length){
		if(a[index] != b[index]) n++;
		index++;
	}

	return(n);
}

unsigned int levenshtein(const char *a, const char *b) {
    unsigned int length = strlen(a);
    unsigned int bLength = strlen(b);
    unsigned int cache[MAX_LEVLENGTH];
	memset(&cache[0],0,sizeof(unsigned int)*MAX_LEVLENGTH);
	/* unsigned int *cache = (unsigned int *)calloc(length, sizeof(unsigned int)); */
	unsigned int index = 0;
    unsigned int bIndex = 0;
    unsigned int distance;
    unsigned int bDistance;
    unsigned int result;
    char code;

    /* Shortcut optimizations / degenerate cases. */
    if (a == b) {
        return 0;
    }

    if (length == 0) {
        return bLength;
    }

    if (bLength == 0) {
        return length;
    }

    /* initialize the vector. */
    while (index < length) {
        cache[index] = index + 1;
        index++;
    }

    /* Loop. */
    while (bIndex < bLength) {
        code = b[bIndex];
        result = distance = bIndex++;
        index = -1;

        while (++index < length) {
            bDistance = code == a[index] ? distance : distance + 1;
            distance = cache[index];

      cache[index] = result = distance > result
        ? bDistance > result
          ? result + 1
          : bDistance
          : bDistance > distance
          ? distance + 1
          : bDistance;
        }
    }

	/* free(cache); */
	
    return result;
}

bool multiple_levenshtein(const char *a, const char *b, unsigned int len, unsigned int min_t){
	return multiple_comp(a,b,len,min_t,levenshtein);
}

bool multiple_ndiff(const char *a, const char *b, unsigned int len, unsigned int min_t){
	return multiple_comp(a,b,len,min_t,ndiff);
}

bool multiple_comp(const char *a, const char *b, unsigned int len, unsigned int min_t,
				   unsigned int (*f)(const char *, const char *)){
    char buffer[MAX_LEVLENGTH];
	unsigned int min_l = 1e9;
    for(unsigned int i=0;i<(unsigned int)(strlen(b));i+=len){
        strncpy(&buffer[0],&b[i],len); buffer[len]='\0';
		unsigned int r = f(a,(const char *)buffer);
		if(r<min_l){ min_l=r; if(min_t!=0 && min_l<min_t){ return false; } }
    }
    return true;
}

#endif // LEVENSHTEIN_H
