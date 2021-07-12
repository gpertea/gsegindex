/*
 * iutil.h
 *
 *  Created on: Jul 1, 2021
 *      Author: gpertea
 */

#ifndef IUTIL_H_
#define IUTIL_H_
#include <zlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "khash.h"
#include "kseq.h"
#include "GBase.h"
#include "GResUsage.h"

//Parse a line of BED file
char *parse_bed(char *s, uint32_t *st_, uint32_t *en_);


#endif /* IUTIL_H_ */
