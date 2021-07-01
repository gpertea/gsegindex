/*
 * iutil.h
 *
 *  Created on: Jul 1, 2021
 *      Author: gpertea
 */

#ifndef IUTIL_H_
#define IUTIL_H_

#include "khash.h"
#include "kseq.h"

//Parse a line of BED file
char *parse_bed(char *s, int32_t *st_, int32_t *en_);



#endif /* IUTIL_H_ */
