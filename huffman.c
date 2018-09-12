
/*-------------------------------------------------------------*/
/*--- Huffman coding low-level stuff                        ---*/
/*---                                             huffman.c ---*/
/*-------------------------------------------------------------*/

/* ------------------------------------------------------------------
   This file is part of bzip2/libbzip2, a program and library for
   lossless, block-sorting data compression.

   bzip2/libbzip2 version 1.0.6 of 6 September 2010
   Copyright (C) 1996-2010 Julian Seward <jseward@bzip.org>

   Please read the WARNING, DISCLAIMER and PATENTS sections in the 
   README file.

   This program is released under the terms of the license contained
   in the file LICENSE.
   ------------------------------------------------------------------ */


#include "bzlib_private.h"

/*---------------------------------------------------*/
#define WEIGHTOF(zz0)  ((zz0) & 0xffffff00)
#define DEPTHOF(zz1)   ((zz1) & 0x000000ff)
#define MYMAX(zz2,zz3) ((zz2) > (zz3) ? (zz2) : (zz3))

#define ADDWEIGHTS(zw1,zw2)                           \
   (WEIGHTOF(zw1)+WEIGHTOF(zw2)) |                    \
   (1 + MYMAX(DEPTHOF(zw1),DEPTHOF(zw2)))

#define UPHEAP(z)                                     \
{                                                     \
   Int32 zz, tmp;                                     \
   zz = z; tmp = heap[zz];                            \
   while (weight[tmp] < weight[heap[zz >> 1]]) {      \
      heap[zz] = heap[zz >> 1];                       \
      zz >>= 1;                                       \
   }                                                  \
   heap[zz] = tmp;                                    \
}

#define DOWNHEAP(z)                                   \
{                                                     \
   Int32 zz, yy, tmp;                                 \
   zz = z; tmp = heap[zz];                            \
   while (True) {                                     \
      yy = zz << 1;                                   \
      if (yy > nHeap) break;                          \
      if (yy < nHeap &&                               \
          weight[heap[yy+1]] < weight[heap[yy]])      \
         yy++;                                        \
      if (weight[tmp] < weight[heap[yy]]) break;      \
      heap[zz] = heap[yy];                            \
      zz = yy;                                        \
   }                                                  \
   heap[zz] = tmp;                                    \
}


/*---------------------------------------------------*/
void BZ2_hbMakeCodeLengths ( Int32 *len, 
                             Int32 *freq,
                             Int32 alphaSize,
                             Int32 maxLen )
{
   /*--
      Nodes and heap entries run from 1.  Entry 0
      for both the heap and nodes is a sentinel.
   --*/
   Int32 nNodes, nHeap, n1, n2, i, j, k;
   Bool  tooLong;

   Int32 heap   [ BZ_MAX_ALPHA_SIZE + 2 ];
   Int32 weight [ BZ_MAX_ALPHA_SIZE * 2 ];
   Int32 parent [ BZ_MAX_ALPHA_SIZE * 2 ]; 

   for (i = 0; i < alphaSize; i++)
      weight[i+1] = (freq[i] == 0 ? 1 : freq[i]) << 8;

   while (True) {

      nNodes = alphaSize;
      nHeap = 0;

      heap[0] = 0;
      weight[0] = 0;
      parent[0] = -2;

      for (i = 1; i <= alphaSize; i++) {
         parent[i] = -1;
         nHeap++;
         heap[nHeap] = i;
         UPHEAP(nHeap);
      }

      AssertH( nHeap < (BZ_MAX_ALPHA_SIZE+2), 2001 );
   
      while (nHeap > 1) {
         n1 = heap[1]; heap[1] = heap[nHeap]; nHeap--; DOWNHEAP(1);
         n2 = heap[1]; heap[1] = heap[nHeap]; nHeap--; DOWNHEAP(1);
         nNodes++;
         parent[n1] = parent[n2] = nNodes;
         weight[nNodes] = ADDWEIGHTS(weight[n1], weight[n2]);
         parent[nNodes] = -1;
         nHeap++;
         heap[nHeap] = nNodes;
         UPHEAP(nHeap);
      }

      AssertH( nNodes < (BZ_MAX_ALPHA_SIZE * 2), 2002 );

      tooLong = False;
      for (i = 1; i <= alphaSize; i++) {
         j = 0;
         k = i;
         while (parent[k] >= 0) { k = parent[k]; j++; }
         len[i-1] = j;
         if (j > maxLen) tooLong = True;
      }
      
      if (! tooLong) break;

      /* 17 Oct 04: keep-going condition for the following loop used
         to be 'i < alphaSize', which missed the last element,
         theoretically leading to the possibility of the compressor
         looping.  However, this count-scaling step is only needed if
         one of the generated Huffman code words is longer than
         maxLen, which up to and including version 1.0.2 was 20 bits,
         which is extremely unlikely.  In version 1.0.3 maxLen was
         changed to 17 bits, which has minimal effect on compression
         ratio, but does mean this scaling step is used from time to
         time, enough to verify that it works.

         This means that bzip2-1.0.3 and later will only produce
         Huffman codes with a maximum length of 17 bits.  However, in
         order to preserve backwards compatibility with bitstreams
         produced by versions pre-1.0.3, the decompressor must still
         handle lengths of up to 20. */

      for (i = 1; i <= alphaSize; i++) {
         j = weight[i] >> 8;
         j = 1 + (j / 2);
         weight[i] = j << 8;
      }
   }
}


/*---------------------------------------------------*/
void BZ2_hbAssignCodes ( Int32 *code,
                         Int32 *length,
                         Int32 minLen,
                         Int32 maxLen,
                         Int32 alphaSize )
{
   Int32 n, vec, i;

   vec = 0;
   for (n = minLen; n <= maxLen; n++) {
      for (i = 0; i < alphaSize; i++)
         if (length[i] == n) { code[i] = vec; vec++; };
      vec <<= 1;
   }
}


/*---------------------------------------------------*/
#define SET_PERM(idx,val) \
   if (length[idx] < BZ_START_DECODE_LEN) { \
      for (j = base[length[idx]] - 1, v = base[length[idx]] - ((1 << BZ_START_DECODE_LEN) >> length[idx]); j >= v; j--) { \
         perm[j] = (val); \
      } \
      base[length[idx]] = v; \
   } else { \
      perm[--base[length[idx]]] = (val); \
   } \

Int32 BZ2_hbCreateDecodeTables ( Int32 *limit,
                                 Int32 *base,
                                 Int16 *perm,
                                 UChar *length,
                                 Int32 EOB )
{
   Int32 i, j, v;
   Int32 minLen;
   Int32 maxLen;

   for (i = 0; i < BZ_MAX_CODE_LEN; i++) base[i] = 0;

   minLen = BZ_MAX_CODE_LEN; maxLen = 0;
   for (i = EOB; i >= 0; i--) {
      if (length[i] < minLen) minLen = length[i];
      if (length[i] > maxLen) maxLen = length[i];
      base[length[i]]++;
   }

   v = 0;
   for (i = 1; i < BZ_START_DECODE_LEN; i++) {
      v += base[i] << (BZ_START_DECODE_LEN - i);
      base[i] = v;
   }
   for (i = BZ_START_DECODE_LEN; i < BZ_MAX_CODE_LEN; i++) {
      v += base[i];
      base[i] = v;
   }

   if (base[BZ_START_DECODE_LEN] > (1 << BZ_START_DECODE_LEN) ||
       (base[BZ_START_DECODE_LEN] < (1 << BZ_START_DECODE_LEN) && maxLen <= BZ_START_DECODE_LEN)) return -1;

   SET_PERM(EOB, 0);
   for (i = EOB - 1; i > BZ_RUNB ; i--) {
      SET_PERM(i, i - 1);
   }
   SET_PERM(BZ_RUNB, 2 << 8);
   SET_PERM(BZ_RUNA, 1 << 8);

   for (i = minLen, j = 0; i <= BZ_START_DECODE_LEN; i++) {
      for (v = base[i + 1]; j < v; j++) length[j] = BZ_START_DECODE_LEN - i;
   }

   v = base[BZ_START_DECODE_LEN + 1];
   limit[0] = v - 1;
   if (maxLen > BZ_START_DECODE_LEN) {
      for (i = BZ_START_DECODE_LEN + 1; i <= maxLen; i++) {
         v = (v << 1) - base[i];
         base[i] = v;
         v += base[i + 1];
         limit[i] = v - 1;
      }
      if (v != (1 << maxLen)) return -1;
   }

   return BZ_START_DECODE_LEN;
}


/*-------------------------------------------------------------*/
/*--- end                                         huffman.c ---*/
/*-------------------------------------------------------------*/
