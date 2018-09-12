
/*-------------------------------------------------------------*/
/*--- Decompression machinery                               ---*/
/*---                                          decompress.c ---*/
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
#define RETURN(rrr)                               \
   { goto save_state_and_return; };

#define RETERR(rrr)                               \
   { retVal = rrr; goto error_return; };

#define GET_BITS(lll,vvv,nnn)                     \
   while (True) {                                 \
      if (bsLive >= nnn) {                        \
         UInt32 v;                                \
         v = (bsBuff >>                           \
             (bsLive-nnn)) & ((1 << nnn)-1);      \
         bsLive -= nnn;                           \
         vvv = v;                                 \
         break;                                   \
      }                                           \
      case lll:                                   \
      if (bsPtr == bsEnd) {                       \
         s->state = lll;                          \
         RETURN(BZ_OK);                           \
      }                                           \
      bsBuff                                      \
         = (bsBuff << 8) | (UInt32)(*bsPtr);      \
      bsLive += 8;                                \
      bsPtr++;                                    \
   }

/* reduced version of GET_BITS when nnn <= 8 */
#define GET_BITS_SMALL(lll,vvv,nnn)               \
   do {                                           \
      if ((bsLive -= nnn) < 0) {                  \
         case lll:                                \
         if (bsPtr == bsEnd) {                    \
            s->state = lll;                       \
            RETURN(BZ_OK);                        \
         }                                        \
         bsBuff                                   \
            = (bsBuff << 8) | (UInt32)(*bsPtr);   \
         bsLive += 8;                             \
         bsPtr++;                                 \
      }                                           \
      {                                           \
         UInt32 v;                                \
         v = (bsBuff >>                           \
             (bsLive)) & ((1 << nnn)-1);          \
         vvv = v;                                 \
      }                                           \
   } while(False);

/* reduced version of GET_BITS when nnn == 8 */
#define GET_UCHAR(lll,vvv)                        \
   do {                                           \
      case lll:                                   \
      if (bsPtr == bsEnd) {                       \
         s->state = lll;                          \
         RETURN(BZ_OK);                           \
      }                                           \
      bsBuff                                      \
         = (bsBuff << 8) | (UInt32)(*bsPtr);      \
      bsPtr++;                                    \
      {                                           \
         UInt32 v;                                \
         v = (bsBuff >> bsLive) & ((1 << 8) - 1); \
         vvv = v;                                 \
         break;                                   \
      }                                           \
   } while(False);

#define GET_BIT(lll,uuu)                          \
   GET_BITS_SMALL(lll,uuu,1)

#define UNGET_BITS(nnn)                           \
   bsLive += (nnn);

/*---------------------------------------------------*/
/* returns 0 on EOB, 256 on BZ_RUNA, 512 on BZ_RUNB, */
/* mtf value of symbol on 1 to EOB-1,                */
#define GET_MTF_VAL(label1,label2,lval)           \
{                                                 \
   if ((groupPos -= 1) < 0) {                     \
      if (--groupNo < 0)                          \
         RETERR(BZ_DATA_ERROR);                   \
      gSel = s->selector[groupNo];                \
      gRec = &(s->groupRec[gSel]);                \
      groupPos = BZ_G_SIZE - 1;                   \
   }                                              \
   do {                                           \
      if ((bsLive -= BZ_START_DECODE_LEN) < 0) {  \
         while (bsPtr == bsEnd) {                 \
            s->state = label1;                    \
            RETURN(BZ_OK);                        \
            case label1:                          \
            gSel = s->selector[groupNo];          \
            gRec = &(s->groupRec[gSel]);          \
         }                                        \
         bsBuff                                   \
            = (bsBuff << 8) | (UInt32)(*bsPtr);   \
         bsLive += 8;                             \
         bsPtr++;                                 \
      }                                           \
      zvec = (bsBuff >> bsLive) &                 \
             ((1 << BZ_START_DECODE_LEN) - 1);    \
   } while(False);                                \
   if ((Int32)zvec <= gRec->limit[0]) {           \
      UNGET_BITS(gRec->len[zvec]);                \
      lval = gRec->perm[zvec];                    \
   } else {                                       \
      zn = BZ_START_DECODE_LEN;                   \
      do {                                        \
         zn++;                                    \
         if (--bsLive < 0) {                      \
            while (bsPtr == bsEnd) {              \
               s->save_zn = zn;                   \
               s->save_zvec = zvec;               \
               s->state = label2;                 \
               RETURN(BZ_OK);                     \
               case label2:                       \
               gSel = s->selector[groupNo];       \
               gRec = &(s->groupRec[gSel]);       \
               zn = s->save_zn;                   \
               zvec = s->save_zvec;               \
            }                                     \
            bsBuff = (UInt32)(*bsPtr);            \
            bsLive = 7;                           \
            bsPtr++;                              \
            zj = (Int)bsBuff >> 7;                \
         } else {                                 \
            zj = ((Int)bsBuff >> bsLive) & 1;     \
         }                                        \
         zvec = (zvec << 1) + zj;                 \
      } while ((Int32)zvec > gRec->limit[zn]);    \
      lval = gRec->perm[zvec - gRec->base[zn]];   \
   };                                             \
}

/*---------------------------------------------------*/
/* additionally saves / restores es and N before /   */
/* after RETURN                                      */
#define GET_MTF_VAL_ES(label1,label2,lval)        \
{                                                 \
   if ((groupPos -= 1) < 0) {                     \
      if (--groupNo < 0)                          \
         RETERR(BZ_DATA_ERROR);                   \
      gSel = s->selector[groupNo];                \
      gRec = &(s->groupRec[gSel]);                \
      groupPos = BZ_G_SIZE - 1;                   \
   }                                              \
   do {                                           \
      if ((bsLive -= BZ_START_DECODE_LEN) < 0) {  \
         while (bsPtr == bsEnd) {                 \
            s->save_es = es;                      \
            s->save_N = N;                        \
            s->state = label1;                    \
            RETURN(BZ_OK);                        \
            case label1:                          \
            gSel = s->selector[groupNo];          \
            gRec = &(s->groupRec[gSel]);          \
            es = s->save_es;                      \
            N = s->save_N;                        \
         }                                        \
         bsBuff                                   \
            = (bsBuff << 8) | (UInt32)(*bsPtr);   \
         bsLive += 8;                             \
         bsPtr++;                                 \
      }                                           \
      zvec = (bsBuff >> bsLive) &                 \
             ((1 << BZ_START_DECODE_LEN) - 1);    \
   } while(False);                                \
   if ((Int32)zvec <= gRec->limit[0]) {           \
      UNGET_BITS(gRec->len[zvec]);                \
      lval = gRec->perm[zvec];                    \
   } else {                                       \
      zn = BZ_START_DECODE_LEN;                   \
      do {                                        \
         zn++;                                    \
         if (--bsLive < 0) {                      \
            while (bsPtr == bsEnd) {              \
               s->save_es = es;                   \
               s->save_N = N;                     \
               s->save_zn = zn;                   \
               s->save_zvec = zvec;               \
               s->state = label2;                 \
               RETURN(BZ_OK);                     \
               case label2:                       \
               gSel = s->selector[groupNo];       \
               gRec = &(s->groupRec[gSel]);       \
               es = s->save_es;                   \
               N = s->save_N;                     \
               zn = s->save_zn;                   \
               zvec = s->save_zvec;               \
            }                                     \
            bsBuff = (UInt32)(*bsPtr);            \
            bsLive = 7;                           \
            bsPtr++;                              \
            zj = (Int)bsBuff >> 7;                \
         } else {                                 \
            zj = ((Int)bsBuff >> bsLive) & 1;     \
         }                                        \
         zvec = (zvec << 1) + zj;                 \
      } while ((Int32)zvec > gRec->limit[zn]);    \
      lval = gRec->perm[zvec - gRec->base[zn]];   \
   };                                             \
}


/*---------------------------------------------------*/
Int32 BZ2_decompress ( DState* s )
{
   UInt32     uc;
   Int32      retVal;
   UChar      *bsPtr, *bsEnd;
   UInt32     bsBuff;
   Int32      bsLive;
   Int        gSel;
   Int        zj;
   UInt       nextSym;
   GRec*      gRec;

   /* stuff that needs to be saved/restored */
   Int    groupNo;
   Int    groupPos;
   Int    nblock;
   Int    es;
   Int    N;
   Int    zn; 
   Int    zvec;

   {
      bz_stream* strm = s->strm;
      bsPtr = (UChar *)strm->next_in;
      bsEnd = bsPtr + strm->avail_in;
      bsBuff = s->save_bsBuff;
      bsLive = s->save_bsLive;
   }
   /*restore from the save area*/
   groupNo     = s->save_groupNo;
   groupPos    = s->save_groupPos;
   nblock      = s->save_nblock;

   retVal = BZ_OK;

   switch (s->state) {

/* aliases of saved variables (nblock, groupPos, groupNo are not used
   before decoding of MTF values) */
#define i nblock
#define j groupPos
#define curr groupNo

/* aliases of DState members (tPos, k0, nblock_used are not used
   before inverse BWT */
#define EOB (s->tPos)
#define nGroups (s->k0)
#define nSelectors (s->nblock_used)

      case BZ_X_MAGIC_0:

      GET_UCHAR(BZ_X_MAGIC_1, uc);
      if ((UChar)uc != BZ_HDR_B) RETERR(BZ_DATA_ERROR_MAGIC);
      GET_UCHAR(BZ_X_MAGIC_2, uc);
      if ((UChar)uc != BZ_HDR_Z) RETERR(BZ_DATA_ERROR_MAGIC);
      GET_UCHAR(BZ_X_MAGIC_3, uc);
      if ((UChar)uc != BZ_HDR_h) RETERR(BZ_DATA_ERROR_MAGIC);

      GET_UCHAR(BZ_X_MAGIC_4, s->blockSize100k)
      if (s->blockSize100k < (BZ_HDR_0 + 1) || 
          s->blockSize100k > (BZ_HDR_0 + 9)) RETERR(BZ_DATA_ERROR_MAGIC);
      s->blockSize100k -= BZ_HDR_0;

      if (s->smallDecompress) {
         bz_stream* strm = s->strm;
         s->ll16 = BZALLOC( s->blockSize100k * 100000 * sizeof(UInt16) );
         s->ll4  = BZALLOC( 
                      ((1 + s->blockSize100k * 100000) >> 1) * sizeof(UChar) 
                   );
         if (s->ll16 == NULL || s->ll4 == NULL) RETERR(BZ_MEM_ERROR);
      } else {
         bz_stream* strm = s->strm;
         s->tt  = BZALLOC( s->blockSize100k * 100000 * sizeof(Int32) );
         if (s->tt == NULL) RETERR(BZ_MEM_ERROR);
      }

      case BZ_X_BLKHDR_0:
      GET_BITS_SMALL(BZ_X_BLKHDR_1, uc, 8)

      if (uc == 0x17) goto endhdr_2;
      if (uc != 0x31) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_BLKHDR_2, uc);
      if ((UChar)uc != 0x41) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_BLKHDR_3, uc);
      if ((UChar)uc != 0x59) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_BLKHDR_4, uc);
      if ((UChar)uc != 0x26) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_BLKHDR_5, uc);
      if ((UChar)uc != 0x53) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_BLKHDR_6, uc);
      if ((UChar)uc != 0x59) RETERR(BZ_DATA_ERROR);

      s->currBlockNo++;
      if (s->verbosity >= 2)
         VPrintf1 ( "\n    [%d: huff+mtf ", s->currBlockNo );
 
      s->storedBlockCRC = 0;
      GET_UCHAR(BZ_X_BCRC_1, uc);
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32)uc);
      GET_UCHAR(BZ_X_BCRC_2, uc);
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32)uc);
      GET_UCHAR(BZ_X_BCRC_3, uc);
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32)uc);
      GET_UCHAR(BZ_X_BCRC_4, uc);
      s->storedBlockCRC = (s->storedBlockCRC << 8) | ((UInt32)uc);

      GET_BIT(BZ_X_RANDBIT, s->blockRandomised);

      s->origPtr = 0;
      GET_UCHAR(BZ_X_ORIGPTR_1, uc);
      s->origPtr = (s->origPtr << 8) | ((Int32)uc);
      GET_UCHAR(BZ_X_ORIGPTR_2, uc);
      s->origPtr = (s->origPtr << 8) | ((Int32)uc);
      GET_UCHAR(BZ_X_ORIGPTR_3, uc);
      s->origPtr = (s->origPtr << 8) | ((Int32)uc);

      if (s->origPtr < 0)
         RETERR(BZ_DATA_ERROR);
      if (s->origPtr > 100000*s->blockSize100k) 
         RETERR(BZ_DATA_ERROR);

      /*--- Receive the mapping table ---*/
      for (i = 0; i < 16; i++) {
         GET_BIT(BZ_X_MAPPING_1, uc);
         if (uc == 1) 
            s->mtfa[i + 256] = True; else 
            s->mtfa[i + 256] = False;
      }

      curr = 0;
      for (j = 0; j < 16; j++)
         if (s->mtfa[j + 256])
            for (i = 0; i < 16; i++) {
               GET_BIT(BZ_X_MAPPING_2, uc);
               if (uc == 1) {
                  s->mtfa[curr] = j * 16 + i;
                  curr++;
               }
            }

      if (curr == 0) RETERR(BZ_DATA_ERROR);
      EOB = curr + 1;

      /*--- Now the selectors ---*/
      GET_BITS(BZ_X_SELECTOR_1, uc, 18);
      nGroups = (uc >> 15) - 1;
      nSelectors = uc & ((1 << 15) - 1); 
      if (nGroups < 1 || nGroups > 5) RETERR(BZ_DATA_ERROR);
      if (nSelectors < 1  || nSelectors > BZ_MAX_SELECTORS) RETERR(BZ_DATA_ERROR);
      for (i = nSelectors - 1, j = nGroups, curr = 0; i >= 0; i--) {
         while (True) {
            GET_BIT(BZ_X_SELECTOR_2, uc);
            if (uc == 0) break;
            if (curr >= j) RETERR(BZ_DATA_ERROR);
            curr++;
         }
         s->selector[i] = curr;
         curr = 0;
      }

      /*--- Undo the MTF values for the selectors. ---*/
      {
         UChar pos[BZ_N_GROUPS], tmp;
         UInt v;
         Int  i;
         for (v = 0; (Int32)v <= nGroups; v++) pos[v] = v;
   
         for (i = nSelectors - 1, tmp = 0; i >= 0; i--) {
            v = s->selector[i];
            if (v > 0) {
               tmp = pos[v];
               do { v--; pos[v + 1] = pos[v]; } while (v > 0);
               pos[0] = tmp;
            }
            s->selector[i] = tmp;
         }
      }

      /*--- Now the coding tables ---*/
      for (j = 0; (Int32)j <= nGroups; j++) {
         GET_BITS_SMALL(BZ_X_CODING_1, curr, 5);
         for (i = 0; (UInt32)i <= EOB; i++) {
            while (True) {
               if (curr < 1 || curr > BZ_MAX_DECODE_LEN) RETERR(BZ_DATA_ERROR);
               GET_BITS_SMALL(BZ_X_CODING_2, uc, 2);
               if ((uc & 2) == 0) break;
               if ((uc & 1) == 0) curr++; else curr--;
            }
            UNGET_BITS(1);
            s->groupRec[j].len[i] = curr;
         }
      }

      /*--- Create the Huffman decoding tables ---*/
      {
         Int   i;
         Int32 startLen;

         for (i = 0; (Int32)i <= nGroups; i++) {
            startLen = BZ2_hbCreateDecodeTables ( 
               &(s->groupRec[i].limit[0]), 
               &(s->groupRec[i].base[0]), 
               &(s->groupRec[i].perm[0]), 
               &(s->groupRec[i].len[0]),
               EOB
             );
             if (startLen < 0) RETERR(BZ_DATA_ERROR);
         }
      }

      /*--- Now the MTF values ---*/

      groupNo  = nSelectors;
      groupPos = 0;
      nblock = 100000 * s->blockSize100k;
      if (s->smallDecompress) {
         s->eblock = s->ll16 + nblock;
      } else {
         s->eblock = s->tt + nblock;
      }
      nblock = -nblock;
      gRec = NULL;

#undef i
#undef j
#undef curr

      {
         Int i;

         for (i = 0; i < 256; i++) s->cftab[i] = 0;
   
         /*-- MTF init --*/
         {
            UChar *p;
            for (i = EOB - 1, p = &(s->mtfa[MTFA_SIZE]); i != 0; ) *--p = s->mtfa[--i];
            for (i = 0; p < &(s->mtfa[MTFA_SIZE]); i++) {
               s->mtfbase[i] = p;
               p += MTFL_SIZE;
            }
         }
         /*-- end MTF init --*/
      }

      while (True) {
         UInt uc;

         GET_MTF_VAL(BZ_X_MTF_1, BZ_X_MTF_2, nextSym);

         if (nextSym >= 256) {

            es = nextSym;
            N = 1;
            while (True) {
               GET_MTF_VAL_ES(BZ_X_MTF_3, BZ_X_MTF_4, nextSym);
               if (nextSym < 256) break;
               /* Check that es doesn't get too big, so that it doesn't
                  go negative. The maximum value that can be
                  RUNA/RUNB encoded is equal to the block size (post
                  the initial RLE), viz, 900k, so bounding es at 1M
                  should guard against overflow without
                  rejecting any legitimate inputs. */
               if (es > 1024*1024*256) RETERR(BZ_DATA_ERROR);
               es += nextSym << N;
               N++;
            }

            uc = *(s->mtfbase[0]);
            es >>= 8;
            s->cftab[uc] += es;
            if ((es += nblock) > 0) RETERR(BZ_DATA_ERROR);
             
            if (s->smallDecompress) {
               UInt16 *p = s->eblock;
               do {
                  p[nblock] = (UInt16)uc;
                  nblock++;
               } while (nblock < es);
            } else {
               UInt32 *p = s->eblock;
               do {
                  p[nblock] = (UInt32)uc;
                  nblock++;
               } while (nblock < es);
            }
         }

         {

            if (nblock == 0) {
               if (nextSym == 0) goto EndOfBlock;
               RETERR(BZ_DATA_ERROR);
            }

            /*-- uc = MTF ( nextSym ) --*/
            {
               UChar *pp, *p0;
               UInt  lno, off;

               if (nextSym < MTFL_SIZE) {
                  /* avoid general-case expense */
                  p0 = s->mtfbase[0];
                  uc = p0[nextSym];
                  switch (nextSym) {
                     case 0:  goto EndOfBlock;
                     case 15: p0[15] = p0[14];
                     case 14: p0[14] = p0[13];
                     case 13: p0[13] = p0[12];
                     case 12: p0[12] = p0[11];
                     case 11: p0[11] = p0[10];
                     case 10: p0[10] = p0[9];
                     case 9:  p0[9]  = p0[8];
                     case 8:  p0[8]  = p0[7];
                     case 7:  p0[7]  = p0[6];
                     case 6:  p0[6]  = p0[5];
                     case 5:  p0[5]  = p0[4];
                     case 4:  p0[4]  = p0[3];
                     case 3:  p0[3]  = p0[2];
                     case 2:  p0[2]  = p0[1];
                     case 1:  p0[1]  = p0[0];
                  }
                  *p0 = uc;
               } else { 
                  /* general case */
                  lno = nextSym / MTFL_SIZE;
                  off = nextSym % MTFL_SIZE;
                  p0 = s->mtfbase[lno];
                  pp = s->mtfbase[0];
                  pp[-1] = p0[off];
                  pp = s->mtfbase[lno - 1];
                  switch (off) {
                     case 15: p0[15] = p0[14];
                     case 14: p0[14] = p0[13];
                     case 13: p0[13] = p0[12];
                     case 12: p0[12] = p0[11];
                     case 11: p0[11] = p0[10];
                     case 10: p0[10] = p0[9];
                     case 9:  p0[9]  = p0[8];
                     case 8:  p0[8]  = p0[7];
                     case 7:  p0[7]  = p0[6];
                     case 6:  p0[6]  = p0[5];
                     case 5:  p0[5]  = p0[4];
                     case 4:  p0[4]  = p0[3];
                     case 3:  p0[3]  = p0[2];
                     case 2:  p0[2]  = p0[1];
                     case 1:  p0[1]  = p0[0];
                  }
                  switch (lno) {
                     case 15: *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[13];
                              s->mtfbase[14] = p0;
                     case 14: *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[12];
                              s->mtfbase[13] = p0;
                     case 13: *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[11];
                              s->mtfbase[12] = p0;
                     case 12: *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[10];
                              s->mtfbase[11] = p0;
                     case 11: *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[9];
                              s->mtfbase[10] = p0;
                     case 10: *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[8];
                              s->mtfbase[9] = p0;
                     case 9:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[7];
                              s->mtfbase[8] = p0;
                     case 8:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[6];
                              s->mtfbase[7] = p0;
                     case 7:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[5];
                              s->mtfbase[6] = p0;
                     case 6:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[4];
                              s->mtfbase[5] = p0;
                     case 5:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[3];
                              s->mtfbase[4] = p0;
                     case 4:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[2];
                              s->mtfbase[3] = p0;
                     case 3:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[1];
                              s->mtfbase[2] = p0;
                     case 2:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              pp = s->mtfbase[0];
                              s->mtfbase[1] = p0;
                     case 0:  /*--not reached--*/
                     case 1:  *p0 = pp[MTFL_SIZE - 1];
                              p0 = pp - 1;
                              s->mtfbase[0] = p0;
                  }
                  uc = *p0;
                  if (p0 == s->mtfa) {
                     lno = (EOB - 2) / MTFL_SIZE;
                     p0 = s->mtfbase[lno];
                     while (lno != 0) {
                        lno--;
                        pp = s->mtfbase[lno] + MTFL_SIZE;
                        do {
                           *--p0 = *--pp;
                        } while (pp > s->mtfbase[lno]);
                        s->mtfbase[lno] = p0;
                     }
                  }
               }
            }
            /*-- end uc = MTF ( nextSym ) --*/

            s->cftab[uc]++;
            if (s->smallDecompress)
               ((UInt16*)(s->eblock))[nblock] = (UInt16)(uc); else
               ((UInt32*)(s->eblock))[nblock] = (UInt32)(uc);
            nblock++;
         }
      }

#undef EOB
#undef nGroups
#undef nSelectors

      EndOfBlock:
      nblock += 100000 * s->blockSize100k;
      /* Now we know what nblock is, we can do a better sanity
         check on s->origPtr.
      */
      if (s->origPtr < 0 || s->origPtr >= nblock)
         RETERR(BZ_DATA_ERROR);

      s->save_bsBuff = bsBuff;
      s->save_bsLive = bsLive;
      {
         bz_stream* strm = s->strm;
         UInt32     processed = bsPtr - (UChar *)strm->next_in;
         UInt32     new_total_in_lo32;
   
         strm->next_in = (char *)bsPtr;
         new_total_in_lo32 = strm->total_in_lo32 + processed;
         if(new_total_in_lo32 < strm->total_in_lo32) strm->total_in_hi32++;
         strm->total_in_lo32 = new_total_in_lo32;
         strm->avail_in -= processed;
      }

      /* let block = (s->smallDecompress ? s->ll16 : s->tt),
         let c[v] = cftab[v]
         now block[i] >= 0 and block[i] < 256 for each i >= 0 and i < nblock,
         cftab[v] is count of value v in block
         thus cftab[v] >= 0 for each v >= 0 and v < 256
         and sum(cftab[v]) for v from 0 to 255 is equal to nblock
      */
      s->nblock_used = nblock;
      s->state_out_len = 0;
      s->state_out_ch  = 0;
      BZ_INITIALISE_CRC ( s->calculatedBlockCRC );
      s->state = BZ_X_OUTPUT;
      if (s->verbosity >= 2) VPrintf0 ( "rt+rld" );

      {
         Int32 i, j;

         /* Actually generate cftab. */
         for (i = 1; i < 256; i++) s->cftab[i] += s->cftab[i-1];
         /* now cftab[0] == c[0], cftab[1] == c[0] + c[1], ...,
            cftab[255] == c[0] + c[1] + ... + c[255]
         */ 

         if (s->smallDecompress) {
            UInt16* c_ll16 = s->ll16;
            UChar*  c_ll4 = s->ll4;
   
            /*-- compute the T vector --*/
            i = nblock - 1;
            if ((i & 1) == 0) {
               uc = (UChar)(c_ll16[i]);
               s->cftab[uc]--;
               c_ll16[i] = (UInt16)(s->cftab[uc] & 0x0000ffff);
               c_ll4[(i) >> 1] = (s->cftab[uc] >> 16);
               i--;
            }
            while (i >= 0) {
               uc = (UChar)(c_ll16[i]);
               s->cftab[uc]--;
               c_ll16[i] = (UInt16)(s->cftab[uc] & 0x0000ffff);
               j = (s->cftab[uc] >> 16) << 4;
               i--;
               uc = (UChar)(c_ll16[i]);
               s->cftab[uc]--;
               c_ll16[i] = (UInt16)(s->cftab[uc] & 0x0000ffff);
               c_ll4[(i) >> 1] = (s->cftab[uc] >> 16) | j;
               i--;
            }
            /* now cftab[0] == 0, cftab[1] == c[0],
               cftab[2] = c[0] + c[1], ...,
               cftab[255] == c[0] + c[1] + ... + c[254];
               
            */ 
   
            /*-- Compute T^(-1) by pointer reversal on T --*/
            i = s->origPtr;
            j = GET_LL_C(i);
            do {
               Int32 tmp = GET_LL_C(j);
               SET_LL_C(j, i);
               i = j;
               j = tmp;
            } while (i != s->origPtr);

            BZ2_indexCreateTab(&(s->cftab[0]), &(s->selector[0]), nblock);

            s->tPos = s->origPtr;
            if (s->blockRandomised) {
               BZ_RAND_INIT_MASK;
               BZ_GET_SMALL(s->k0);
               BZ_RAND_UPD_MASK; s->k0 ^= BZ_RAND_MASK; 
            } else {
               BZ_GET_SMALL(s->k0);
            }
   
         } else {
            UInt32* c_tt = s->tt;
   
            /*-- compute the T^(-1) vector --*/
            for (i = nblock - 1; i >= 0; i--) {
               uc = (UChar)(c_tt[i] & 0xff);
               s->cftab[uc]--;
               c_tt[s->cftab[uc]] |= (i << 8);
            }
            /* now (s->tt[i] >> 8) >= 0 and (s->tt[i] >> 8) < nblock
               for each i >= 0 and i < nblock
               and each value (s->tt[i] >> 8) is unique and assigned once
            */
   
            s->tPos = c_tt[c_tt[s->origPtr] >> 8];
            if (s->blockRandomised) {
               BZ_RAND_INIT_MASK;
               BZ_GET_FAST(s->k0);
               BZ_RAND_UPD_MASK; s->k0 ^= BZ_RAND_MASK; 
            } else {
               BZ_GET_FAST(s->k0);
            }
   
         }
      }

      return(BZ_OK);


    endhdr_2:

      GET_UCHAR(BZ_X_ENDHDR_2, uc);
      if ((UChar)uc != 0x72) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_ENDHDR_3, uc);
      if ((UChar)uc != 0x45) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_ENDHDR_4, uc);
      if ((UChar)uc != 0x38) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_ENDHDR_5, uc);
      if ((UChar)uc != 0x50) RETERR(BZ_DATA_ERROR);
      GET_UCHAR(BZ_X_ENDHDR_6, uc);
      if ((UChar)uc != 0x90) RETERR(BZ_DATA_ERROR);

      s->storedCombinedCRC = 0;
      GET_UCHAR(BZ_X_CCRC_1, uc);
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32)uc);
      GET_UCHAR(BZ_X_CCRC_2, uc);
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32)uc);
      GET_UCHAR(BZ_X_CCRC_3, uc);
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32)uc);
      GET_UCHAR(BZ_X_CCRC_4, uc);
      s->storedCombinedCRC = (s->storedCombinedCRC << 8) | ((UInt32)uc);

//      RETERR(BZ_STREAM_END);
      retVal = BZ_STREAM_END;
      error_return:
      s->state = BZ_X_IDLE;
   
      s->save_bsBuff = bsBuff;
      s->save_bsLive = bsLive;
      {
         bz_stream* strm = s->strm;
         UInt32     processed = bsPtr - (UChar *)strm->next_in;
         UInt32     new_total_in_lo32;
   
         strm->next_in = (char *)bsPtr;
         new_total_in_lo32 = strm->total_in_lo32 + processed;
         if(new_total_in_lo32 < strm->total_in_lo32) strm->total_in_hi32++;
         strm->total_in_lo32 = new_total_in_lo32;
         strm->avail_in -= processed;
      }
   
      return retVal;   

      default: AssertH ( False, 4001 );
   }

   AssertH ( False, 4002 );

   save_state_and_return:

   s->save_groupNo     = groupNo;
   s->save_groupPos    = groupPos;
   s->save_nblock      = nblock;

   s->save_bsBuff = bsBuff;
   s->save_bsLive = bsLive;
   {
      bz_stream* strm = s->strm;
      UInt32     processed = bsPtr - (UChar *)strm->next_in;
      UInt32     new_total_in_lo32;

      strm->next_in = (char *)bsPtr;
      new_total_in_lo32 = strm->total_in_lo32 + processed;
      if(new_total_in_lo32 < strm->total_in_lo32) strm->total_in_hi32++;
      strm->total_in_lo32 = new_total_in_lo32;
      strm->avail_in -= processed;
   }

   return BZ_OK;
}


/*-------------------------------------------------------------*/
/*--- end                                      decompress.c ---*/
/*-------------------------------------------------------------*/
