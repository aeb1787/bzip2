
/*-------------------------------------------------------------*/
/*--- Compression machinery (not incl block sorting)        ---*/
/*---                                            compress.c ---*/
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


/* CHANGES
    0.9.0    -- original version.
    0.9.0a/b -- no changes in this file.
    0.9.0c   -- changed setting of nGroups in sendMTFValues() 
                so as to do a bit better on small files
*/

#include "bzlib_private.h"


/*---------------------------------------------------*/
/*--- Bit stream I/O                              ---*/
/*---------------------------------------------------*/

/*---------------------------------------------------*/
void BZ2_bsInitWrite ( EState* s )
{
   s->bsLive = sizeof(s->bsBuff) * 8;
   s->bsBuff = 0;
}


/*---------------------------------------------------*/
static
void bsFinishWrite ( EState* s )
{
   while (s->bsLive < sizeof(s->bsBuff) * 8) {
      s->zbits[s->numZ] = (UChar)(s->bsBuff >> (sizeof(s->bsBuff) * 8 - 8));
      s->numZ++;
      s->bsBuff <<= 8;
      s->bsLive += 8;
   }
}


/*---------------------------------------------------*/
#define bsNEEDW(nz)                                           \
{                                                             \
   while (s->bsLive <= (sizeof(s->bsBuff) * 8 - 8)) {         \
      s->zbits[s->numZ]                                       \
         = (UChar)(s->bsBuff >> (sizeof(s->bsBuff) * 8 - 8)); \
      s->numZ++;                                              \
      s->bsBuff <<= 8;                                        \
      s->bsLive += 8;                                         \
   }                                                          \
}


/*---------------------------------------------------*/
static
__inline__
void bsW ( EState* s, Int32 n, UInt v )
{
   bsNEEDW ( n );
   s->bsBuff |= (v << (s->bsLive -= n));
}


/*---------------------------------------------------*/
static
void bsPutUInt32 ( EState* s, UInt u )
{
   bsW ( s, 16, (u >> 16) );
   bsW ( s, 16, (u & 0xffff) );
}

/*---------------------------------------------------*/
/*--- The back end proper                         ---*/
/*---------------------------------------------------*/

/*---------------------------------------------------*/
static
void generateMTFValues ( EState* s )
{
   UChar   yy[256];
   UChar*  mtfbase[255 / MTFL_SIZE + 1];
   Int     zPend;
   Int     i;
   Int     EOB;

   /* 
      After sorting (eg, here),
         s->arr1 [ 0 .. s->nblock-1 ] holds BWT values,

      The first thing to do is generate the MTF values,
      and put them in
         ((UInt16*)s->arr1) [ 0 .. s->nblock-1 ].
      Because there are strictly fewer or equal MTF values
      than block values, ptr values in this area are overwritten
      with MTF values only when they are no longer needed.

      The final compressed bitstream is generated into the
      area starting at
         (UChar*) (&((UInt16*)s->arr1)[s->blockSize100k*100000])

      These storage aliases are set up in bzCompressInit().
   */
   UInt32 *ptr;
   UInt16 *mtfv;
   UInt    yy0, ll_i, tmp, lno;
   UChar  *pc, *pe;

   for (pc = (UChar*)(s->ftab) + MTFA_SIZE, i = 255; i >= 0; i--) {
      if (s->inUse[i]) {
         *--pc = i;
      }
   }
   s->nInUse = (UChar*)(s->ftab) + MTFA_SIZE - pc;

   for (i = 0; pc < (UChar*)(s->ftab) + MTFA_SIZE; i++) {
      mtfbase[i] = pc;
      pe = pc + MTFL_SIZE;
      if (pe > (UChar*)(s->ftab) + MTFA_SIZE)
         pe = (UChar*)(s->ftab) + MTFA_SIZE;
      while (pc < pe) {
         yy[*pc] = i;
         pc++;
      }
   }

   for (i = s->nInUse; i >= 0; i--) {
      s->mtfFreq[i] = 0;
   }

   ptr = s->ptr; mtfv = s->mtfv; yy0 = mtfbase[0][0];
   ptr[s->nblock] = 256;
   while (True) {

      AssertD ( mtfv - s->mtfv <= ptr - s->ptr, "generateMTFValues(1)" );

      ll_i = *ptr; ptr++;

      if (yy0 == ll_i) { 
         zPend = 0;
         while (ll_i = *ptr, ptr++, yy0 == ll_i) {
            zPend++;
         }
         while (True) {
            i = (zPend & 1) + BZ_RUNA;

            *mtfv = i; mtfv++; 
            s->mtfFreq[i]++; 

            if ((zPend -= 2) < 0) break;
            zPend >>= 1;
         }
      }
      if (ll_i > 255) break;
      AssertD ( s->inUse[ll_i], "generateMTFValues(2a)" );

      lno = yy[ll_i];
      yy[yy0] = 0;
      if (lno == 0) {
         pc = mtfbase[0];
         tmp   = yy0;
         yy0   = pc[1];
         pc[1] = tmp;
         i = 2;
         while ( ll_i != yy0 ) {
            tmp    = yy0;
            yy0    = pc[i];
            pc[i]  = tmp;
            i++;
         }
         *mtfv = i; mtfv++;
         s->mtfFreq[i]++;
      } else {
         mtfbase[0][0] = yy0;
         pc = mtfbase[lno];
         yy0 = pc[0];
         i = 1;
         while (ll_i != yy0) {
            tmp   = yy0;
            yy0   = pc[i];
            pc[i] = tmp;
            i++;
         }
         i += lno * MTFL_SIZE;
         *mtfv = i; mtfv++;
         s->mtfFreq[i]++;
         do {
            tmp = mtfbase[lno - 1][MTFL_SIZE - 1];
            yy[tmp] = lno;
            mtfbase[lno][0] = tmp;
            mtfbase[lno - 1]--;
            lno--;
         } while (lno != 0);
         if (mtfbase[0] == (UChar*)(s->ftab)) {
            lno = (s->nInUse - 1) / MTFL_SIZE;
            pe = mtfbase[lno];
            while (lno != 0) {
               lno--;
               pc = mtfbase[lno] + MTFL_SIZE;
               do {
                  *--pe = *--pc;
               } while (pc > mtfbase[lno]);
               mtfbase[lno] = pe;
            }
         }
      }
   }
   AssertD ( ptr == s->ptr + s->nblock + 1, "generateMTFValues(2b)" );

   EOB = s->nInUse + 1;
   *mtfv = EOB; mtfv++; s->mtfFreq[EOB] = 1;

   s->nMTF = mtfv - s->mtfv;
}


/*---------------------------------------------------*/
#define BZ_LESSER_ICOST  0
#define BZ_GREATER_ICOST 15

typedef
   struct {
#ifdef USE_INT64
      UInt64   len_pack[BZ_MAX_ALPHA_SIZE];
#else
      UInt32   len_pack[BZ_MAX_ALPHA_SIZE * 2];
#endif
      Int32    code[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
      Int32    len[BZ_N_GROUPS][BZ_MAX_ALPHA_SIZE];
      Int32    cost[BZ_N_GROUPS];
      Int32    fave[BZ_N_GROUPS];
   }
   ETabs;

static
void sendMTFValues ( EState* s )
{
   Int   t;
   Int   nSelectors, alphaSize;
   Int   nGroups;
   ETabs *ftab = (ETabs *)(s->ftab);
   Int32 iter, nBytes;

   if (s->verbosity >= 3)
      VPrintf3( "      %d in block, %d after MTF & 1-2 coding, "
                "%d+2 syms in use\n", 
                s->nblock, s->nMTF, s->nInUse );

   /*--- Decide how many coding tables to use ---*/
   AssertH ( s->nMTF > 0, 3001 );
   if (s->nMTF < 200)  nGroups = 2; else
   if (s->nMTF < 600)  nGroups = 3; else
   if (s->nMTF < 1200) nGroups = 4; else
   if (s->nMTF < 2400) nGroups = 5; else
                       nGroups = 6;

   alphaSize = s->nInUse+2;
   if (alphaSize <= 0) return;

   /*--- Generate an initial set of coding tables ---*/
   {
      Int   t, v;
      Int32 nPart, remF, tFreq, aFreq;
      Int32 gs, ge;
   
      for (t = 0; t < nGroups; t++)
         for (v = 0; v < alphaSize; v++)
            ftab->len[t][v] = BZ_GREATER_ICOST;
   
      nPart = nGroups;
      remF  = s->nMTF;
      gs = 0;
      while (nPart > 0) {
         tFreq = remF / nPart;
         ge = gs - 1;
         aFreq = 0;
         while (aFreq < tFreq && ge < alphaSize - 1) {
            ge++;
            aFreq += s->mtfFreq[ge];
         }
   
         if (ge > gs 
             && nPart != nGroups && nPart != 1 
             && ((nGroups - nPart) % 2 == 1)) {
            aFreq -= s->mtfFreq[ge];
            ge--;
         }
   
         if (s->verbosity >= 3)
            VPrintf5( "      initial group %d, [%d .. %d], "
                      "has %d syms (%4.1f%%)\n",
                      nPart, gs, ge, aFreq, 
                      (100.0 * (float)aFreq) / (float)(s->nMTF) );
   
         while (gs <= ge) {
            ftab->len[nPart - 1][gs] = BZ_LESSER_ICOST;
            gs++;
         }
   
         nPart--;
         gs = ge + 1;
         remF -= aFreq;
      }
   }

   /*--- 
      Iterate up to BZ_N_ITERS times to improve the tables.
   ---*/
   for (iter = 0; iter < BZ_N_ITERS; iter++) {
      UChar  *sel = (UChar *)(s->arr2);
      UInt16 *emtfv = s->mtfv + s->nMTF;
      UInt16 *gs = s->mtfv;
      UInt16 *ge = s->mtfv + BZ_G_SIZE;
      Int    t, v;
      Int32  totc = 0;
   
      for (t = 0; t < nGroups; t++) ftab->fave[t] = 0;
   
      for (t = 0; t < nGroups; t++)
         for (v = 0; v < alphaSize; v++)
            ftab->code[t][v] = 0;
   
      if (nGroups == 6) {
         /*---
           Set up an auxiliary length table which is used to fast-track
           the common case (nGroups == 6). 
         ---*/
         for (v = 0; v < alphaSize; v++) {
#ifdef USE_INT64
            ftab->len_pack[v] = ((Int64)(ftab->len[5][v]) << 50) | ((Int64)(ftab->len[4][v]) << 40) |
                                ((Int64)(ftab->len[3][v]) << 30) | ((Int64)(ftab->len[2][v]) << 20) |
                                ((Int64)(ftab->len[1][v]) << 10) | (Int64)(ftab->len[0][v]);
#else
            ftab->len_pack[v * 2] = (ftab->len[2][v] << 20) | (ftab->len[1][v] << 10) | ftab->len[0][v];
            ftab->len_pack[v * 2 + 1] = (ftab->len[5][v] << 20) | (ftab->len[4][v] << 10) | ftab->len[3][v];
#endif
         }
         while (ge <= emtfv) {
            Int    bt;
            Int32  bc;
#ifdef USE_INT64
            Int64 cost05 = 0;
   
#           define BZ_ITER(nn) (cost05 += ftab->len_pack[gs[nn]])
#else
            UInt16 *gi = gs;
            Int32 cost02 = 0, cost35 = 0;
   
#           define BZ_ITER(nn)                                     \
               cost02 += ((Int32 *)(ftab->len_pack))[*gi * 2];     \
               cost35 += ((Int32 *)(ftab->len_pack))[*gi * 2 + 1]; \
               gi++;                                               \

#endif

#if BZ_G_SIZE != 50
#error BZ_G_SIZE != 50: Next code is incorrect
#endif

            BZ_ITER(0);  BZ_ITER(1);  BZ_ITER(2);  BZ_ITER(3);  BZ_ITER(4);
            BZ_ITER(5);  BZ_ITER(6);  BZ_ITER(7);  BZ_ITER(8);  BZ_ITER(9);
            BZ_ITER(10); BZ_ITER(11); BZ_ITER(12); BZ_ITER(13); BZ_ITER(14);
            BZ_ITER(15); BZ_ITER(16); BZ_ITER(17); BZ_ITER(18); BZ_ITER(19);
            BZ_ITER(20); BZ_ITER(21); BZ_ITER(22); BZ_ITER(23); BZ_ITER(24);
            BZ_ITER(25); BZ_ITER(26); BZ_ITER(27); BZ_ITER(28); BZ_ITER(29);
            BZ_ITER(30); BZ_ITER(31); BZ_ITER(32); BZ_ITER(33); BZ_ITER(34);
            BZ_ITER(35); BZ_ITER(36); BZ_ITER(37); BZ_ITER(38); BZ_ITER(39);
            BZ_ITER(40); BZ_ITER(41); BZ_ITER(42); BZ_ITER(43); BZ_ITER(44);
            BZ_ITER(45); BZ_ITER(46); BZ_ITER(47); BZ_ITER(48); BZ_ITER(49);
   
#           undef BZ_ITER
   
            /*-- 
               Find the coding table which is best for this group,
               and record its identity in the selector table.
            --*/
#ifdef USE_INT64
            bc = cost05 & 0x3ff; bt = 0;
            for (t = 1; t < nGroups; t++) {
               cost05 >>= 10;
               if ((Int32)(cost05 & 0x3ff) < bc) { bc = (Int32)(cost05 & 0x3ff); bt = t; };
            }
#else
            bc = cost02 & 0x3ff; ftab->cost[1] = (cost02 >> 10) & 0x3ff; ftab->cost[2] = cost02 >> 20;
            ftab->cost[3] = cost35 & 0x3ff; ftab->cost[4] = (cost35 >> 10) & 0x3ff; ftab->cost[5] = cost35 >> 20;
            bt = 0;
            for (t = 1; t < nGroups; t++)
               if (ftab->cost[t] < bc) { bc = ftab->cost[t]; bt = t; };
#endif
    
            totc += bc;
            ftab->fave[bt]++;
            *sel++ = bt;
   
            /*-- 
               Increment the symbol frequencies for the selected table.
             --*/
            {
               /*--- fast track the common case ---*/
               Int32 *ifr = &(ftab->code[bt][0]);
   
#              define BZ_ITUR(nn) ifr[gs[nn]]++

#if BZ_G_SIZE != 50
#error BZ_G_SIZE != 50: Next code is incorrect
#endif
   
               BZ_ITUR(0);  BZ_ITUR(1);  BZ_ITUR(2);  BZ_ITUR(3);  BZ_ITUR(4);
               BZ_ITUR(5);  BZ_ITUR(6);  BZ_ITUR(7);  BZ_ITUR(8);  BZ_ITUR(9);
               BZ_ITUR(10); BZ_ITUR(11); BZ_ITUR(12); BZ_ITUR(13); BZ_ITUR(14);
               BZ_ITUR(15); BZ_ITUR(16); BZ_ITUR(17); BZ_ITUR(18); BZ_ITUR(19);
               BZ_ITUR(20); BZ_ITUR(21); BZ_ITUR(22); BZ_ITUR(23); BZ_ITUR(24);
               BZ_ITUR(25); BZ_ITUR(26); BZ_ITUR(27); BZ_ITUR(28); BZ_ITUR(29);
               BZ_ITUR(30); BZ_ITUR(31); BZ_ITUR(32); BZ_ITUR(33); BZ_ITUR(34);
               BZ_ITUR(35); BZ_ITUR(36); BZ_ITUR(37); BZ_ITUR(38); BZ_ITUR(39);
               BZ_ITUR(40); BZ_ITUR(41); BZ_ITUR(42); BZ_ITUR(43); BZ_ITUR(44);
               BZ_ITUR(45); BZ_ITUR(46); BZ_ITUR(47); BZ_ITUR(48); BZ_ITUR(49);
   
#              undef BZ_ITUR
   
            }
   
            gs = ge;
            ge += BZ_G_SIZE;
         }
      }
   
      /*--
         Groups of size BZ_G_SIZE when nGroups != 6
      --*/
      while (ge <= emtfv) {
         UInt16 *gi;
         Int    bt;
         Int32  bc;
   
         /*-- 
            Calculate the cost of this group as coded
            by each of the coding tables.
         --*/
         for (t = 0; t < nGroups; t++) ftab->cost[t] = 0;
         for (gi = gs; gi < ge; gi++) { 
            Int icv = *gi;
            for (t = 0; t < nGroups; t++) ftab->cost[t] += ftab->len[t][icv];
         }
    
         /*-- 
            Find the coding table which is best for this group,
            and record its identity in the selector table.
         --*/
         bc = ftab->cost[0]; bt = 0;
         for (t = 1; t < nGroups; t++)
            if (ftab->cost[t] < bc) { bc = ftab->cost[t]; bt = t; };
         totc += bc;
         ftab->fave[bt]++;
         *sel++ = bt;
   
         /*-- 
            Increment the symbol frequencies for the selected table.
          --*/
         {
            /*--- fast track the common case ---*/
            Int32 *ifr = &(ftab->code[bt][0]);
   
#           define BZ_ITUR(nn) ifr[gs[nn]]++

#if BZ_G_SIZE != 50
#error BZ_G_SIZE != 50: Next code is incorrect
#endif
   
            BZ_ITUR(0);  BZ_ITUR(1);  BZ_ITUR(2);  BZ_ITUR(3);  BZ_ITUR(4);
            BZ_ITUR(5);  BZ_ITUR(6);  BZ_ITUR(7);  BZ_ITUR(8);  BZ_ITUR(9);
            BZ_ITUR(10); BZ_ITUR(11); BZ_ITUR(12); BZ_ITUR(13); BZ_ITUR(14);
            BZ_ITUR(15); BZ_ITUR(16); BZ_ITUR(17); BZ_ITUR(18); BZ_ITUR(19);
            BZ_ITUR(20); BZ_ITUR(21); BZ_ITUR(22); BZ_ITUR(23); BZ_ITUR(24);
            BZ_ITUR(25); BZ_ITUR(26); BZ_ITUR(27); BZ_ITUR(28); BZ_ITUR(29);
            BZ_ITUR(30); BZ_ITUR(31); BZ_ITUR(32); BZ_ITUR(33); BZ_ITUR(34);
            BZ_ITUR(35); BZ_ITUR(36); BZ_ITUR(37); BZ_ITUR(38); BZ_ITUR(39);
            BZ_ITUR(40); BZ_ITUR(41); BZ_ITUR(42); BZ_ITUR(43); BZ_ITUR(44);
            BZ_ITUR(45); BZ_ITUR(46); BZ_ITUR(47); BZ_ITUR(48); BZ_ITUR(49);
   
#           undef BZ_ITUR
            }
      
         gs = ge;
         ge += BZ_G_SIZE;
      }
   
      /*--
         If last group is not full.
      --*/
      if (gs < emtfv) {
         UInt16 *gi;
         Int    bt;
         Int32  bc;
   
         /*-- 
            Calculate the cost of this group as coded
            by each of the coding tables.
         --*/
         for (t = 0; t < nGroups; t++) ftab->cost[t] = 0;
         for (gi = gs; gi < emtfv; gi++) { 
            Int icv = *gi;
            for (t = 0; t < nGroups; t++) ftab->cost[t] += ftab->len[t][icv];
         }
    
         /*-- 
            Find the coding table which is best for this group,
            and record its identity in the selector table.
         --*/
         bc = ftab->cost[0]; bt = 0;
         for (t = 1; t < nGroups; t++)
            if (ftab->cost[t] < bc) { bc = ftab->cost[t]; bt = t; };
         totc += bc;
         ftab->fave[bt]++;
         *sel++ = bt;
   
         /*-- 
            Increment the symbol frequencies for the selected table.
          --*/
         {
            /*--- slow version which correctly handles all situations ---*/
            Int32 *ifr = &(ftab->code[bt][0]);
   
            for (; gs < emtfv; gs++) ifr[*gs]++;
         }
      }
   
   
      if (s->verbosity >= 3) {
         VPrintf2 ( "      pass %d: size is %d, grp uses are ", 
                   iter + 1, totc / 8 );
         for (t = 0; t < nGroups; t++)
            VPrintf1 ( "%d ", ftab->fave[t] );
         VPrintf0 ( "\n" );
      }
   
      nSelectors = sel - (UChar *)(s->arr2);
   
      /*--
        Recompute the tables based on the accumulated frequencies.
      --*/
      /* maxLen was changed from 20 to 17 in bzip2-1.0.3.  See 
         comment in huffman.c for details. */
      for (t = 0; t < nGroups; t++)
         BZ2_hbMakeCodeLengths ( &(ftab->len[t][0]), &(ftab->code[t][0]), alphaSize, BZ_LIMIT_CODE_LEN /*20*/ );
   }

   AssertH( nGroups < 8, 3002 );
   AssertH( nSelectors < 32768 &&
            nSelectors <= BZ_MAX_SELECTORS,
            3003 );

   /*--- Compute MTF values for the selectors. ---*/
#if BZ_N_GROUPS > 6
#error BZ_N_GROUPS > 6: Next code is incorrect
#endif

   {
      UChar  *sel = (UChar *)(s->arr2);
      UInt   pos0 = 0;
      UInt   pos1 = 1;
      UInt   pos2 = 2;
      UInt   pos3 = 3;
      UInt   pos4 = 4;
      UInt   pos5 = 5;
   
      sel[nSelectors] = 8;
      while(1) {
         UInt ll_i = *sel++;
         if (ll_i == pos0) {
            sel[BZ_MAX_SELECTORS - 1] = 1;
            continue;
         }
         if (ll_i == pos1) {
            sel[BZ_MAX_SELECTORS - 1] = 2;
            pos1 = pos0;
            pos0 = ll_i;
            continue;
         }
         if (ll_i == pos2) {
            sel[BZ_MAX_SELECTORS - 1] = 3;
            pos2 = pos1;
            pos1 = pos0;
            pos0 = ll_i;
            continue;
         }
         if (ll_i == pos3) {
            sel[BZ_MAX_SELECTORS - 1] = 4;
            pos3 = pos2;
            pos2 = pos1;
            pos1 = pos0;
            pos0 = ll_i;
            continue;
         }
         if (ll_i == pos4) {
            sel[BZ_MAX_SELECTORS - 1] = 5;
            pos4 = pos3;
            pos3 = pos2;
            pos2 = pos1;
            pos1 = pos0;
            pos0 = ll_i;
            continue;
         }
         if (ll_i == pos5) {
            sel[BZ_MAX_SELECTORS - 1] = 6;
            pos5 = pos4;
            pos4 = pos3;
            pos3 = pos2;
            pos2 = pos1;
            pos1 = pos0;
            pos0 = ll_i;
            continue;
         }
         break;
      }
   }

    /*--- Assign actual codes for the tables. --*/
   {
      Int32  *clen_row = &(ftab->len[0][0]);
      Int    t, v;
   
      for (t = 0; t < nGroups; t++) {
         Int32 minLen = 32;
         Int32 maxLen = 0;
   
         for (v = 0; v < alphaSize; v++) {
            if (clen_row[v] > maxLen) maxLen = clen_row[v];
            if (clen_row[v] < minLen) minLen = clen_row[v];
         }
         AssertH ( !(maxLen > BZ_LIMIT_CODE_LEN /*20*/ ), 3004 );
         AssertH ( !(minLen < 1),  3005 );
         BZ2_hbAssignCodes ( clen_row - sizeof(ftab->code) / sizeof(Int32), clen_row, minLen, maxLen, alphaSize );
         clen_row += BZ_MAX_ALPHA_SIZE;
      }
   }

   /*--- Transmit the mapping table. ---*/
   nBytes = s->numZ;
   { 
      Int   i, j;
      Int   v, t;

      for (i = 0, t = 0; i < 16; i++) {
          t <<= 1;
          for (j = 0; j < 16; j++)
             t |= s->inUse[i * 16 + j];
      }
      bsW(s, 16, t);

      for (i = 0; i < 16; i++) {
         if (t & 0x8000) {
            for (j = 0, v = 0; j < 16; j++) {
               v = (v << 1) | s->inUse[i * 16 + j];
            }
            bsW(s, 16, v);
         }
         t <<= 1;
      }
   }
   if (s->verbosity >= 3) 
      VPrintf1( "      bytes: mapping %d, ", s->numZ - nBytes );

   /*--- Now the selectors. ---*/
   nBytes = s->numZ;
   bsW ( s, 3, nGroups );
   bsW ( s, 15, nSelectors );
   {
      UChar *sel, *esel;

      for(sel = (UChar *)(s->arr2) + BZ_MAX_SELECTORS, esel = sel + nSelectors; sel < esel; sel++) {
         bsW(s, *sel, (1 << *sel) - 2);
      }
   }
   if (s->verbosity >= 3)
      VPrintf1( "selectors %d, ", s->numZ - nBytes );

   /*--- Now the coding tables. ---*/
   nBytes = s->numZ;

   for (t = 0; t < nGroups; t++) {
      Int   v;
      Int32 curr = ftab->len[t][0];

      bsW ( s, 5, curr );
      for (v = 0; v < alphaSize; v++) {
         while (curr < ftab->len[t][v]) { bsW(s,2,2); curr++; /* 10 */ };
         while (curr > ftab->len[t][v]) { bsW(s,2,3); curr--; /* 11 */ };
         bsW ( s, 1, 0 );
      }
   }

   if (s->verbosity >= 3)
      VPrintf1 ( "code lengths %d, ", s->numZ - nBytes );

   /*--- And finally, the block data proper ---*/

   nBytes = s->numZ;
   {
      UInt16 *ge = s->mtfv + BZ_G_SIZE;
      UInt16 *emtfv = s->mtfv + s->nMTF;
      UChar  *sel = (UChar *)(s->arr2);

#if BZ_LIMIT_CODE_LEN == 17
      UChar  *czbits = &(s->zbits[s->numZ]);
      UInt   bsBuff = s->bsBuff;
      Int32  bsLive = s->bsLive;
#ifdef USE_INT64
#define bsWc(n,v)                             \
{                                             \
   if (bsLive <= 16) {                        \
      czbits[0] = (UChar)(bsBuff >> 56);      \
      czbits[1] = (UChar)(bsBuff >> 48);      \
      czbits[2] = (UChar)(bsBuff >> 40);      \
      czbits[3] = (UChar)(bsBuff >> 32);      \
      czbits[4] = (UChar)(bsBuff >> 24);      \
      czbits[5] = (UChar)(bsBuff >> 16);      \
      czbits += 6;                            \
      bsBuff <<= 48;                          \
      bsLive += 48;                           \
   }                                          \
   bsBuff |= (UInt)(v) << (bsLive -= (n));    \
}
#else
#define bsWc(n,v)                             \
{                                             \
   if (bsLive <= 16) {                        \
      czbits[0] = (UChar)(bsBuff >> 24);      \
      czbits[1] = (UChar)(bsBuff >> 16);      \
      if (bsLive > 8) {                       \
         czbits += 2;                         \
         bsBuff <<= 16;                       \
         bsLive += 16;                        \
      } else {                                \
         czbits[2] = (UChar)(bsBuff >> 8);    \
         czbits += 3;                         \
         bsBuff <<= 24;                       \
         bsLive += 24;                        \
      }                                       \
   }                                          \
   bsBuff |= (UInt)(v) << (bsLive -= (n));    \
}
#endif
#else
#define bsWc(n,v) bsW(s,(n),(v))
#endif
   
      //--- fast version for full group ---
      while (ge <= emtfv) {
         Int32* code_row;
         AssertH ( *sel < nGroups, 3006 );
         code_row = &(ftab->code[*sel][0]);
   
#        define BZ_ITAH(nn) bsWc(code_row[ge[nn - BZ_G_SIZE] + sizeof(ftab->code) / sizeof(Int32)], code_row[ge[nn - BZ_G_SIZE]])

#if BZ_G_SIZE != 50
#error BZ_G_SIZE != 50: Next code is incorrect
#endif
   
         BZ_ITAH(0);  BZ_ITAH(1);  BZ_ITAH(2);  BZ_ITAH(3);  BZ_ITAH(4);
         BZ_ITAH(5);  BZ_ITAH(6);  BZ_ITAH(7);  BZ_ITAH(8);  BZ_ITAH(9);
         BZ_ITAH(10); BZ_ITAH(11); BZ_ITAH(12); BZ_ITAH(13); BZ_ITAH(14);
         BZ_ITAH(15); BZ_ITAH(16); BZ_ITAH(17); BZ_ITAH(18); BZ_ITAH(19);
         BZ_ITAH(20); BZ_ITAH(21); BZ_ITAH(22); BZ_ITAH(23); BZ_ITAH(24);
         BZ_ITAH(25); BZ_ITAH(26); BZ_ITAH(27); BZ_ITAH(28); BZ_ITAH(29);
         BZ_ITAH(30); BZ_ITAH(31); BZ_ITAH(32); BZ_ITAH(33); BZ_ITAH(34);
         BZ_ITAH(35); BZ_ITAH(36); BZ_ITAH(37); BZ_ITAH(38); BZ_ITAH(39);
         BZ_ITAH(40); BZ_ITAH(41); BZ_ITAH(42); BZ_ITAH(43); BZ_ITAH(44);
         BZ_ITAH(45); BZ_ITAH(46); BZ_ITAH(47); BZ_ITAH(48); BZ_ITAH(49);
   
#        undef BZ_ITAH
   
         ge += BZ_G_SIZE;
         sel++;
      }

      ge -= BZ_G_SIZE;
      if (ge < emtfv) {
         Int32* code_row;
         AssertH ( *sel < nGroups, 3006 );
         code_row = &(ftab->code[*sel][0]);

	 //--- slow version for last group ---
         for (; ge < emtfv; ge++) {
            bsWc(code_row[*ge + sizeof(ftab->code) / sizeof(Int32)], code_row[*ge]);
         }
         sel++;
      }
   
      AssertH( sel - (UChar *)(s->arr2) == nSelectors, 3007 );
   
#if BZ_LIMIT_CODE_LEN == 17
      s->bsLive = bsLive;
      s->bsBuff = bsBuff;
      s->numZ = czbits - s->zbits;
#endif
   }

#undef bsWc

   if (s->verbosity >= 3)
      VPrintf1( "codes %d\n", s->numZ - nBytes );
}


/*---------------------------------------------------*/
void BZ2_compressBlock ( EState* s, Bool is_last_block )
{
   if (s->nblock > 0) {

      BZ_FINALISE_CRC ( s->blockCRC );
      s->combinedCRC = (s->combinedCRC << 1) | (s->combinedCRC >> 31);
      s->combinedCRC ^= s->blockCRC;
      if (s->blockNo > 1) s->numZ = 0;

      if (s->verbosity >= 2)
         VPrintf4( "    block %d: crc = 0x%08x, "
                   "combined CRC = 0x%08x, size = %d\n",
                   s->blockNo, s->blockCRC, s->combinedCRC, s->nblock );

      BZ2_blockSort ( s );
      generateMTFValues ( s );
   }

   /*-- If this is the first block, create the stream header. --*/
   if (s->blockNo == 1) {
      BZ2_bsInitWrite ( s );
      bsW( s, 24, (((BZ_HDR_B << 8) | BZ_HDR_Z) << 8) | BZ_HDR_h);
      bsW( s, 8, BZ_HDR_0 + s->blockSize100k );
   }

   if (s->nblock > 0) {

      bsW( s, 24, 0x314159);
      bsW( s, 24, 0x265359);

      /*-- Now the block's CRC, so it is in a known place. --*/
      bsPutUInt32 ( s, s->blockCRC );

      /*-- 
         Now a single bit indicating (non-)randomisation. 
         As of version 0.9.5, we use a better sorting algorithm
         which makes randomisation unnecessary.  So always set
         the randomised bit to 'no'.  Of course, the decoder
         still needs to be able to handle randomised blocks
         so as to maintain backwards compatibility with
         older versions of bzip2.
      --*/
      bsW(s,1,0);

      bsW ( s, 24, s->origPtr );
      sendMTFValues ( s );
   }


   /*-- If this is the last block, add the stream trailer. --*/
   if (is_last_block) {

      bsW( s, 24, 0x177245);
      bsW( s, 24, 0x385090);
      bsPutUInt32 ( s, s->combinedCRC );
      if (s->verbosity >= 2)
         VPrintf1( "    final combined CRC = 0x%08x\n   ", s->combinedCRC );
      bsFinishWrite ( s );
   }
}


/*-------------------------------------------------------------*/
/*--- end                                        compress.c ---*/
/*-------------------------------------------------------------*/
