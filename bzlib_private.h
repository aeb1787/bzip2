
/*-------------------------------------------------------------*/
/*--- Private header file for the library.                  ---*/
/*---                                       bzlib_private.h ---*/
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


#ifndef _BZLIB_PRIVATE_H
#define _BZLIB_PRIVATE_H

#include <stdlib.h>

#ifndef BZ_NO_STDIO
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#endif

#include "bzlib.h"



/*-- General stuff. --*/

#define BZ_VERSION  "1.0.6x, 16-Aug-2018"

typedef char            Char;
typedef unsigned char   Bool;
typedef unsigned char   UChar;
typedef int             Int32;
typedef unsigned int    UInt32;
typedef short           Int16;
typedef unsigned short  UInt16;

#ifdef __GNUC__
#ifdef __x86_64__
#define USE_INT64 1
typedef long long int          Int64;
typedef unsigned long long int UInt64;
#endif
#endif

#ifdef _MSC_VER
#ifdef _M_AMD64
#define USE_INT64 1
typedef long long int          Int64;
typedef unsigned long long int UInt64;
#endif
#endif

#ifdef USE_INT64
typedef Int64           Int;
typedef UInt64          UInt;
#else
typedef Int32           Int;
typedef UInt32          UInt;
#endif

#define True  ((Bool)1)
#define False ((Bool)0)

#ifndef __GNUC__
#define __inline__  /* */
#endif 

#ifndef BZ_NO_STDIO

extern void BZ2_bz__AssertH__fail ( int errcode );
#define AssertH(cond,errcode) \
   { if (!(cond)) BZ2_bz__AssertH__fail ( errcode ); }

#if BZ_DEBUG
#define AssertD(cond,msg) \
   { if (!(cond)) {       \
      fprintf ( stderr,   \
        "\n\nlibbzip2(debug build): internal error\n\t%s\n", msg );\
      exit(1); \
   }}
#else
#define AssertD(cond,msg) /* */
#endif

#define VPrintf0(zf) \
   fprintf(stderr,zf)
#define VPrintf1(zf,za1) \
   fprintf(stderr,zf,za1)
#define VPrintf2(zf,za1,za2) \
   fprintf(stderr,zf,za1,za2)
#define VPrintf3(zf,za1,za2,za3) \
   fprintf(stderr,zf,za1,za2,za3)
#define VPrintf4(zf,za1,za2,za3,za4) \
   fprintf(stderr,zf,za1,za2,za3,za4)
#define VPrintf5(zf,za1,za2,za3,za4,za5) \
   fprintf(stderr,zf,za1,za2,za3,za4,za5)

#else

extern void bz_internal_error ( int errcode );
#define AssertH(cond,errcode) \
   { if (!(cond)) bz_internal_error ( errcode ); }
#define AssertD(cond,msg)                do { } while (0)
#define VPrintf0(zf)                     do { } while (0)
#define VPrintf1(zf,za1)                 do { } while (0)
#define VPrintf2(zf,za1,za2)             do { } while (0)
#define VPrintf3(zf,za1,za2,za3)         do { } while (0)
#define VPrintf4(zf,za1,za2,za3,za4)     do { } while (0)
#define VPrintf5(zf,za1,za2,za3,za4,za5) do { } while (0)

#endif


#define BZALLOC(nnn) (strm->bzalloc)(strm->opaque,(nnn),1)
#define BZFREE(ppp)  (strm->bzfree)(strm->opaque,(ppp))


/*-- Header bytes. --*/

#define BZ_HDR_B 0x42   /* 'B' */
#define BZ_HDR_Z 0x5a   /* 'Z' */
#define BZ_HDR_h 0x68   /* 'h' */
#define BZ_HDR_0 0x30   /* '0' */
  
/*-- Constants for the back end. --*/

#define BZ_MAX_ALPHA_SIZE 258
#define BZ_LIMIT_CODE_LEN  17
#define BZ_MAX_DECODE_LEN  20
#define BZ_START_DECODE_LEN 8
#define BZ_MAX_CODE_LEN (BZ_MAX_DECODE_LEN + 3)

#define BZ_RUNA 0
#define BZ_RUNB 1

#define BZ_N_GROUPS 6
#define BZ_G_SIZE   50
#define BZ_N_ITERS  4

#define BZ_MAX_SELECTORS (2 + (900000 / BZ_G_SIZE))

/*-- Constants for the fast MTF coder / decoder. --*/

#define MTFA_SIZE 4096
#define MTFL_SIZE 16

/*-- Stuff for randomising repetitive blocks. --*/

extern Int32 BZ2_rNums[512];

#define BZ_RAND_DECLS                          \
   Int32 rNToGo;                               \
   Int32 rTPos                                 \

#define BZ_RAND_INIT_MASK                      \
   s->rNToGo = 0;                              \
   s->rTPos  = 0                               \

#define BZ_RAND_MASK ((s->rNToGo == 1) ? 1 : 0)

#define BZ_RAND_UPD_MASK                       \
   if (s->rNToGo == 0) {                       \
      s->rNToGo = BZ2_rNums[s->rTPos];         \
      s->rTPos++;                              \
      if (s->rTPos == 512) s->rTPos = 0;       \
   }                                           \
   s->rNToGo--;



/*-- Stuff for doing CRCs. --*/

extern UInt32 BZ2_crc32Table[256];

#define BZ_INITIALISE_CRC(crcVar)              \
{                                              \
   crcVar = 0xffffffffL;                       \
}

#define BZ_FINALISE_CRC(crcVar)                \
{                                              \
   crcVar = ~(crcVar);                         \
}

#define BZ_UPDATE_CRC(crcVar,cha)              \
{                                              \
   crcVar = (((UInt32)crcVar) << 8) ^          \
            BZ2_crc32Table[(crcVar >> 24) ^    \
                           cha];               \
}



/*-- States and modes for compression. --*/

#define BZ_M_IDLE      1
#define BZ_M_RUNNING   2
#define BZ_M_FLUSHING  3
#define BZ_M_FINISHING 4

#define BZ_S_OUTPUT    1
#define BZ_S_INPUT     2

#define BZ_N_RADIX 2
#define BZ_N_QSORT 12
#define BZ_N_SHELL 18
#define BZ_N_OVERSHOOT (BZ_N_RADIX + BZ_N_QSORT + BZ_N_SHELL + 2)




/*-- Structure holding all the compression-side stuff. --*/

typedef
   struct {
      /* pointer back to the struct bz_stream */
      bz_stream* strm;

      /* for doing the block sorting */
      UInt32*  arr1;
      UChar*   arr2;
      UInt32*  ftab;

      /* aliases for arr1 and arr2 */
      UInt32*  ptr;
      UChar*   block;
      UInt16*  mtfv;
      UChar*   zbits;

      /* the buffer for bit stream creation */
      UInt     bsBuff;
      Int32    bsLive;

      /* mode this stream is in, and whether inputting */
      /* or outputting data */
      Int32    mode;
      Int32    state;

      /* remembers avail_in when flush/finish requested */
      UInt32   avail_in_expect;

      /* for doing the block sorting */
      Int32    origPtr;

      /* for deciding when to use the fallback sorting algorithm */
      Int32    workFactor;

      /* run-length-encoding of the input */
      UInt32   state_in_ch;
      UInt32   state_in_len;
      BZ_RAND_DECLS;

      /* input and output limits and current posns */
      Int32    nblock;
      Int32    nblockMAX;
      Int32    numZ;
      Int32    state_out_pos;

      /* block and combined CRCs */
      UInt32   blockCRC;
      UInt32   combinedCRC;

      /* misc administratium */
      Int32    verbosity;
      Int32    blockNo;
      Int32    blockSize100k;

      /* stuff for coding the MTF values */
      Int32    nMTF;
      Int32    nInUse;

      Bool     inUse[256];
      Int32    mtfFreq    [BZ_MAX_ALPHA_SIZE];

   }
   EState;



/*-- externs for compression. --*/

extern void 
BZ2_blockSort ( EState* );

extern void 
BZ2_compressBlock ( EState*, Bool );

extern void 
BZ2_bsInitWrite ( EState* );

extern void 
BZ2_hbAssignCodes ( Int32*, Int32*, Int32, Int32, Int32 );

extern void 
BZ2_hbMakeCodeLengths ( Int32*, Int32*, Int32, Int32 );



/*-- states for decompression. --*/

#define BZ_X_IDLE        1
#define BZ_X_OUTPUT      2

#define BZ_X_MAGIC_0     10
#define BZ_X_MAGIC_1     11
#define BZ_X_MAGIC_2     12
#define BZ_X_MAGIC_3     13
#define BZ_X_MAGIC_4     14
#define BZ_X_BLKHDR_0    15
#define BZ_X_BLKHDR_1    16
#define BZ_X_BLKHDR_2    17
#define BZ_X_BLKHDR_3    18
#define BZ_X_BLKHDR_4    19
#define BZ_X_BLKHDR_5    20
#define BZ_X_BLKHDR_6    21
#define BZ_X_BCRC_1      22
#define BZ_X_BCRC_2      23
#define BZ_X_BCRC_3      24
#define BZ_X_BCRC_4      25
#define BZ_X_RANDBIT     26
#define BZ_X_ORIGPTR_1   27
#define BZ_X_ORIGPTR_2   28
#define BZ_X_ORIGPTR_3   29
#define BZ_X_MAPPING_1   30
#define BZ_X_MAPPING_2   31
#define BZ_X_SELECTOR_1  32
#define BZ_X_SELECTOR_2  33
#define BZ_X_SELECTOR_3  34
#define BZ_X_CODING_1    35
#define BZ_X_CODING_2    36
#define BZ_X_CODING_3    37
#define BZ_X_MTF_1       38
#define BZ_X_MTF_2       39
#define BZ_X_MTF_3       40
#define BZ_X_MTF_4       41
#define BZ_X_ENDHDR_2    42
#define BZ_X_ENDHDR_3    43
#define BZ_X_ENDHDR_4    44
#define BZ_X_ENDHDR_5    45
#define BZ_X_ENDHDR_6    46
#define BZ_X_CCRC_1      47
#define BZ_X_CCRC_2      48
#define BZ_X_CCRC_3      49
#define BZ_X_CCRC_4      50



/*-- Data to decode huffman codes of one group. --*/

typedef
   struct {
      Int32    limit  [BZ_MAX_CODE_LEN];
      Int32    base   [BZ_MAX_CODE_LEN];
      Int16    perm   [BZ_MAX_ALPHA_SIZE + (1 << BZ_START_DECODE_LEN)];
      UChar    len    [BZ_MAX_ALPHA_SIZE];
   }
   GRec;

/*-- Structure holding all the decompression-side stuff. --*/

typedef
   struct {
      /* pointer back to the struct bz_stream */
      bz_stream* strm;

      /* for undoing the Burrows-Wheeler transform (FAST) */
      UInt32   *tt;

      /* for undoing the Burrows-Wheeler transform (SMALL) */
      UInt16   *ll16;
      UChar    *ll4;

      /* state indicator for this stream */
      Int32    state;

      /* for doing the final run-length decoding */
      UChar    state_out_ch;
      Int32    state_out_len;
      Bool     blockRandomised;
      BZ_RAND_DECLS;

      /* save area for bit stream reading */
      UInt32   save_bsBuff;
      Int32    save_bsLive;

      /* misc administratium */
      Int32    blockSize100k;
      Bool     smallDecompress;
      Int32    currBlockNo;
      Int32    verbosity;

      /* for undoing the Burrows-Wheeler transform */
      Int32    origPtr;
      UInt32   tPos;
      Int32    k0;
      Int32    nblock_used;

      /* stored and calculated CRCs */
      UInt32   storedBlockCRC;
      UInt32   storedCombinedCRC;
      UInt32   calculatedBlockCRC;
      UInt32   calculatedCombinedCRC;

      void    *eblock;

      UChar    *mtfbase[256 / MTFL_SIZE];
      Int32    cftab[256];

      /* save area for scalars in the main decompress code */
      Int32    save_groupNo;
      Int32    save_groupPos;
      Int32    save_nblock;
      Int32    save_es;
      Int32    save_N;
      Int32    save_zn; 
      Int32    save_zvec;

      /* for decoding the MTF values */
      UChar    mtfa[MTFA_SIZE];
      UChar    selector[BZ_MAX_SELECTORS];
      GRec     groupRec[BZ_N_GROUPS];

   }
   DState;



/*-- Macros for decompression. --*/

#define BZ_GET_FAST(cccc)                    \
    cccc = (UChar)(s->tPos & 0xff);          \
    s->tPos = c_tt[s->tPos >> 8];            \

#define BZ_GET_FAST_C(cccc)                  \
    cccc = (UChar)(c_tPos & 0xff);           \
    c_tPos = c_tt[c_tPos >> 8];              \

#define GET_LL4(i)                           \
   ((((UInt32)(s->ll4[(i) >> 1])) >> (((i) << 2) & 0x4)) & 0xF)

#define GET_LL(i)                            \
   (((UInt32)c_ll16[i]) | (GET_LL4(i) << 16))

#define GET_LL4_C(i)                         \
   ((((UInt32)(c_ll4[(i) >> 1])) >> (((i) << 2) & 0x4)) & 0xF)

#define GET_LL_C(i)                          \
   (((UInt32)c_ll16[i]) | (GET_LL4_C(i) << 16))

#define SET_LL4_C(i,n)                                            \
   { if (((i) & 0x1) == 0)                                        \
        c_ll4[(i) >> 1] = (c_ll4[(i) >> 1] & 0xf0) | (n); else    \
        c_ll4[(i) >> 1] = (c_ll4[(i) >> 1] & 0x0f) | ((n) << 4);  \
   }

#define SET_LL_C(i,n)                        \
   { c_ll16[i] = (UInt16)(n & 0x0000ffff);   \
     SET_LL4_C(i, n >> 16);                  \
   }

#define BZ_GET_SMALL(cccc)                   \
    cccc = BZ2_indexIntoF ( s, s->tPos );    \
    s->tPos = GET_LL(s->tPos);

#define BZ_GET_SMALL_C(cccc)                 \
    cccc = BZ2_indexIntoF ( s, c_tPos );     \
    c_tPos = GET_LL_C(c_tPos);


/*-- externs for decompression. --*/

extern void
BZ2_indexCreateTab ( Int32*, UChar*, Int32 );

extern Int32 
BZ2_indexIntoF ( DState*, Int32 );

extern Int32 
BZ2_decompress ( DState* );

extern Int32
BZ2_hbCreateDecodeTables ( Int32*, Int32*, Int16*, UChar*, Int32 );


#endif


/*-- BZ_NO_STDIO seems to make NULL disappear on some platforms. --*/

#ifdef BZ_NO_STDIO
#ifndef NULL
#define NULL 0
#endif
#endif


/*-------------------------------------------------------------*/
/*--- end                                   bzlib_private.h ---*/
/*-------------------------------------------------------------*/
