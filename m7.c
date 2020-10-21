
/*********************************/

/****************************************
(C) Author: J.M. Montolio A. 2020/8.
****************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.H>
#include <assert.h>
#include <conio.h>
#include <time.h>
#include <omp.h> 

#define EQ ==
#define NE !=
#define GE >=
#define LE <=
#define GT >
#define LT <

#define FALSE 0
#define TRUE 1
#define AND &&
#define OR ||
#define NOT !

#define bitand &
#define bitor |
#define bitxor ^
#define bitnot ~

/****** PARAM ************/
#define DIMPAR 4
#define DCX    7900000
#define RVALX  2000000
#define XPEE 6
#define XPEf 720
#define XPBX 64
/****** modo 64b **********/
#define QFLX 65
#define SK64 64
#define XPWX 72
/**************************/
#define RULEMX 70
#define STRMX 40

typedef uint_fast16_t U16;
typedef int_fast32_t  I32;
typedef uint_fast32_t U32;
typedef uint_fast64_t U64;
typedef unsigned __int128 VAL128;
typedef uint_fast64_t VAL64;
typedef uint_fast64_t WCX;

#pragma pack(16)

static I32 iC=0;
static I32 PEf,RVAL=0;
static I32 QORBDUA=0;
static I32 QBASE=0;
static U16 biT=0,PEE=0,PBX;
static bool ESPOT=FALSE;
static int NUMTHR=DIMPAR;

static VAL64 PW2[XPWX]__attribute__((aligned));
static VAL64 APO[DCX];
static WCX CSET[DCX]__attribute__((aligned));
static bool ISDUPLI[DCX]__attribute__((aligned));
static I32 RSET[RVALX]__attribute__((aligned));
static I32 RBASE[RVALX];
static I32 INDQ[QFLX]__attribute__((aligned));
static I32 INDF[QFLX];
static I32 INDL[QFLX];
static I32 CENTER[DCX]__attribute__((aligned));
static I32 ORBSZ[DCX];
static I32 VALDUAL[DCX]__attribute__((aligned));
static I32 ORBDUAL[DCX];
static I32 QLE[DCX]__attribute__((aligned));
static I32 QGE[DCX];

static WCX WTF[DIMPAR][DCX]
__attribute__((aligned));
static I32 WHEIS[DIMPAR][DCX]
__attribute__((aligned));
static VAL64 PRODJK[DIMPAR][DCX];
static U64 RULNUM[RULEMX];
static char RULES[RULEMX][RULEMX];
static U16 PERMsmall[XPEf][XPEE];
static U16 PERMbig[XPEf][XPBX];
static char BUFCOM[STRMX];

/**************************************/
static inline void newl(void)__attribute__((cold));
static inline void newl(void)
{printf("\n");
}

static inline void afterR(void)
__attribute__((cold));
static inline void afterR(void)
{
printf("\r                           ");
printf("\r");
}

static inline void presskey(void)
__attribute__((cold));
static inline void presskey(void)
{
char c;
printf("-- (space) or doexit -- ");
c=getch();
if(c NE 32)exit(0);
newl();
}

static inline U16 istrlen(const char *p)
__attribute__((hot));
static inline U16 istrlen(const char *p)
{
U16 r=strlen(p);
return r;
}

static inline void today(void)
__attribute__((cold));
static inline void today(void)
{
time_t *a,b;

a=&b;b=time(a);
printf("  CLK %s \a\n",ctime(a));
}

static inline void pribin(const WCX p)
__attribute__((hot));
static inline void pribin(const WCX p)
{
WCX n;char ch;int L=biT-1;
n=p;
for(int b=L;b GE 0;b--)
  {
  ch='0';if((n&PW2[b]) NE 0)ch='1';
  printf("%c",ch);
  if(b%4 EQ 0)printf(" ");
  }
printf(";");
}


static inline void BUFROM128(const VAL128 p)
__attribute__((cold));
static inline void BUFROM128(const VAL128 p)
{
VAL128 n;
const char CHZ='0';
int k,m;

n=p;
for(int i=0;i<STRMX;i++)BUFCOM[i]=CHZ;
k=(STRMX-1);
while(n NE 0)
  {
  assert(k GE 2);
  m=n%10;
  if(k%5 EQ 0)BUFCOM[k--]='.';
  BUFCOM[k--]=(CHZ+m);
  n=(n-m);n=n/10;
  }
assert(k GE 2);
for(int i=0;i<k;i++)BUFCOM[i]=' ';
}

static void configOMP(void)__attribute__((cold));
static void configOMP(void)
{
int nuc;

nuc=omp_get_num_procs();
if(nuc GE DIMPAR)NUMTHR=DIMPAR;
else NUMTHR=nuc;
newl();
printf("OpenMP.    DIMPAR     : %4d \n",DIMPAR);
printf("  Detected NodeCores  : %4d \n",nuc);
printf("  Setting  Num.Threads: %4d \n",NUMTHR);
newl();

omp_set_num_threads(NUMTHR);
omp_set_dynamic(0);
omp_set_nested(0);
omp_set_max_active_levels( 1 );
omp_set_schedule(omp_sched_static,1000);
}

static void genrules(const int pn)
__attribute__((cold));
static void genrules(const int pn)
{
const int dkn=pn;
char l[dkn+4];
int a,b,k;

assert(dkn LE RULEMX);
a=3;b=3;
for(int r=0;r< dkn;r++)
  {
  if(r EQ 0){a=3;b=3;}
  for(int i=0;i<dkn;i++)l[i]='n';
  l[dkn]='f';l[dkn+1]=0;k=dkn-1-r;
  for(int i=a;i<=b;i++)
    {
    if(a EQ (a&i))l[k]='1';else l[k]='x';
    k++;
    }
    l[dkn-1-r]='e';
    l[dkn+1]=0;
    strcpy(RULES[r],l);
    RULES[0][dkn-1]='e';
    {
    WCX rvn;char ch;
    int cox,uder;

    strcpy(l,RULES[r]);uder=0;
    for(int b=dkn-1;
      (l[b] EQ '1') AND (b>=0);b--)uder++;

    rvn=0;cox=0;
    for(int b=dkn;b>=0;b--)
      {
      assert(b LE RULEMX);ch=l[b];
      if(ch EQ '1')rvn++;
      else if(ch EQ 'x')cox++;
      else if(ch EQ 'e'){rvn++;break;}
      rvn<<=1;
      }
    RULNUM[r]=rvn;
    }
  a--;if(a<0){a=b;b=2*b+1;}
  }
}

static bool createset(const int pn)
__attribute__((cold));
static bool createset(const int pn)
{
WCX tmp,rnm,pwm,mit;
I32 grovI,kC;
const int nk=pn;

assert(nk GE 1);
CSET[0]=0;CSET[1]=1;grovI=2;kC=grovI;biT=1;
for(int k=biT;k<nk;k++)
  {
  biT++;pwm=PW2[k];rnm=RULNUM[k];
  for(int i=1;i<grovI;i++)
    {
    tmp=CSET[i];tmp=tmp|pwm;mit=tmp&rnm;
    if(mit NE rnm)continue;
    if(kC GE DCX)return FALSE;
    CSET[kC++]=tmp;
    }
  grovI=kC;
  }
iC=grovI;
for(int i=0;i LE (iC-2);i++)
  assert(CSET[i] LT CSET[i+1]);
/* Checked set m.ascending */
return TRUE;
}

/******** hot functions **************/
static inline bool nbiti(
const WCX n,const U16 bit)__attribute__((hot));
static inline bool nbiti(
const WCX n,const U16 bit)
{
bool v;WCX w;

v=FALSE;
w=n&PW2[bit];if(w>0)v=TRUE;
return v;
}


static int factorial(const int e)
__attribute__((cold));
static int factorial(const int e)
{
int f=1;
for(int i=2;i<=e;i++)f*=i;
return f;
}

/* modo 64b */
#define TSUXEN 4
typedef struct
  {U16 E[TSUXEN];
  } TSUXE;
typedef union
  {WCX L;TSUXE P;  
  } TU64;

/*** mode WCX 64bit ***/
static inline bool CLEC(
const I32 pa,const I32 pb)
__attribute__((hot)) __attribute__((pure));
static inline bool CLEC(
const I32 pa,const I32 pb)
{
TU64 ua,ub,tmp;
/* 64b mode leq binary */

ua.L=CSET[pa];ub.L=CSET[pb];
for(int i=0;i LT TSUXEN;i++)
  {
  tmp.P.E[i]=(ua.P.E[i] & ub.P.E[i]);
  if(tmp.P.E[i] NE ua.P.E[i])
    return FALSE;
  }
return TRUE;
}

static void fillQLE(
const I32 ini,const I32 las)
__attribute__((hot));
static void fillQLE(
const I32 ini,const I32 las)
{

for(I32 ii=ini;ii<las;ii++)
  {
  const I32 igk=RSET[ii];
  I32 cc=0;

  cc=0;
  for(I32 k=0;k LE igk;k++)if(CLEC(k,igk))cc++;
  QLE[igk]=cc;

  if(ii%1000 EQ 0)printf("\rQ %09d",ii);
  }
afterR();
}


static void fillQGE(
const I32 ini,const I32 las)
__attribute__((hot));
static void fillQGE(
const I32 ini,const I32 las)
{
const I32 iCk=iC;

for(I32 ii=ini;ii<las;ii++)
  {
  const I32 igk=RSET[ii];
  I32 cc=0;

  cc=0;
  for(I32 k=igk;k LT iCk;k++)if(CLEC(igk,k))cc++;
  QGE[igk]=cc;

  if((ii%1000) EQ 0)printf("\rQ %09d",ii);
  }
afterR();
}

static void fillQQ(void)__attribute__((hot));
static void fillQQ(void)
{
const I32 iCk=iC;
const I32 RVALk=RVAL;
I32 q1,q2,q3;

/* Q precalculated.  ********************/
for(I32 i=0;i<DCX;i++){QLE[i]=0;QGE[i]=0;}
q1=RVALk/4;q2=2*q1;q3=3*q1;

#pragma omp parallel sections
{
#pragma omp section
  fillQLE(0,q1);

#pragma omp section
  fillQLE(q1,q2);

#pragma omp section
  fillQLE(q2,q3);

#pragma omp section
  fillQLE(q3,RVALk);

#pragma omp section
  fillQGE(0,q1);

#pragma omp section
  fillQGE(q1,q2);

#pragma omp section
  fillQGE(q2,q3);

#pragma omp section
  fillQGE(q3,RVALk);
}

/* Propiedad probada. */
for(I32 i=0;i<iCk;i++)
  {
  if(ISDUPLI[i])
    {QLE[i]=QLE[CENTER[i]];QGE[i]=QGE[CENTER[i]];
    }
  else
    {assert(QLE[i] GE 1);assert(QGE[i] GE 1);
    }
  }
}


static inline I32 findCi(
const int capoin,const bool mbe,const I32 ini,
const I32 fin,const U64 wtf)
__attribute__((hot))
__attribute__((pure));
static inline I32 findCi(
const int capoin,const bool mbe,const I32 ini,
const I32 fin,const WCX wtf)
{
const WCX wtfk=wtf;
const int k64=SK64;
I32 fi,la,mi,di;
int lc;

/* binary search */
if(wtfk EQ 0)lc=k64;
else lc=__builtin_clzll(wtfk);
fi=INDF[lc];if(ini GT fi)fi=ini;
la=INDL[lc];if(fin LT la)la=fin;
di=la-fi;
while(di GE 0) 
  {
  WCX Cmi;

  di>>=1; 
  mi=fi+di;Cmi=CSET[mi];
  if((wtfk LT Cmi))la=mi-1;
  else if((wtfk GT Cmi))fi=mi+1;    
  else return mi;
  di=(la-fi);
  }
if(mbe)assert(mbe EQ FALSE);
if(capoin EQ capoin)return -1;/* unused */
return -1;
}


static inline VAL64 unapo(const I32 ii)
__attribute__((hot));
static inline VAL64 unapo(const I32 ii)
{
const I32 iCk=iC;
const I32 iCk1=iC-1;
const int ntk=omp_get_thread_num();
const I32 iik=ii;
const WCX Cik=CSET[ii];
const WCX Cik4[4]={Cik,Cik,Cik,Cik};
VAL64 sapo;
const I32 L4=(iC-(iC%4));

for(I32 h=0;h<iCk;h++)WHEIS[ntk][h]=-1;

for(I32 h=0;h<L4;h+=4)
  {
  const I32 hk=h;

#pragma omp simd simdlen(4)
  for(I32 i=0;i< 4;i++)
    {WTF[ntk][hk+i]=(Cik4[i] bitand CSET[hk+i]);}
  }

for(I32 h=L4;h<iCk;h++)
  {WTF[ntk][h]=Cik bitand CSET[h];
  }

for(I32 h=0;h<iCk;h++)
  {
  I32 pista;

  pista=iik;if(h LT pista)pista=h;
  WHEIS[ntk][h]
  =findCi(1,TRUE,0,pista,WTF[ntk][h]);
  }

for(I32 h=0;h<iCk;h++)
  PRODJK[ntk][h]=QLE[WHEIS[ntk][h]];

for(I32 h=0;h<iCk;h++)WHEIS[ntk][h]=-1;

for(I32 h=0;h<L4;h+=4)
  {
  const I32 hk=h;

#pragma omp simd simdlen(4)
  for(I32 i=0;i< 4;i++)
    {WTF[ntk][hk+i]=(Cik4[i] bitor CSET[hk+i]);}
  }

for(I32 h=L4;h<iCk;h++)
  {WTF[ntk][h]=Cik bitor CSET[h];
  }

for(I32 h=0;h<iCk;h++)
  {
  I32 pista;

  pista=iik;if(h GT pista)pista=h;
  WHEIS[ntk][h]=
  findCi(2,TRUE,pista,iCk1,WTF[ntk][h]);
  }

for(I32 h=0;h<iCk;h++)
  {PRODJK[ntk][h]*=QGE[WHEIS[ntk][h]];}

sapo=0;
for(I32 h=0;h<iCk;h++)sapo+=PRODJK[ntk][h];
assert(sapo GT 0);return sapo;
}

static VAL128 x4Cn(void)__attribute__((hot));
static VAL128 x4Cn(void)
{
const I32 iCk=iC;
const I32 RVALk=RVAL;
VAL128 apo8;
I32 lcrum;

for(I32 i=0;i<DCX;i++)APO[i]=0;
for(I32 i=0;i<RVALX;i++)RBASE[i]=0;
lcrum=0;
for(I32 j=0;j<RVALk;j++)
  {
  I32 ii;
  ii=RSET[j];
  if(ORBDUAL[ii] LT 0)RBASE[lcrum++]=ii;
  }
QBASE=lcrum;
printf("DBG:  x4 ESPOT:%d QBASE= %d \n"
,ESPOT,QBASE);
if(!ESPOT)assert(QBASE EQ RVAL);

/* The HOT code */
  {
  const I32 QBASEk=QBASE;

#pragma omp parallel for schedule(static,1000) 
  for(I32 j=0;j<QBASEk;j++)
    {
    const I32 jk=j;
    const I32 iik=RBASE[jk];

    APO[iik]=unapo(iik);

#pragma omp critical
    if((jk%100) EQ 0)printf("\rx4 %06d",jk);
    }
  }
afterR();today();
if(!ESPOT)assert(QORBDUA EQ 0);

#pragma omp parallel for
for(I32 j=0;j<RVALk;j++)
  {
  I32 ii;

  ii=RSET[j];
  if(ORBDUAL[ii] GE 0)
    APO[ii]=APO[ORBDUAL[ii]];
  }

if(ESPOT)
  {
  for(int i=0;i<iCk;i++)
    {
    const I32 ik=i;
    const VAL64 APOk=APO[i];

    if(APOk EQ 0)continue;
#pragma omp parallel for 
    for(int j=ik+1;j<iCk;j++)
      {
      if(APO[j] NE APOk)continue;
      if(ORBDUAL[j] EQ ik)continue;
      printf("INTERNALERROR! APOeq: %d=%d %I64d \n"
        ,ik,j,APOk);
      exit(0);
      }
    }
  }

apo8=0;
for(I32 j=0;j<RVALk;j++)
  {
  VAL128 poa;I32 ig;
  
  ig=RSET[j];
  assert(APO[ig] GT 0);
  poa=(VAL128)(ORBSZ[ig]*APO[ig]);
  apo8+=poa;
  }

return apo8;
}


static void fillQFL(void)__attribute__((hot));
static void fillQFL(void)
{
const I32 iCk=iC;
const I32 iCk1=iC-1;
const int k64=SK64;

for(int i=0;i LT QFLX;i++)
  {INDQ[i]=0;INDF[i]=0;INDL[i]=0;
  }

for(I32 i=0;i LT iCk;i++)
  {
  WCX Ci=CSET[i];
  int lc;

  if(Ci EQ 0)lc=k64;
  else lc=__builtin_clzll(CSET[i]);

  assert(lc GE 0);assert(lc LE k64);
  if(INDQ[lc] EQ 0)INDF[lc]=i;
  INDQ[lc]++;
  }

for(int i=0;i<QFLX;i++)
  {INDL[i]=INDF[i]+INDQ[i]-1;}

  {
  I32 s=0;
  for(int i=0;i<QFLX;i++)s+=INDQ[i];
  assert(s EQ iCk);

  /* Check index integrity */
  s=biT;
  for(I32 i=0;i<iCk;i+=s)
    assert(findCi(0,FALSE,0,iCk1,CSET[i]) EQ i);
  }
}

static void fillDUA(void)__attribute__((hot));
static void fillDUA(void)
{
const int iCk=iC;
const int iCk1=iC-1;
const int biTk=biT;
const int biTk1=biT-1;
bool ba[SK64];

for(int i=0;i<iCk;i++)VALDUAL[i]=-1;
if(!ESPOT)return;

for(int i=0;i<RVAL;i++)
  {
  WCX w;I32 j;

  j=RSET[i];if(ISDUPLI[j])continue;
  w=CSET[j];
  for(int b=0;b<SK64;b++)ba[b]=0;
  for(int b=0;b<biTk;b++)
    {
    WCX f;bool e;
    e=FALSE;f=w&PW2[b];if(f>0)e=TRUE;
    ba[biTk1-b]= NOT e ;
    }
  w=0;/* ~w   */
  for(int b=0;b<biTk;b++)w+=ba[b]*PW2[b];
  VALDUAL[j]=findCi(0,FALSE,0,iCk1,w);
  printf("\rDUA %9d",j);
  }
afterR();
}

static void fillORBDUA(void)
__attribute__((hot));
static void fillORBDUA(void)
{
const I32 iCk=iC;

QORBDUA=0;
for(I32 i=0;i<iCk;i++)ORBDUAL[i]=-1;
if(!ESPOT)return;

for(I32 i=0;i<iCk;i++)
  {
  I32 di;

  if(ISDUPLI[i])continue;
  if(ORBDUAL[i] NE -1)continue;
  di=VALDUAL[i];if(di EQ i)continue;
  if(CENTER[di] EQ i)continue;
  if(ORBSZ[di] NE ORBSZ[i])continue;
  if(ORBDUAL[CENTER[di]] NE -1)continue;
  ORBDUAL[CENTER[di]]=i;QORBDUA++;
  }
}

static bool fillR(const int pn)
__attribute__((cold));
static bool fillR(const int pn)
{
I32  GLBPPX[XPEf];
const I32 iCk=iC;
const I32 iCk1=iC-1;
const int pnk=pn;
const int PEfk=PEf;
const int PBXk=PBX;
int impa,npermut=0;
bool BITPBX[XPBX];
bool GLBNUL[XPEf];
bool PEVALE[XPEf];

impa=pn;ESPOT=FALSE;
while(impa>1 AND (impa%2) EQ 0)impa/=2;
if(impa EQ 1)ESPOT=TRUE;

for(I32 i=0;i<DCX;i++)
  {ISDUPLI[i]=FALSE;ORBSZ[i]=-1;CENTER[i]=-1;
  }
for(U16 i=0;i<PEfk;i++)
  PEVALE[i]=TRUE;
if(!ESPOT)
 for(U16 p=0;p<PEfk;p++)
  {
  bool vale=TRUE;

  for(U16 b=0;b<pn;b++)
    {U16 d;d=PERMbig[p][b];
    if(d GE pnk)vale=FALSE;
    }
  PEVALE[p]=vale;
  }

npermut=0;
for(U16 p=0;p<PEfk;p++)if(PEVALE[p])npermut++;
printf("DBG: npermut %d \n",npermut);
assert(npermut GT 0);

/* Big for fillR */
for(I32 iii=0;iii<iCk;iii++)
  {
  I32 sze;WCX Ciii;

  if(ISDUPLI[iii])continue;
  CENTER[iii]=iii;Ciii=CSET[iii];

  for(U16 i=0;i<PEfk;i++)
    {GLBPPX[i]=-1;GLBNUL[i]=FALSE;
    }

/* for LLenar vectores */
  for(U16 p=0;p<PEfk;p++)
    {
    WCX Cp=0;

    if(!ESPOT AND !PEVALE[p])
      {
      GLBNUL[p]=TRUE;GLBPPX[p]=-1;
      continue;
      }

    for(U16 b=0;b<PBXk;b++)
      {U16 d;
      d=PERMbig[p][b];
      BITPBX[d]=nbiti(Ciii,b);
      }
    for(U16 b=0;b<PBXk;b++)Cp+=BITPBX[b]*PW2[b];
    GLBPPX[p]=findCi(0,FALSE,0,iCk1,Cp);
    if(GLBPPX[p] LT 0)GLBNUL[p]=TRUE;
    } /** Final llenar **/

  sze=0;
  for(U16 p=0;p<PEfk;p++)
    {
    I32 Cp;
 
    if(GLBNUL[p])continue;
    Cp=GLBPPX[p];sze++;
    for(U16 q=(p+1);q<PEfk;q++)
      if(GLBPPX[q] EQ Cp)GLBNUL[q]=TRUE;
    }

  assert(sze GE 1);assert(PEf%sze EQ 0);
  ORBSZ[iii]=sze;

  for(I32 p=0;p<PEfk;p++)
    {
    I32 j;

    if(GLBNUL[p])continue;
    j=GLBPPX[p];if(j EQ iii)continue;
    ISDUPLI[j]=TRUE;CENTER[j]=iii;
    ORBSZ[j]=sze;
    }
  if((iii%1000) EQ 0)printf("\rR:%9d",iii);  
  }

  afterR();/* fuera gran for */

  {
  I32 sorb,ctd;

  for(I32 i=0;i<RVALX;i++)RSET[i]=0;
  ctd=0;
  for(I32 i=0;i<iCk;i++)
    {if(!ISDUPLI[i])
      {
      if(ctd GE RVALX)return FALSE;
      RSET[ctd++]=i;
      }
    }
  RVAL=ctd;sorb=0;
  for(I32 i=0;i<iCk;i++)
    {if(!ISDUPLI[i])sorb+=ORBSZ[i];}
  assert(sorb EQ iCk);
  }
return TRUE;
}

static void fillPP(const int e)
__attribute__((cold));
static void fillPP(const int e)
{
U16 PPcifra[XPEE],BITPBX[XPBX];
U16 d,dig;U16 perm;
bool ok;I32 s,t;

PEE=e;PEf=factorial(e);PBX=PW2[e];
assert(biT LE PBX);assert(PBX LE XPBX);
assert(PEE LE XPEE);assert(PEf LE XPEf);
perm=0;
for(I32 n=0;perm LT PEf;n++)
  {
  t=n;ok=TRUE;
  for(U16 c=0;c<PEE;c++)
    {dig=t%10;
    if(dig GE PEE){ok=FALSE;break;}
    PPcifra[c]=dig;t/=10;
    }
  if(!ok)continue;

  for(U16 c=0;c<PEE;c++)
   for(U16 f=(c+1);f<PEE;f++)
     if(PPcifra[c] EQ PPcifra[f])
       {ok=FALSE;break;
       }
  if(!ok)continue;

  for(U16 c=0;c<PEE;c++)PERMsmall[perm][c]=PPcifra[c];
  perm++;
  }
assert(perm EQ PEf);
for(U16 p=0;p<PEf;p++)
  {
  for(U16 n=0;n<PBX;n++)
    {PERMbig[p][n]=0;
    for(U16 b=0;b<PEE;b++)
      {d=PERMsmall[p][b];BITPBX[d]=nbiti(n,b);
      }
    s=0;
    for(U16 b=0;b<PEE;b++)if(BITPBX[b])s+=PW2[b];
    PERMbig[p][n]=s;
    }
  }
}    

static void titles(void)
{
printf("Series OEIS A132581. Dedekind Numbers. \n");
printf("Author (C) JM Montolio A. Spain. 2020. ");
newl();
}

int main(void)
{
VAL128 FF;
int e,n4;
FILE *FLOG;
bool st;

titles();

for(int e=0;e<64;e++)
  {
/* dont touch PW2 dim */
  PW2[e]=1;
  for(int i=0;i<e;i++)PW2[e]*=2;;
  }

FLOG=fopen("DEDEKIND.log","a");
e= 6;   

for(int n=32;n LE 32;n++)
  {
  assert((unsigned)n LE PW2[e]);
  n4=4*n;
  printf("==== n %d, 4n %d ====\n",n,n4);
  configOMP();

  printf("PGCQFL -- n %2d -------\n",n);
  today();
  fillPP(e);genrules(n);
  st=createset(n);if(!st)continue;
  fillQFL();

  printf("RDO -- n %2d biT %2d iC %8d \n"
  ,n,biT,iC);
  today();
  st=fillR(n);if(!st)continue;
  fillDUA();
  fillORBDUA();

  printf(
"Q -- n %d biT %d iC %d RVAL %d DUAL %d:%d \n"
,n,biT,iC,RVAL,ESPOT,QORBDUA);
  today();

#pragma omp flush
  fillQQ();
#pragma omp flush

  printf(
"x4 -- n %d biT %d iC %d RVAL %d QDUA %d n4 %d \n"
  ,n,biT,iC,RVAL,QORBDUA,n4);
  today();

#pragma omp flush
  FF=x4Cn();
#pragma omp flush

  printf(
"ESPOT %d iC %d RVAL %d QORBDUA %d QBASE %d \n"
,ESPOT,iC,RVAL,QORBDUA,QBASE);

  BUFROM128(FF);
  printf("F(%4d)= %s \n\n",n4,BUFCOM);
  fprintf(FLOG,"%4d %s \n",n4,BUFCOM);
  fflush(FLOG);
  }

fflush(FLOG);
fclose(FLOG);
return 0;
}

