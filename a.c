

/*********************************/
/** Author JM Montolio A 2020  ***/
/** Computes OEIS A132581       **/
/** Related to Dedekind numbers **/
/** License MIT                 **/
/*********************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <stdint.h>
#include <limits.h>
#include <stdbool.h>
#include <string.H>
#include <assert.h>
#include <conio.h>
#include <math.h>
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

typedef int_fast8_t I8;
typedef uint_fast8_t U8;
typedef int_fast16_t  I16;
typedef uint_fast16_t U16;
typedef int_fast32_t  I32;
typedef uint_fast32_t U32;
typedef int_fast64_t I64;
typedef uint_fast64_t U64;
typedef __int128       I128;
typedef unsigned __int128 U128;

typedef unsigned __int128 VAL128;
typedef uint_fast64_t VAL64;
typedef uint_fast64_t WCX;

typedef struct
  {U32 W[2];
  } LEQUW;
typedef union
  {U64 L;LEQUW V;  
  } LEQLL;

#define TASKMX 64
#define SLENX 40
#define XPWX 72
static VAL64 PW2[XPWX];
static I32 iC;
static U8 wC;
static int LLEL;

#define DCX   8000000
static WCX C[DCX];
static I32 QLE[DCX];
static I32 QGE[DCX];

#define DRX 64
static U8   RUD[DRX];
static bool RTU[DRX];
static U8   RTX[DRX];
static WCX RN[DRX];
static char m[DRX][DRX];

#define DK32MX 2147483647
#define DK64MX 9223372036854775807

#define KOEX 81
#define X257 257
static VAL128 F[X257];
static bool K[X257];
static bool New[X257];
static VAL64 KOEVAL[KOEX]={1,3,6,14,20,53,84,148,168
,483,849,1681,2008,4414,5573,7413,7581
,22574,42825,92864,118462,286776,392247
,574114,595649,1563955,2314534,3710775
,3914177,6755518,7232705,7820773
,  7828354, 23477480, 46344318
,106040393,144293310,372009434,550589743
,884509066,937590104,2612894837,4200753533
,7401315524,8019492260,15265307010
,16893485613
,19009503981,19048569022
,55317412710,95000049832
,180535994846,201589169126
,415479721161,477614529319
,563195559435,565350642844,1243047970527
,1485569924436,1838302321392,1849331398154
,2377522493058,2395633471976,2414674212644
,2414682040998,7244038294639,14469012371130
,33692086389603,47448009054385
,125212050618508,194010762604008
,330477252598028,361255999828826
,1031737462826481,1752280502579682
,3305726594853050,3736260740451679
,7700396680569478,9007800679818571
,10859823582036583,10915524932266679};


U8 istrlen(char *p)
{
U8 r=strlen(p);
return r;
}

void newl(void)
{printf("\n");
}

void Dp(char *p)
{printf("%s. \n",p);
}

void retuline(int n)
{
printf("\r %d",n);
}

void clreol(void)
{
char buf[]="                       ";
printf("\r%s\r",buf);
}

int intmin(int a,int b)
{
int v;
v=a;if(b<a)v=b;
return v;
}

int intmax(int a,int b)
{
int v;
v=a;if(b>a)v=b;
return v;
}

void sprint128(char *buf,VAL128 p)
{
const int BX=SLENX;
const char CHZ='0';int k,m;
VAL128 n=p;

for(int i=0;i<BX;i++)buf[i]=CHZ;
k=BX-1;buf[k]=0;
while(n NE 0)
  {assert(k GE 1);
  m=n%10;buf[--k]=CHZ+m;n=(n-m);n=n/10;
  }
}

VAL128 str128(char *buf)
{
int sl,dig;
VAL128 n;
const char CHZ='0';

n=0;sl=istrlen(buf);
for(int k=0;k<sl;k++)
  {dig=buf[k]-CHZ;n*=10;n+=dig;
  }
return n;
}

void inits(void)
{
int v;

assert(DRX EQ 64);
for(int e=0;e LE 63;e++)
  {/* dont touch PW2 dim */
  PW2[e]=1;for(int i=0;i LT e;i++)PW2[e]*=2;;
  }

for(int n=0;n<X257;n++){K[n]=0;New[n]=0;}
for(int i=0;i<KOEX;i++)
  {K[2*i]=1;F[2*i]=(U128)KOEVAL[i];}

F[256]=str128("56130437228687557907788");
K[256]=TRUE;

omp_set_num_threads(1);
v=omp_get_num_procs();LLEL=v;
printf("OpenMP on N.Procs: %02d.",v);
printf("  Settings nonested, 1-level. \n");
omp_set_nested(0);omp_set_dynamic(0);
omp_set_schedule(omp_sched_static,2048);
omp_set_max_active_levels(1);
omp_set_num_threads(LLEL);
}

void DKgen(int pn)
{
char l[DRX];int dkn,a,b,k;

dkn=pn;
a=3;b=3;assert(pn LE DRX);
for(int r=0;r< dkn;r++)
  {
  assert(r LT DRX);
  if(r EQ 0){a=3;b=3;}
  for(int i=0;i<dkn;i++)l[i]='n';
  l[dkn]='f';l[dkn+1]=0;k=dkn-1-r;
  for(int i=a;i<=b;i++)
    {
    if(a EQ (a&i))l[k]='1';else l[k]='x';
    k++;
    }
  l[dkn-1-r]='e';l[dkn+1]=0;strcpy(m[r],l);
  a--;if(a<0){a=b;b=2*b+1;}
  }
  m[0][dkn-1]='e';
  {
  WCX rvn;char ch;

  for(int r=0;r<DRX;r++)RTU[r]=0;
  for(int r=0;r<DRX;r=2*r+1)RTU[r]=1;
  for(int r=0;r<dkn;r++)
    {
    int cox,uder;

    assert(r LE DRX);strcpy(l,m[r]);uder=0;
    for(int b=dkn-1;
      (l[b] EQ '1') AND (b>=0);b--)uder++;
    RUD[r]=uder;rvn=0;cox=0;
    for(int b=dkn;b>=0;b--)
      {
      assert(b LE DRX);ch=l[b];
      if(ch EQ '1')rvn++;
      else if(ch EQ 'x')cox++;
      else if(ch EQ 'e'){rvn++;break;}
      rvn<<=1;
      }
    RN[r]=rvn;RTX[r]=cox;
    }
  }
}

void DKexpanC(int pn)
{
WCX tmp,rnm,pwm,mit;
I32 lowC,kC;int n;

n=pn;assert(n GE 3);
DKgen(n);printf("Created %04d rules.",n);
C[0]=0;C[1]=1;C[2]=3;C[3]=5;C[4]=7;
lowC=5;kC=lowC;wC=3;
for(int m=wC;m<n;m++)
  {
  wC++;pwm=PW2[m];rnm=RN[m];
  for(int i=1;i<lowC;i++)
    {
    tmp=C[i];tmp=tmp|pwm;mit=tmp&rnm;
    if(mit NE rnm)continue;
    assert(kC LT DCX);C[kC++]=tmp;
    }
  lowC=kC;
  }
iC=lowC;
printf(" C(%d)=[%d*%d] \n",n,wC,iC);
if(n%2 EQ 0)
  {
  assert(iC EQ F[n]);
  }
}

bool aboutFn(int p,bool trz)
{
int nn=p;
bool some;

some=FALSE;
for(int e=1;e LE 8;e++)
  {
  int n,dx,a,b;

  n=PW2[e];
  if(e%2 EQ 0)dx=PW2[e/2];
  else        dx=PW2[(e-1)/2];
  for(int d=2;d LE dx;d*=2)
    {
    a=n-d;b=n-(n/d);

    if((nn EQ n) OR (nn EQ a) OR (nn EQ b))
      {
      some=TRUE;
      if(trz)printf("  F(%d)=F(%d)+F(%d). \n",n,a,b); 
      }
    }
  }
return some;
}


bool Dleq(WCX pa,WCX pb)
{
LEQLL ua,ub,und;/* leq binary */

assert(DRX EQ 64);
ua.L=pa;ub.L=pb;
und.V.W[0]=(ua.V.W[0] & ub.V.W[0]);
if(und.V.W[0] NE ua.V.W[0])return 0;
und.V.W[1]=(ua.V.W[1] & ub.V.W[1]);
if(und.V.W[1] NE ua.V.W[1])return 0;
return 1;
}

I32 DfindCi(I32 ini,I32 ult,WCX wtf)
{
WCX se,Cmi;
I32 n,nm1,fi,la,mi,di;
/* binary search */
se=wtf;n=iC;nm1=n-1;fi=ini;la=ult;
while(fi LE la) 
  {
  di=(la-fi);di=di>>1;mi=fi+di;Cmi=C[mi];
  if(se GT Cmi){fi=mi;fi++;}    
  else if(se LT Cmi){la=mi;la--;}
  else
    {assert(C[mi] EQ wtf);return mi;
    }
  }
assert(se GE C[0]);assert(se LE C[nm1]);
printf("EXIT ABORT INTERNAL DfindCi \n");
exit(0);
}

I32 DRskanup(I32 pp)
{
WCX vab;I32 p;I32 c;

c=0;p=pp;vab=C[p];
for(I32 k=0;k<=p;k++)
  {if(Dleq(C[k],vab))c++;
  }
return c;
}

I32 DRskandown(I32 pp)
{
WCX vab;I32 p;I32 c;

c=0;p=pp;vab=C[p];
for(I32 k=p;k LT iC;k++)
  {if(Dleq(vab,C[k]))c++;
  }
return c;
}

void fillQQ(void)
{
/* Q precalculated */

for(I32 i=0;i<DCX;i++){QLE[i]=0;QGE[i]=0;}
printf("\r Generating vector Q.");
printf("\r   scanning until card(C)");
#pragma omp parallel for schedule(static) 
for(I32 i=0;i<iC;i++)
  {
  I32 vl,vg;

  vl=DRskanup(i);vg=DRskandown(i);
#pragma omp critical
    {
    QLE[i]=vl;QGE[i]=vg;
    if(i%1000 EQ 0)retuline(i);
    }
  }

for(I32 i=0;i<iC;i++)
  {
  assert(QLE[i] GT 0);assert(QGE[i] GT 0);
  }
printf("\r Vector Q ready.");
clreol();
}


VAL64 DandorQQ(I32 pp)
{
WCX aa,bb,mitab,jonab;
VAL64 apo,rs,tmp;
I32 ii,mj,mk,n,nm1;
I32 rmit,sjon,iijmin,iijmax;

apo=0;ii=pp;aa=C[ii];n=iC;nm1=n-1;
/* era    0  */
for(I32 j=ii;j<iC;j++)
  {
  bb=C[j];
  iijmin=intmin(ii,j);
  iijmax=intmax(ii,j);
/* adjust search interval */
  mitab=(aa&bb);
  mj=DfindCi(0,iijmin,mitab);
  rmit=QLE[mj];

  jonab=(aa|bb);
  mk=DfindCi(iijmax,nm1,jonab);
  sjon=QGE[mk];

  rs=(VAL64)rmit;rs*=(VAL64)sjon;
  if(j NE ii)rs*=2;
/* no habia *2 */

  tmp=apo;apo+=rs;assert(apo GT tmp);
  }
return apo;
}

VAL128 Fx4(int n)
{
VAL128 FF;char buf[SLENX];
int n4=4*n;

DKexpanC(n);
fillQQ();

FF=0;
printf("\r  Adding up until card(C)");
#pragma omp parallel for schedule(static)
for(I32 i=0;i<iC;i++)
  {
  VAL64 P;
  P=DandorQQ(i);
#pragma omp critical
    {
    int t;
    FF+=P;t=omp_get_thread_num();
    if(t EQ 0)if(i%1000 EQ 0)retuline(i);
    }
  }

clreol();
printf("F(%4d)=",n4);
sprint128(buf,FF);printf("%s \n",buf);
if(F[n4] GT 0)assert(FF EQ F[n4]);
return FF;
}


void titles(void)
{
Dp("Compute OEIS A132581");
Dp("Author: JM Montolio A. 2020");
Dp("  gcc7+Openmp. MS Windows x64");
Dp("  MIT Licensed ");
newl();
Dp("Finds values chained to other two");
Dp("by the GENERAL RELATION. See Docs");
Dp("First the ruled set C is created");
Dp("After the standard x4 method is applied to it");
Dp("Note values above F(128) requires COMPUTERPOWER");
Dp("  If you have a Cray, no problem");
Dp("See Docs or Source for Details");
}

void endes(void)
{
Dp("Derived from Author works");
Dp("Contact: JOSMANCOLT@protonmail.com");
Dp("Ref's & Ack's");
Dp("Peter Koehler 2017");
newl();
Dp("Program ready to go F(256)");
Dp("If you reach it, let us know");
}

void presskey(void)
{
char c;
printf("\r  -- space or leave -- ");
c=getch();
clreol();
if(c NE 32)exit(0);
}

int main(void)
{
VAL128 FF;bool binded;int n4;
FILE *f;char buf[SLENX];

titles();endes();
presskey();newl();

inits();

f=fopen("harmor.log","w");
for(int n=4;n LE 64;n++)
  {
  n4=4*n;
  binded=aboutFn(n4,FALSE);
  if(binded)
    {
    clreol();
    newl();
    printf("Going for F(%d) because\n",n4);
    aboutFn(n4,TRUE);
    FF=Fx4(n);
    sprint128(buf,FF);
    fprintf(f,"%4d %s\n",n4,buf);
    if(0)presskey();
    }
  else 
    printf(
"\r Unrelated: %04d,",n4);
  }

endes();
return 1;
}

/************* MINGW BAT *****************
SET OMP_PROC_BIND="TRUE"
SET OMP_PLACES="CORES"
gcc -Wall -march=x86-64 -O2 -o %1.exe %1.c 
*****************************************/


