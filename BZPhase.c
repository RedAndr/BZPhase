
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glaux.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <WinCon.h>
#include <malloc.h>
#include <WinBase.h>
//#include <float.h>

#define pmax 1000000                                       // points count

#define n  11                                              // variables count
#define nk 15                                              // consts count
#define PI (3.14159265358979323846)

//#define double float

void CALLBACK myReshape(GLsizei w, GLsizei h);

int     NumPoints;                                        // begin conditions
double  BegRange,DeltaSolve,k[nk],x[n];                   // begin conditions
double  Delta,psize=1;
float   mx=0,my=0,mz=0;

int     pixcount,glmode=0;
int     dd=1,pl=1,axiz=1,persp=0;
double  zs275,xs2,ys2,zs2,zs22;
float   *pm,px[5][3]={{0,0,0},{1,1,0},{1,1,1},{1,0,1},{0,0,0}};

GLfloat roma[16] = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
};

#define p(i,j) (*(pm+i*4+j))

double ax=0, ay=0, az=0, maxcor;

#define inca 10
#define incm 0.02

GLsizei ow,oh;

GLfloat zx1[4] = { 0.9, 0.1, 0.1, 1.0 },
        zx2[4] = { 0.1, 0.1, 0.9, 1.0 },
        light_color1[4]=  { 1.0, 1.0, 1.0, 1.0 },
        light_color2[4]=  { 1.0, 1.0, 1.0, 1.0 },
        light_color3[4]=  { 1.0, 1.0, 1.0, 1.0 },
        light_color4[4]=  { 1.0, 1.0, 1.0, 1.0 },
//        mat_ambient1[4] = { 0.19225, 0.19225, 0.19225, 1.0 }, // Silver   //   { 2.0, 2.0, 2.0, 2.0 },
        mat_ambient1[4] = { 3.0, 3.0, 3.0, 1.0 },
        mat_diffuse1[4] = { 0.50754, 0.50754, 0.50754, 1.0 },
        mat_specular1[4]= { 0.508273,0.508273,0.508273,1.0 },             //   { 1.0, 1.0, 1.0, 1.0 },
//        mat_specular1[4]= { 0.51, 0.51, 0.51, 1.0 },
        lm_ambient1[4] =     { 0.2, 0.2, 0.2, 0.5 },
        light_position1[4]=  {  1,  1,  1, 1},
        light_position2[4]=  { -1, -1, -1, 1};
        light_position3[4]=  {  0,  0,  0, 0},
        light_position4[4]=  {  0,  0,  0, 1};

//  Initialize lighting and other values.
int myinit(void)
{

  light_position1[0]=-maxcor/2;  light_position1[1]=-maxcor/2;  light_position1[2]= maxcor;
  light_position2[0]= maxcor/2;  light_position2[1]= maxcor/2;  light_position2[2]= maxcor;
  light_position3[0]=-maxcor/2;  light_position3[1]= maxcor/2;  light_position3[2]= maxcor;
  light_position4[0]= maxcor/2;  light_position4[1]=-maxcor/2;  light_position4[2]= maxcor;
    
  glMaterialfv(GL_FRONT, GL_AMBIENT,  mat_ambient1);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,  mat_diffuse1);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular1);
  glMaterialf (GL_FRONT, GL_SHININESS, 51.2);

  glLightfv(GL_LIGHT0, GL_POSITION, light_position1);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_color1);
  glLightf (GL_LIGHT0, GL_SHININESS, 50.0);

  glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light_color2);
  glLightf (GL_LIGHT1, GL_SHININESS, 50.0);

  glLightfv(GL_LIGHT2, GL_POSITION, light_position3);
  glLightfv(GL_LIGHT2, GL_SPECULAR, light_color3);
  glLightf (GL_LIGHT2, GL_SHININESS, 50.0);

  glLightfv(GL_LIGHT3, GL_POSITION, light_position4);
  glLightfv(GL_LIGHT3, GL_SPECULAR, light_color4);
  glLightf (GL_LIGHT3, GL_SHININESS, 50.0);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lm_ambient1);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
  glEnable(GL_LIGHT3);

/*  glDepthFunc(GL_LESS);
  glEnable(GL_DEPTH_TEST);*/
  glEnable(GL_AUTO_NORMAL);

/*    glEnable (GL_LINE_SMOOTH);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);*/

  glLineWidth(psize);
  glPointSize(psize);

  return 0;
}

void rotxyz (int a, int b, int c, double theta, GLfloat A[4][4]);

void CALLBACK dis1    (void) {  maxcor/=1.2;   myReshape(ow,oh); }
void CALLBACK dis2    (void) {  maxcor*=1.2;   myReshape(ow,oh); }
void CALLBACK dis3    (void) {  persp=1-persp; myReshape(ow,oh); }
void CALLBACK psize1  (void) {  psize+=0.5;}
void CALLBACK psize2  (void) {  psize-=0.5;}
void CALLBACK glmodesw(void) {  glmode++; }
void CALLBACK move0   (void) {  mx=0; my=0; mz=0; ax=0; ay=0; az=0; }
void CALLBACK movex1  (void) {  mx+=incm;}
void CALLBACK movex2  (void) {  mx-=incm;}
void CALLBACK movey1  (void) {  my+=incm;}
void CALLBACK movey2  (void) {  my-=incm;}
void CALLBACK movez1  (void) {  mz+=incm;}
void CALLBACK movez2  (void) {  mz-=incm;}
void CALLBACK rotx1 (void) { rotxyz(1,0,0, inca,roma); }
void CALLBACK rotx2 (void) { rotxyz(1,0,0,-inca,roma); }
void CALLBACK roty1 (void) { rotxyz(0,1,0, inca,roma); }
void CALLBACK roty2 (void) { rotxyz(0,1,0,-inca,roma); }
void CALLBACK rotz1 (void) { rotxyz(0,0,1, inca,roma); }
void CALLBACK rotz2 (void) { rotxyz(0,0,1,-inca,roma); }


int MouseUp = 0;

static void/*GLenum*/ CALLBACK Mouse_leftup(AUX_EVENTREC *event) {
    MouseUp = 0;
//    return GL_FALSE;
}

GLint mouseX, mouseY, mouseS, mouseXo, mouseYo, mouseSo;

static void/*GLenum*/ CALLBACK Mouse_move(AUX_EVENTREC *event) {
      if(MouseUp==0){
        mouseX = event->data[AUX_MOUSEX];     
        mouseY = event->data[AUX_MOUSEY];     
        mouseS = event->data[AUX_MOUSESTATUS];
        MouseUp = 1;
      } else {
        mouseXo = mouseX;
        mouseYo = mouseY;
        mouseSo = mouseS;
        mouseX  = event->data[AUX_MOUSEX];
        mouseY  = event->data[AUX_MOUSEY];
        mouseS  = event->data[AUX_MOUSESTATUS];
        rotxyz(0,1,0, (mouseXo-mouseX),  roma);
        rotxyz(1,0,0, (mouseYo-mouseY),  roma);
      };

//      return GL_FALSE;
}


// Rozenbrok-type method of 3 order              }
// Written by Novikov E.A. on Fortran            }
// Translated by Ryzhkov A.B. for Borland Pascal }
// And then to C                                 }

int     hf=2,qf=20,nr,fun=0,jac=0,slt=0,stp=0,lum=0,mi[n],mt[n],ie,ls1,
        im=0,isa=0,ms[12]={99, 0, 0, 000, 0, 0, 0, 0, 0, 0, 0, 0};
double  h=1e-2,h2,t=0,tk,hm=1e-12,hmax=10,ep,tor,yy[n],ep1=0,
        jacb[n][n],jacbls[n][n],f[n],fnm[n],fn2[n],fn3[n],fn4[n],fn5[n],fn6[n],
        his,h2,h1,ah,sts,br,t1,cr,gr,p1,p2,p3,p4,bb1,bb2,pp1,h3,rt1,rt2,
        xflow,flowc0[n];

// the right-hand side of the ODE system
void rp(double *y, double *f){

  int i;
  double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15;

  w1  = k[ 0]*y[ 0]*y[ 2]  ;//  ! BrO3(-) + Br(-) + 2H(+) -> HBrO2 + HOBr               | 2.1      +1
  w2  = k[ 1]*y[ 1]*y[ 3]  ;//  ! HBrO2 + HOBr            -> BrO3(-) + Br(-) + 2H(+)    | 1.0E+4   -1
  w3  = k[ 2]*y[ 1]*y[ 2]  ;//  ! HBrO2 + Br(-) + H(+)    -> 2HOBr                      | 3.0E+6   +2
  w4  = k[ 3]*y[ 0]*y[ 1]  ;//  ! BrO3(-) + HBrO2 + H(+)  -> 2BrO2' + H2O               | 42       +3
  w5  = k[ 4]*y[ 4]*y[ 4]  ;//  ! 2BrO2' + H2O            -> BrO3(-) + HBrO2 + H(+)     | 4.2E+7   -3
  w6  = k[ 5]*y[ 4]*y[ 5]  ;//  ! BrO2' + Me(+) + H(+)    -> HBrO2 + Me(++)             | 8.0E+4   +4
  w7  = k[ 6]*y[ 1]*y[ 6]  ;//  ! HBrO2 + Me(++)          -> BrO2' + Me(+) + H(+)       | 8.9E+3   -4
  w8  = k[ 7]*y[ 1]*y[ 1]  ;//  ! 2HBrO2                  -> BrO3(-) + HOBr + H(+)      | 3.0E+3   +5
  w9  = k[ 8]*y[ 2]*y[ 3]  ;//  ! HOBr + Br(-) + H(+)     -> Br2 + H2O                  | 8.0E+9   +6
  w10 = k[ 9]*y[ 7]        ;//  ! Br2 + H2O               -> HOBr + Br(-) + H(+)        | 110      +7
  w11 = k[10]*y[ 7]*y[ 8]  ;//  ! RH + Br2                -> RBr + Br(-) + H(+)         | 4.6E-3   -7
  w12 = k[11]*y[ 3]*y[10]  ;//  ! HOBr + R'               -> ROH + Br'                  | 1E6      +8
  w13 = k[12]*y[ 8]*y[ 9]  ;//  ! RH + Br'                -> Br(-) + H(+) + R'          | 1E6      +9
  w14 = k[13]*y[ 6]*y[ 8]  ;//  ! RH + Me(++)             -> Me(+) + H(+) + R'          | 0.2      +10
  w15 = k[14]*y[10]*y[10]  ;//  ! 2R' + H2O               -> RH + ROH                   | 3.2E+9   +11
  
// w13=w12

  f[ 0] = -w1+w2-w4+w5+w8             ;// ! BrO3m
  f[ 1] = +w1-w2-w3-w4+w5+w6-w7-w8-w8 ;// ! HBrO2
  f[ 2] = -w1+w2-w3-w9+w10+w11+w13    ;// ! Brm
  f[ 3] = +w1-w2+w3+w3+w8-w9+w10-w12  ;// ! HOBr
  f[ 4] = +w4+w4-w5-w5-w6+w7          ;// ! BrO2'
  f[ 5] = -w6+w7+w14                  ;// ! Me+
  f[ 6] = +w6-w7-w14                  ;// ! Me++
  f[ 7] = +w9-w10-w11                 ;// ! Br2
  f[ 8] = -w11-w13-w14+w15            ;// ! RH
  f[ 9] = +w12-w13                    ;// ! Br'
  f[10] = -w12+w13+w14-w15-w15        ;// ! R'
  
  for (i=0;i<11;i++) {
    f[i] += xflow*(flowc0[i]-y[i]);
//    printf("y[%2u]=%E f[%2u]=%E\n",i,y[i],i,f[i]);
  }

}

//  Creation the matrix D = E - ahf' and calculation norm of matrix f'
void node04() {
int  i,j,k;
double ss1;
  sts=0;
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      jacb[i][j]=-ah*(jacbls[i][j]); 
      ss1=fabs(jacbls[i][j]); if(sts<ss1) sts=ss1;
    };
    jacb[i][i]++;
  };
}

/*{
c      The program NODE uses two procedures NODE18 and NODE19. They solve
c the double precision system Ax = b using the LU - factorization of
c the matrix A. The factorization can be written A = L * U. The calls
c to the procedures NODE18 and NODE19 have the form:
c
c            call NODE18(jacb, mi, ie)
c            call NODE19(jacb, fn2, mi)
c
c      The meanings of the parameters are:
c  jacb - the matrix to be factored of order (N,N);
c  mi   - work array of length N;
c  IE   - integer parameter. If IE = 0, the factorization succeedes;
c  fn2  - the right - hand side of the system Ax = b.
}*/
void node18(){

  int    i,j,k,m;
  double t;

  ms[9]++;
  ie=0;
  mi[n-1]=1;
  /*if (n!=1)*/ 
  {
    for( k=1; k<=n-1; k++ ) {

      m = k;

      for ( i=k; i<n; i++ ) 
         if ( fabs(jacb[i][k-1]) > fabs(jacb[m-1][k-1]) ) 
            m = i;

      mi[k-1]=m;

      t = jacb[m-1][k-1];

      if (m!=k) { 
         mi[n-1]=-mi[n-1]; 
         jacb[m-1][k-1]=jacb[k-1][k-1]; 
         jacb[k-1][k-1]=t;
      };

      if ( t==0 ) { 
         ie=k; 
         mi[n-1]=0; 
         return; 
      }

      t=1/t;

      for( i=k; i<n; i++)  
        jacb[i][k-1] = -jacb[i][k-1]*t;

      for(j=k;j<n;j++) {

          t=jacb[m-1][j];  
          jacb[m-1][j]=jacb[k-1][j];  
          jacb[k-1][j]=t;

          if(t!=0) 
             for ( i=k; i<n; i++)  
                jacb[i][j] += jacb[i][k-1]*t;
      };
    };
  };

  if ( jacb[n-1][n-1] == 0 ) { 
    ie = n; 
    mi[n-1] = 0; 
  };

}

void node19(double *b){
  
  int    k,m,i,kb,km1;
  double t;

  ms[10]++;
//  if(n==1) { b[0]/=jacb[0][0]; return;};
  
  for(k=1;k<=n-1;k++){
    m=mi[k-1];
    t=b[m-1]; b[m-1]=b[k-1]; b[k-1]=t;
    for(i=k;i<n;i++) b[i]+=jacb[i][k-1]*t;
  }

  for(kb=1;kb<=n-1;kb++){
    km1=n-kb;
    k=km1+1;
    b[k-1]/=jacb[k-1][k-1];
    t=-b[k-1];
    for (i=0;i<km1;i++) 
      b[i]+=jacb[i][k-1]*t;
  }

  b[0]/=jacb[0][0];
}

/* L-stable Rozenbrok-type method 3-order of precision
   with freezing of Jacobi matrix */
void node06(){

  int  i,j;
  double yn;

  if(im!=1){
    im=1;
    rp(yy,fnm); ms[7]++;
    for(i=0;i<n;i++) f[i]=1/(fabs(yy[i])+tor);
    isa=0;
    mt[5]=1;
    ep1=ep*2;
  };
  his=h;

  h2=tk-t;

  if(h>h2) h=h2;
  
  if(fabs(his-h)>1e-12) mt[5]=1;

  if (mt[5]!=0) {
l40:
    if(ms[3]==0) {
      for(i=0;i<n;i++){
        h1=1e-7*fabs(yy[i]);
        if (h1<1e-12) h1=1e-12;
        yn=yy[i];  yy[i]=yy[i]+h1;   rp(yy,fn2);
        for(j=0;j<n;j++) jacbls[j][i]=(fn2[j]-fnm[j])/h1;
        yy[i]=yn;
      };
      ms[7]+=n; 
      ms[8]++;
    } else /*drp()*/return;
    isa=0;
l90:
    do{
      ah=0.39735167201943117*h;
      node04();         //(ah, sts, jacbls, jacb)
      sts=sts*h;
      node18();         //(jacb, mi, ie)
      mt[5]=0;
      if(ie!=0) h*=0.8;
    } while (ie!=0);
  };
  for(i=0;i<n;i++){
    fn2[i]=fnm[i];
    fn6[i]=2.1100365925712e-2*fnm[i];
  };
  node19(fn2);         //(jacb, fn2, mi)
  br=h*0.79470334403886235;
  t1=t+br;
  for(i=0;i<n;i++){
    fn4[i]=yy[i]+br*fn2[i];
    fn6[i]=fn6[i]+9.5779926814858e-1*fn2[i];
  }
  rp(fn4, fn3);  ms[7]++;
  for(i=0;i<n;i++){
    fn6[i]=fn6[i]-fn3[i];
    fn3[i]=fn3[i]-0.97889963407428840*fn2[i];
  }
  node19(fn3);         //(jacb, fn2, mi)
  cr= h*0.23616433228820091;
  gr=-h*0.16118733973112458;
  for(i=0;i<n;i++){
    fn5[i]=yy[i]+cr*fn2[i]+gr*fn3[i];
    fn6[i]=fn6[i]+fn3[i];
  }
  t1=t+h*0.913461538460975e-1;
  rp(fn5, fn4);  ms[7]++;
  for(i=0;i<n;i++) fn4[i]=fn4[i]-0.10747713570105259e1*fn2[i]+0.15470514246623036*fn3[i];
  node19(fn4);         //(jacb, fn2, mi)
  for(i=0;i<n;i++) fn5[i]=fn4[i]+0.79957730742082261e-1*fn2[i]-0.40050058265625429*fn3[i];
  node19(fn5);         //(jacb, fn2, mi)
  rt1=0;
  p1=h*0.98991591230618822;
  p2=h*0.59256424028676605;
  p3=h*0.033832126509584816;
  p4=h*0.45560848220695990;
  rt2=0;
  for(i=0;i<n;i++){
    bb1=fabs(fn3[i]-2.1100365925712e-2*fn2[i]);
    bb2=0;
    if(bb1>1e-14) bb2=fabs(fn6[i])/bb1;
    if(rt2<bb2) rt2=bb2;
    bb1=-0.79957730742012588e-1*fn2[i]+0.40050058265586903*fn3[i]+fn5[i]-fn4[i];
    fn2[i]=yy[i]+p1*fn2[i]+p2*fn3[i]+p3*fn4[i]+p4*fn5[i];
    fn3[i]=bb1;
    bb1=fabs(bb1)*f[i];
    if(rt1<bb1) rt1=bb1;
  }
//  rt2=rt2/0.39735167201943117;
  rt2=sts;
  rt1=rt1*h;
  pp1=rt1;
  if(!((h<=hm)||(rt1<=ep1))){
    rt1=0;
    node19(fn3);                        //(jacb, fn2, mi)
    for(i=0;i<n;i++){
      bb1=fabs(fn3[i])*f[i];
      if(rt1<bb1) rt1=bb1;
    }
    rt1=rt1*h;
    if(!((h<=hm)||(rt1<ep1))){
      do{
        rt1=rt1/1.331;
        h=h/1.1;
      } while(rt1>=ep1);
      h=h/1.1;
      if(h<hm) h=hm;
      ms[11]++;
      if(isa==0) goto l90;
      goto l40;
    }
  }
  for(i=0;i<n;i++){
    f[i]=1/(fabs(fn2[i])+tor);
    yy[i]=fn2[i];
  }
  rp(yy, fnm); ms[7]++;
  t+=h;
  ms[6]++;
  isa++;
  h3=h;
  if(h3<hm) h3=hm;
  if(pp1<=ep1){
    for(;;){
      h3=1.1*h3;
      rt1=rt1*1.331;
      rt2=rt2*1.1;
      if(!((rt1<ep1)&&(h3<h2))) break;
    }
    h3=h3/1.1;
    if((h3<hf*h)&&(isa<qf)) return;
  }
  h=h3;
  if(h>hmax) h=hmax;
  mt[5]=1;
}

void stepx(){
  tk=t+Delta;
  do{
    node06();
  }while (fabs(t-tk)>1e-12);
  stp=ms[6];  fun=ms[7];  jac=ms[8];  lum=ms[9];  slt=ms[10];
}

int oldxx,oldyy;

GLfloat  bondmat[] = { 0.27, 0.57, 0.73, 1.0 };

void rotxyz (int a, int b, int c, double theta, GLfloat A[4][4])
{
      GLfloat    ct,st,B[9];

      ct = cos(PI*theta/180.0);
      st = sin(PI*theta/180.0);

      B[0] = A[0][0]*(a+ct*(1-a)) + A[0][1]*(  st*c)     + A[0][2]*( -st*b);
      B[1] = A[0][0]*( -st*c)     + A[0][1]*(b+ct*(1-b)) + A[0][2]*(  st*a);
      B[2] = A[0][0]*(  st*b)     + A[0][1]*( -st*a)     + A[0][2]*(c+ct*(1-c));

      B[3] = A[1][0]*(a+ct*(1-a)) + A[1][1]*(  st*c)     + A[1][2]*( -st*b);
      B[4] = A[1][0]*( -st*c)     + A[1][1]*(b+ct*(1-b)) + A[1][2]*(  st*a);
      B[5] = A[1][0]*(  st*b)     + A[1][1]*( -st*a)     + A[1][2]*(c+ct*(1-c));

      B[6] = A[2][0]*(a+ct*(1-a)) + A[2][1]*(  st*c)     + A[2][2]*( -st*b);
      B[7] = A[2][0]*( -st*c)     + A[2][1]*(b+ct*(1-b)) + A[2][2]*(  st*a);
      B[8] = A[2][0]*(  st*b)     + A[2][1]*( -st*a)     + A[2][2]*(c+ct*(1-c));

      A[0][0] = B[0];
      A[0][1] = B[1];
      A[0][2] = B[2];

      A[1][0] = B[3];
      A[1][1] = B[4];
      A[1][2] = B[5];

      A[2][0] = B[6];
      A[2][1] = B[7];
      A[2][2] = B[8];
}

void CALLBACK display(void)
{
   int c;
 
   glLineWidth(psize);
   glPointSize(psize);

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glShadeModel(GL_SMOOTH);

   glPushMatrix ();
   glMultMatrixf(roma);
 
        if (glmode==0) glBegin    (GL_POINTS);
   else if (glmode==1) glBegin    (GL_LINE_STRIP);
   else if (glmode==2) glBegin    (GL_LINES);
   /*else if (glmode==3) glBegin    (GL_TRIANGLES);
   else if (glmode==4) glBegin    (GL_TRIANGLE_STRIP);
   else if (glmode==5) glBegin    (GL_QUADS);
   else if (glmode==6) glBegin    (GL_QUAD_STRIP);
   else if (glmode==7) glBegin    (GL_POLYGON);*/
   else               {glBegin    (GL_POINTS); glmode=0;}

   bondmat[0]=0.5;   bondmat[1]=1-p(0,3);   bondmat[2]=p(0,3);
//   glMaterialfv(GL_FRONT, GL_DIFFUSE, bondmat);
   glMaterialfv(GL_FRONT, GL_AMBIENT,  bondmat);
   glMaterialfv(GL_FRONT, GL_DIFFUSE,  bondmat);
   glMaterialfv(GL_FRONT, GL_SPECULAR, bondmat);
   glMaterialf (GL_FRONT, GL_SHININESS, 51.2);
   if (glmode==7){
     glNormal3f (mx,my,mz);
     glVertex3f (mx,my,mz);
   } else {
     glNormal3f (p(0,0)+mx,p(0,1)+my,p(0,2)+mz);
     glVertex3f (p(0,0)+mx,p(0,1)+my,p(0,2)+mz);
   }

   glMaterialf (GL_FRONT, GL_SHININESS, 51.2);

   for (c=1;c<pixcount;c++) {
      bondmat[0]=0.8*4;   bondmat[1]=(1-p(c,3))*4;   bondmat[2]=p(c,3)*4;
      glMaterialfv(GL_FRONT, GL_AMBIENT,  bondmat);
      glMaterialfv(GL_FRONT, GL_DIFFUSE,  bondmat);
      glMaterialfv(GL_FRONT, GL_SPECULAR, bondmat);
//      glNormal3f (p(c,0)+mx-(p(c-1,0)+p(c+1,0))/2, p(c,1)+mx-(p(c-1,1)+p(c+1,1))/2, p(c,2)+mx-(p(c-1,2)+p(c+1,2))/2 );
      glVertex3f (p(c,0)+mx,p(c,1)+my,p(c,2)+mz);
//      glNormal3f (p(c,0)+mx,p(c,1)+my,p(c,2)+mz);
   }

   glEnd ();

   glPopMatrix ();

   glFlush();
   auxSwapBuffers();
}

void CALLBACK myReshape(GLsizei w, GLsizei h)
{
  ow=w; oh=h;
  h = (h == 0) ? 1 : h;
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  if (persp==0) {
    if (w <= h) glOrtho (-maxcor    , maxcor    , -maxcor*h/w, maxcor*h/w, -maxcor, maxcor);
    else        glOrtho (-maxcor*w/h, maxcor*w/h, -maxcor    , maxcor    , -maxcor, maxcor);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated (0.0, 0.0, 0.0);  // viewing transform
  } else {
    gluPerspective(60.0, (double)w/(double)h, 1.0, 20.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated (0.0, 0.0, -2*maxcor);  // viewing transform
  }
}

float minpx, maxpx, minpy, maxpy, minpz, maxpz, maxpw, minpw;
int cj = 0, mj,mmj;

void CALLBACK Calc() {

  stepx();

  p(cj,0)= (float)yy[2];                // Br-  2
  p(cj,1)= (float)yy[5];                // Me   5
  p(cj,2)= (float)yy[7];                // Br2  7
  p(cj,3)= (float)yy[9];                // Br'  9
  
  p(cj,0)=(p(cj,0)-minpx)/(maxpx-minpx)-0.5;
  p(cj,1)=(p(cj,1)-minpy)/(maxpy-minpy)-0.5;
  p(cj,2)=(p(cj,2)-minpz)/(maxpz-minpz)-0.5;
  p(cj,3)=(p(cj,3)-minpw)/(maxpw-minpw);

  cj++; 
  if (cj==mj) cj=0;

//  printf("\r %6.2f %%",t*100/mmj); fflush(stdout);

  if (cj % 100000 == 0) display();

}


void CALLBACK scale() {

   int c;
   float sx/*, minpx2, maxpx2, minpy2, maxpy2, minpz2, maxpz2, maxpw2, minpw2*/;

   for(c=0;c<pixcount;c++) {                                    // inverse normalize
     p(c,0)=(p(c,0)+0.5)*(maxpx-minpx)+minpx;
     p(c,1)=(p(c,1)+0.5)*(maxpy-minpy)+minpy;
     p(c,2)=(p(c,2)+0.5)*(maxpz-minpz)+minpz;
     p(c,3)=(p(c,3)    )*(maxpw-minpw)+minpw;
   }

   maxpx=minpx=p(0,0);                                          // calculate new min and max
   maxpy=minpy=p(0,1);
   maxpz=minpz=p(0,2);
   maxpw=minpw=p(0,3);

   for(c=1;c<pixcount;c++) {
           sx=p(c,0); if(sx>maxpx) maxpx=sx; else if(sx<minpx) minpx=sx;
           sx=p(c,1); if(sx>maxpy) maxpy=sx; else if(sx<minpy) minpy=sx;
           sx=p(c,2); if(sx>maxpz) maxpz=sx; else if(sx<minpz) minpz=sx;
           sx=p(c,3); if(sx>maxpw) maxpw=sx; else if(sx<minpw) minpw=sx;
   }

   for(c=0;c<pixcount;c++) {                                    // new normalize
     p(c,0)=(p(c,0)-minpx)/(maxpx-minpx)-0.5;
     p(c,1)=(p(c,1)-minpy)/(maxpy-minpy)-0.5;
     p(c,2)=(p(c,2)-minpz)/(maxpz-minpz)-0.5;
     p(c,3)=(p(c,3)-minpw)/(maxpw-minpw);
   }

}


void CALLBACK print() {
  int i;
  printf("%20.12E\n",t);
  for(i=0;i<n ;i++) printf("%20.12E\n",yy[i]);
}


void CALLBACK par1() {  xflow += 0.100000E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par2() {  xflow += 0.010000E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par3() {  xflow += 0.001000E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par4() {  xflow += 0.000100E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par5() {  xflow += 0.000010E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par6() {  xflow += 0.000001E-6; printf("Flow parameter = %14.8E\n", xflow); }

void CALLBACK par7() {  xflow -= 0.100000E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par8() {  xflow -= 0.010000E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK par9() {  xflow -= 0.001000E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK para() {  xflow -= 0.000100E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK parb() {  xflow -= 0.000010E-6; printf("Flow parameter = %14.8E\n", xflow); }
void CALLBACK parc() {  xflow -= 0.000001E-6; printf("Flow parameter = %14.8E\n", xflow); }

void help(){
  puts("\nUsage of BZPhase is 'BZPhase [DataFile.DAT]'");
  puts("if DataFile.DAT isn't specified then by default is used BZPhase.DAT");
  puts("\nWork keys:");
  puts("        Left,Right,Up,Down,z,x - to rotate attractor");
  puts("        1,2 - to scale");
  puts("        3 - to perspective");
  puts("        4 - to change view mode (pixel or line)");
  puts("        c,v - to change pixel or line size");
  puts("        a,d; s,w; q,e - to move attractor");
}


int main(int pn, char **ps){
  int i,j,pxc,c;
  double sx,second[4];
  FILE *ic;
  char ss[60];

  puts(" ###############################################################################");
  puts(" # BZPhaseFlow -  Phase Portraits Builder of the Belousov-Zhabotinsky reaction #");
  puts(" # Copyright (C) Andrew B. Ryzhkov and Arcady V. Antipin, 1997-2006. Ver. 2.00 #");
  puts(" # Ufa,      Institute of Organic Chemistry,   Laboratory of Chemical Kinetics #");
  puts(" # Montreal, McGill University, Department of Oceanic and Atmospheric Sciences #");
  puts(" # E-Mail: RedAndr@GMail.com                         WWW: http://RedAndr.ca/bz #");
  puts(" ###############################################################################");
                                                                                      
  SetConsoleTitle("BZPhaseFlow 2.00");                                                    
                                                                                      
  /* for Borland */
//  _clear87();
//  _control87(MCW_EM, MCW_EM);  /* defined in float.h */

  if(pn<2) {                                                                          
    ic=fopen("BZPhase.dat","rt");                                                     
    if(ic!=NULL) {
ReadData:
      do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&BegRange);
      do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%d" ,&NumPoints);
      do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&DeltaSolve);
      do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&ep); tor=ep;

      for(i=0;i<nk;i++) { do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&k[i]);}
      for(i=0;i<n ;i++) { do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&x[i]);}
      for(i=0;i<n ;i++) { do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&flowc0[i]);}

      printf("BegRange = %5.0f, NumPoints = %6d, DeltaSolve = %2.2f, Precision = %9.2E\n",  BegRange,NumPoints,DeltaSolve,ep);
      for(i=0;i<n ;i++) printf("k[%2d] = %E, c0[%2d] = %E, flowc0[%2d] = %E\n", i+1, k[i], i+1, x[i], i+1, flowc0[i]);
      for(i=n;i<nk;i++) printf("k[%2d] = %E\n", i+1, k[i]);

      do fgets(ss,80,ic); while(ss[0]==';'); sscanf(ss,"%le",&xflow);                           // flow parametr
      printf("Flow parameter = %14.8E\n", xflow);

      fclose(ic);
    } else if(errno==2) {
      help();
      return 1;
    } else {
      printf("Error opening file 'bzphase.dat' #%d\n",errno);
      return 2;
    }
  } else {
    if(pn>2) {
      puts("Opening file");
      ic=fopen(ps[1],"rt");
      if(ic==NULL) { printf("Error opening file '%s' #%d",ps[1],errno); return 2; }
      pm=malloc(pmax*sizeof(float)*4);
      if( pm == NULL ) {
        puts( "Unable to allocate memory\n" );
        return -1;
      }
      puts("Reading file");
      j=0;
      do {
        fscanf(ic,"%le %le %le %le",&t,&p(j,0),&p(j,1),&p(j,2));
        j++;
        if (j%100==0) printf("\r %d",j);
        fflush(stdout);
        if (j==pmax) break;
      } while (!feof(ic));
      fclose(ic);
      printf("\n%d points is readout\n",j-1);
      goto Show;
    } else {
      ic=fopen(ps[1],"rt");
      if(ic!=NULL) goto ReadData;
      printf("Error opening file '%s' #%d\n",ps[1],errno);
      return 2;
    }
  }

  tor=ep;
  for(i=0;i<n;i++) yy[i]=x[i];          // begin concs

  printf("Calculate begining range: %10.2f",BegRange);  fflush(stdout);
  Delta=BegRange;
  second[0]=(double)(GetTickCount())/1000;
  stepx();
  second[1]=(double)(GetTickCount())/1000;
  printf("\n");

  for(i=0;i<n ;i++) printf("c[%2d] = %20.12E\n",i+1,yy[i]);

  i=NumPoints*sizeof(float)*4;
  printf("Allocate %d bytes memory for %d points\n",i,NumPoints);
  pm=malloc(i);
  if( pm == NULL ) {
    puts( "Unable to allocate memory\n" );
    return -1;
  }

  Delta=DeltaSolve;  mj=NumPoints;  mmj=(int)(t+mj*Delta);
  printf("Calculate until %d \n",mmj);
  second[2]=(double)(GetTickCount())/1000;
  for (j=0;j<mj;j++) {
        p(j,0)= (float)yy[2];                // Br-  2
        p(j,1)= (float)yy[5];                // Me   5
        p(j,2)= (float)yy[7];                // Br2  7
        p(j,3)= (float)yy[9];                // Br'  9
//printf("%e %e %e %e\n",t,yy[2],yy[5],yy[7],yy[9]);
    if (j%1000==0) {printf("\r %6.2f %%",t*100/mmj); fflush(stdout);}
    stepx();
  }
  second[3]=(double)(GetTickCount())/1000;
  printf("\nCalculation complete, elapsed time: %10.3f and %10.3f seconds\n",second[1]-second[0],second[3]-second[2]);
  stp=ms[6];  fun=ms[7];  jac=ms[8];  lum=ms[9];  slt=ms[10];
  printf("Funs: %d Jacs: %d LUm: %d Slt: %d Steps: %d\n",fun,jac,lum,slt,stp);

Show:
   pixcount=j-1;
   pxc=pixcount;

   puts("Search Max&Min");
   maxpx=minpx=p(0,0);
   maxpy=minpy=p(0,1);
   maxpz=minpz=p(0,2);
   maxpw=minpw=p(0,3);
   for(c=1;c<pxc;c++) {
           sx=p(c,0); if(sx>maxpx) maxpx=sx; else if(sx<minpx) minpx=sx;
           sx=p(c,1); if(sx>maxpy) maxpy=sx; else if(sx<minpy) minpy=sx;
           sx=p(c,2); if(sx>maxpz) maxpz=sx; else if(sx<minpz) minpz=sx;
           sx=p(c,3); if(sx>maxpw) maxpw=sx; else if(sx<minpw) minpw=sx;
   }
   printf("Min[2]=%E Max[2]=%E\n",minpx,maxpx);
   printf("Min[5]=%E Max[5]=%E\n",minpy,maxpy);
   printf("Min[7]=%E Max[7]=%E\n",minpz,maxpz);
   printf("Min[9]=%E Max[9]=%E\n",minpw,maxpw);
   
   if(fabs(maxpx-minpx)<1e-10 && fabs(maxpy-minpy)<1e-10 && fabs(maxpz-minpz)<1e-10) {
     puts("I am sorry, but you have the attracting point only.");
     return 1;
   }

   puts("Stretching");
   for(c=0;c<pxc;c++) {
     p(c,0)=(p(c,0)-minpx)/(maxpx-minpx)-0.5;
     p(c,1)=(p(c,1)-minpy)/(maxpy-minpy)-0.5;
     p(c,2)=(p(c,2)-minpz)/(maxpz-minpz)-0.5;
     p(c,3)=(p(c,3)-minpw)/(maxpw-minpw);
   }

   maxcor=1;

   puts("Go to graph");
   help();
   printf("Flow parameter = %14.8E\n", xflow);
   fflush(stdout);
   
   auxInitDisplayMode (AUX_DOUBLE | AUX_RGB | AUX_ACCUM | AUX_DEPTH24);
   auxInitPosition (0, 0, 700, 700);
   auxInitWindow ("BZPhase");

   if (myinit()!=0) {
     puts("Error OpenGL initialization.");
     return 2;
   };

  auxReshapeFunc (myReshape);
  auxKeyFunc (AUX_UP,      rotx1);
  auxKeyFunc (AUX_DOWN,    rotx2);
  auxKeyFunc (AUX_LEFT,    roty1);
  auxKeyFunc (AUX_RIGHT,   roty2);
  auxKeyFunc (AUX_SPACE,   move0);
  auxKeyFunc (AUX_x,       rotz1);
  auxKeyFunc (AUX_z,       rotz2);
  auxKeyFunc (AUX_1,        dis1);
  auxKeyFunc (AUX_2,        dis2);
  auxKeyFunc (AUX_3,        dis3);
  auxKeyFunc (AUX_4,    glmodesw);
  auxKeyFunc (AUX_v,      psize1);
  auxKeyFunc (AUX_c,      psize2);

  auxKeyFunc (AUX_d,      movex1);
  auxKeyFunc (AUX_a,      movex2);
  auxKeyFunc (AUX_w,      movey1);
  auxKeyFunc (AUX_s,      movey2);
  auxKeyFunc (AUX_q,      movez1);
  auxKeyFunc (AUX_e,      movez2);

  auxKeyFunc (AUX_m,      scale);
  auxKeyFunc (AUX_p,      print);

  auxKeyFunc (AUX_r,      par1);
  auxKeyFunc (AUX_t,      par2);
  auxKeyFunc (AUX_y,      par3);
  auxKeyFunc (AUX_u,      par4);
  auxKeyFunc (AUX_i,      par5);
  auxKeyFunc (AUX_o,      par6);

  auxKeyFunc (AUX_f,      par7);
  auxKeyFunc (AUX_g,      par8);
  auxKeyFunc (AUX_h,      par9);
  auxKeyFunc (AUX_j,      para);
  auxKeyFunc (AUX_k,      parb);
  auxKeyFunc (AUX_l,      parc);

  auxMouseFunc (AUX_LEFTBUTTON , AUX_MOUSEUP,   Mouse_leftup);
  auxMouseFunc (AUX_LEFTBUTTON , AUX_MOUSELOC,  Mouse_move);
  auxIdleFunc  (Calc);

  auxMainLoop(display);

  return(0);
}
