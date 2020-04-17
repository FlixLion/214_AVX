    #ifndef __WIN32__
    #ifdef WIN32
    #define __WIN32__
    #endif
    #ifdef _WIN32
    #define __WIN32__
    #endif
    #ifdef __WIN32
    #define __WIN32__
    #endif
    #ifdef _WIN64
    #define __WIN32__
    #endif
    #ifdef WIN64
    #define __WIN32__
    #endif
    #ifdef _WINDOWS
    #define __WIN32__
    #endif
#endif

#ifdef __WIN32__
#include <windows.h>     
#include <float.h>
#else
          #define _POSIX_C_SOURCE = 199309L
            #define _POSIX_TIMERS 1  
#endif
        
#include "mex.h"
#include <math.h>
#include <stdio.h>

#ifndef _PORTABLEFPU_H_INCLUDED
#include "portableFPU.h"
#endif

#ifndef _PSTDINT_H_INCLUDED
#include "pstdint_new.h"
#endif


//extern __declspec(dllimport) int _imp_pthread_join();
//extern __declspec(dllimport) int _imp_pthread_create();
/*#ifdef __cplusplus
extern "C" {
#endif
 __declspec(dllimport) int __imp_pthread_join();

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllimport) __imp_pthread_create();

#ifdef __cplusplus
}
#endif*/


//enable/disable debugging output
#define addsig2vol_debug
//#undef addsig2vol_debug
        
// enable/disable pthreads
#define p_threads

//enable/disable C-CODE version (disabled is Asm-code)
#undef C_CODE

#ifdef p_threads
#include "pthread.h"
#define NUMCORES 32 //potential maximum value for init structs etc.
#else
#define NUMCORES 1
#endif

//define interpol-ratio of AScan (for asm fix)
#define INTERP_RATIO 5    // resizes ascan from 3000x1 -> 15000x1 (FIXED SIZE lin. interp.)

//define min-Voxel Number for X (parallel voxels in X_pipe; for 64bit code = 4voxel, for 32 bit = 2voxel)
#define MIN_VOXEL 4

//define matlab in and out 
#define out       plhs[0]
#define out2      plhs[1]

#define AScan     prhs[0]
#define pix_vect  prhs[1]
#define rec_pos   prhs[2]
#define send_pos  prhs[3]
#define speed     prhs[4]
#define res       prhs[5]
#define timeint   prhs[6]
#define IMAGE_XYZ prhs[7]
#define IMAGE_SUM prhs[8]

//global variable 
unsigned int addsig2vol_mode=0; //mode = single ascan, constant speed, 1=blocked, constant speed, 2= unblocked with soundmap
 

//for MT handling
static int64_t latency[NUMCORES]={0}; //uint64_t nSec
static int64_t throughput[NUMCORES]={0}; //uint64 MVoxel/s
static uint32_t nCores_bench = -1;
static uint32_t nCores = NUMCORES; //used value might be reduced by imagesize, benchmark etc

//TEST-CASE 0
// tic; x=100; image_1=addsig2vol_3(rand([3000 1]),single(ones(3,1)),single(ones(3,1)),single(ones(3,1)),single(ones(1,1)),single(ones(1,1)),single(ones(3,1)),uint32([x,x,x]),zeros([x,x,x]),2); toc, j=reshape(image_1,[x x x]); imagesc(j(:,:,1));
// inter_p test
// ascan=[repmat(1000,[1 100]) repmat(0,[1 2900])]';
// tic; x=100; image_1=addsig2vol_3(ascan,single(ones(3,1)),single(ones(3,1)),single(ones(3,1)),single(ones(1,1)),single(ones(1,1)),single(ones(3,1)),uint32([x,x,x]),zeros([x,x,x]); toc, j=reshape(image_1,[x x x]); imagesc(j(:,:,1));

//TEST-CASE 1
//x=128; addsig2vol_3(4), image_1=zeros([x,x,x]); rand('seed',0); x=128; for i=1:2 ascan=rand([1 3000]); [image_1,kkk]=addsig2vol_3(ascan,rand([3 1],'single'),10.*rand([3 1],'single'),400.*rand([3 1],'single'),rand([1 1],'single'),rand([1 1],'single'),rand([1 1],'single'),uint32([x,x,x]),image_1); end, j=reshape(image_1,[x x x]); imagesc(j(:,:,1));
//j=reshape(image_1-image_2,[x x x]); figure; imagesc(j(:,:,128));
 
//// TEST--CASE Blocked
// rng('default');
// count=2; senderPos = 0.01.*rand(3,count); receiverPos = 0.01.*rand(3,count); IMAGE_STARTPOINT = [0,0,0]; IMAGE_RESOLUTION= 0.001; Speed=1500; TimeInterval=1e-7; DataLength=3000; Data=zeros(3000,count); Data(floor(DataLength.*rand(count,1)),1:count)=1;
// x=100; bild=addsig2vol_3(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));

// 4x times same parameter to be have compatible win64&linux64 calling convention
 typedef   struct   /*struct reordered to circumvent alignment problems, 64pointer then 32bit values*/
        {
        double*outz;
        double*AScanz;
        double*out_complexz;
        double*bufferz;
        float*pix_vectz;  
		double*buffer_complexz;
        float*rec_posz; 
        float*send_posz;
        float*speedz;
        float*resz;
        float*timeintz;
        double*AScan_complexz;
		    double *IMAGE_SUMz;
        double *IMAGE_SUM_complexz;
        unsigned int n_Yz;
        unsigned int n_Zz;
        unsigned int n_AScanz;
        unsigned int n_Xz;
        unsigned int qwb0;
        unsigned int qwb1;
        unsigned int qwb2;
        unsigned int qwb3;
		    } Addsig2vol_param;
		    
  /* typedef   struct   //struct reorded to circumvent alignment problems, 64pointer then 32bit values
        {
        double*outz;
        double*out_complexz;
        double*AScanz;
        double*AScan_complexz;
        double*bufferz;
        double*buffer_complexz;
        double *IMAGE_SUMz;
        double *IMAGE_SUM_complexz;
        float*pix_vectz;  
        float*rec_posz; 
        float*send_posz;
        float*speedz;
        float*resz;
        float*timeintz;
        unsigned int n_AScanz;
		    unsigned int n_Xz;
        unsigned int n_Yz;
        unsigned int n_Zz;
        } Addsig2vol_param;*/

//CPUcount 		    
uint64_t CPUCount(void);
uint64_t TimeCounter(void);

//fpu
void fpu_check(void);
         
//ellipsoide backprojections		    
void as2v_complex(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*); 
void as2v_complex_sm(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*); 
void as2v_c(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);
void as2v_MT(double*outz, double*AScanz, unsigned int n_AScanz, double*bufferz, float*pix_vectz,
		    unsigned int n_Xz, float*rec_posz, float*send_posz, float*speedz, float*resz,
		    float*timeintz, double*AScan_complexz,
		    double*buffer_complexz, double*out_complexz, unsigned int n_Yz, unsigned int n_Zz,
		    double* IMAGE_SUMz, double* IMAGE_SUM_complexz);

//Xsum and interpol
void xsum_complex(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*); 
void xsum_c(Addsig2vol_param *, Addsig2vol_param*, Addsig2vol_param*, Addsig2vol_param*);

//thread function		    
void *thread_function(void *arg);

//thread benchmark
void as2v_bench( uint64_t lat[],  uint64_t through[]);

///////////////////End declarations



void mexFunction (int nlhs, mxArray*plhs[],
		  int nrhs, const mxArray*prhs[]) {

 
	//int width_test;
  const mwSize* mwSizePtr;
	unsigned int i;
	uint32_t n_AScan;
	uint32_t n_AScan_block;
	double *AScan_pi;
	  
	float *rec_vec_ptr=NULL;
    float *send_vec_ptr=NULL;    
    float *speed_vec_ptr=NULL;
	unsigned int n_rec_vec_block;
	unsigned int n_send_vec_block;
  	unsigned int n_speed_vec_block;
	unsigned int n_X=1;
	unsigned int n_Y=1;
	unsigned int n_Z=1;
	unsigned int n_IMAGE=0;
	//int*IMAGE_ptr;
	mwSize*BUFFER_ptr;
	mwSize setImageDim[3]={0};
	//mwSize setBufferDim[2];
    mwSize setBufferDim[2]={0};

	double*pr;
	double*pi;

	mxArray *buffer;

//     mxArray* ascan_output_buffer;
//     mxArray* recpos_output_buffer;
//     mxArray* recpos_org = rec_pos;
	if (nlhs > 2) mexErrMsgTxt("Too many output arguments.");

	switch (nrhs) {
	default:
		mexErrMsgTxt("Incorrect number of arguments.\n");
		break;
	case 0:
		mexPrintf("\naddSig2Vol_2 SSE1 Assembler Optimized 64bit LINUX&Windows v3.1 (Multiple Rec-AScan Vers.)\n\n Calculate the ellip. backprojection.\nUses SSE.\n\n#define out       plhs[0] (Double(1:end))\n#define out2      plhs[1] DEBUG (Double(1:end))\n#define AScan     prhs[0] (Double(NxM))\n#define pix_vect  prhs[1] (Single(1:3))\n#define rec_pos   prhs[2] (Single(1:3xM) or Single(1:3x1))\n#define send_pos  prhs[3] (Single(1:3xM) or Single(1:3x1))\n#define speed     prhs[4] (Single (1x1 or 1xM))\n#define res       prhs[5] (Single)\n#define timeint   prhs[6] (Single)\n#define IMAGE_XYZ prhs[7] (UINT32(1:3))\n#define IMAGE_SUM prhs[8] (Double(1:end))\nFeatures: Win&Linux, 32&64bit version, SSE1-Assembler&C-Implementation, Multithreaded by PosixThreads (under windows pthreadsVC2.dll might be needed)\n\t %s. M.Zapf KIT-IPE\n\n",__DATE__);
		
		//benchmark
		as2v_bench(&throughput[0], &latency[0]);
		break;
        
    case 1:
        if (((int)*((double*)mxGetPr(prhs[0]))<=NUMCORES) & ((int)*((double*)mxGetPr(prhs[0]))>=1))
        {
            mexPrintf("Clear Benchmark results, manual set used Cores to %i\n", (int)*((double*)mxGetPr(prhs[0])));
            nCores_bench = (uint32_t) ceil(*((double*)mxGetPr(prhs[0])));
            
            //fill 
            for (i=NUMCORES;i>0;i--)
            { throughput[i-1]=1;
              latency[i-1]=10000000;
            }
            throughput[nCores_bench-1]=299;
            latency[nCores_bench-1]=1;
        }
		break;
	case 9:
    
    //init variables
		//IMAGE_ptr       = *(uint32_t*) mxGetPr(IMAGE_XYZ);
		n_X              = *((int*) mxGetPr(IMAGE_XYZ));//*IMAGE_ptr;
		n_Y              = *((int*) mxGetPr(IMAGE_XYZ)+1);//*(IMAGE_ptr + 1);
		 
		//% workaround for matlab behaviour for size(IMAGE(1 x1x1))->reduced to 1 x1 (some kind of squeeze)
		if (mxGetNumberOfElements(IMAGE_XYZ) < 3)
			n_Z = 1;
		else
			n_Z  = *(((int*) mxGetPr(IMAGE_XYZ)+2));//*(IMAGE_ptr + 2);
            
    	//check if X_Dim >= MIN_VOXEL
		if (n_X < MIN_VOXEL)
		{
		out     = mxCreateDoubleMatrix(0, 0, mxCOMPLEX); //out  == empty
	  mexPrintf("Error: X-Dim has to be at least %d (X is %u, Y is %u, Z %u). Use Y-dim for such small sizes.\n",MIN_VOXEL,n_X,n_Y,n_Z);
	  return;
	  }
    //number of voxel
		n_IMAGE   = n_X * n_Y * n_Z;
				
		//check code if parameter are coherent, if not exit
	  if (mxGetNumberOfElements(IMAGE_SUM)!=n_IMAGE)
	  {
	  out     = mxCreateDoubleMatrix(0, 0, mxCOMPLEX); //out  == empty
	  mexPrintf("Error: Size mismatch of [X Y Z] and size of image, break\n");
	  return;
	  }
	    #ifdef addsig2vol_debug  
		mexPrintf("mxGetNumberOfElements(IMAGE_SUM):  %i\n", mxGetNumberOfElements(IMAGE_SUM));
		mexPrintf("n_IMAGE:  %i\n", nCores_bench);
		#endif
                
	 #ifdef addsig2vol_debug 
    //fpu_status
     fpu_check();
	 #endif
                
	 ////benchmark of performance, selecting number of threads
    if (nCores_bench==-1) as2v_bench(&throughput[0], &latency[0]); 
    
    //select if use or not use multithreading
     #ifdef addsig2vol_debug 
      mexPrintf("selectedCore perf: %f , %f\n",((double)n_IMAGE/((double)throughput[nCores_bench-1] * 1000000))+((double)latency[nCores_bench-1]/1000000000),(((double)n_IMAGE/((double)throughput[0]*1000000))+((double)latency[0]/1000000000)));
      mexPrintf("nimage %i, throughput: %f , latency %f, %f, %f\n",n_IMAGE, (double)throughput[nCores_bench-1],(double)latency[nCores_bench-1],(double)throughput[0],(double)latency[0]);
     #endif     
     
    if ( ( ((double)n_IMAGE/((double)throughput[nCores_bench-1] *1000000))+((double)latency[nCores_bench-1]/1000000000)) <= (((double)n_IMAGE/((double)throughput[0] *1000000))+((double)latency[0]/1000000000)))
         {  nCores = nCores_bench;}
    else {
      #ifdef addsig2vol_debug 
      mexPrintf("Overhead to big, switch to single thread.\n");
      #endif
      nCores = 1;}
      
    #ifdef addsig2vol_debug  
		mexPrintf("selectedNumCores:  %i\n", nCores);
		mexPrintf("savedNumCORE:  %i\n", nCores_bench);
		mexPrintf("perf_MT:  %e\n", ( (throughput[nCores_bench] * n_IMAGE)+latency[nCores_bench]));
		mexPrintf("perf_single:  %e\n", ((throughput[1]*n_IMAGE)+latency[1]));
	#endif
	   
	  
	  //////schlauchimage code (old)
		//setImageDim[0]  = n_IMAGE;
		//setImageDim[1]  = 1;        //z.b: 400000x1
					
		//////NOT Schlauchimage Code
		setImageDim[0]  = n_X;
		setImageDim[1]  = n_Y;        
		setImageDim[2]  = n_Z;    

    ////check for BLOCKED/UNBlocked Parameters
		BUFFER_ptr      = (mwSize*) mxGetDimensions(AScan);
		n_AScan         = (*BUFFER_ptr);        //gesamtanzahl elemente IN EINEM ASCAN!!!
		n_AScan_block   = *(BUFFER_ptr + 1);    //2 dim %number of parallel

    mwSizePtr       = mxGetDimensions(rec_pos); //dimension of Rec-vec
    n_rec_vec_block = (unsigned int) mwSizePtr[1];//*(int*)(mwPtr + 1); second DIM
    mwSizePtr       = mxGetDimensions(send_pos); //dimension of send-vec
    n_send_vec_block= (unsigned int) mwSizePtr[1];//*(int*)(send_vec_ptr + 1);
    mwSizePtr       = mxGetDimensions(speed); //dimension of send-vec
    n_speed_vec_block= (unsigned int) mwSizePtr[1];//*(int*)(send_vec_ptr + 1);
    
    //soundmap version ?
//mexPrintf("Info: Soundmap versionrrr\n");
    if ( ((unsigned int) mwSizePtr[0]==n_X) & ((unsigned int) mwSizePtr[1]==n_Y) & ((unsigned int) mwSizePtr[2]==n_Z | mxGetNumberOfDimensions(speed)==2) )
    { mexPrintf("Info: Soundmap version\n");
    addsig2vol_mode = 2;
    }    
    else{ //not soundmapversion; blocked version? 
        
        if ( ((*(unsigned int*)mxGetDimensions(rec_pos))!=3) | ((*(unsigned int*)mxGetDimensions(send_pos))!=3)) //check first dimension
        {
        #ifdef addsig2vol_debug
        mexPrintf("Dim1 send_vec: %d, Dim1 rec_vec: %d",*(int*)rec_vec_ptr,*(int*)send_vec_ptr);
        mexPrintf("ascan_block: %d, blocksize rec_vec: %d, blocksize send_vec: %d",n_AScan_block,n_rec_vec_block,n_send_vec_block);
         #endif
        out     = mxCreateDoubleMatrix(0, 0, mxCOMPLEX); //out  == empty
        mexPrintf("Error: 3-d vectors needed for emitter & receiver positions or transposed blocked pos (1x3 instead of 3x1), break\n");
        return;}
    

         #ifdef addsig2vol_debug
            mexPrintf("nascan: %d,ascan_block: %d, blocksize rec_vec: %d, blocksize send_vec: %d",n_AScan,n_AScan_block,n_rec_vec_block,n_send_vec_block);
         #endif

         if (!( (((n_AScan_block == n_rec_vec_block) & (n_AScan_block == n_send_vec_block)) | ((1 == n_rec_vec_block) & (n_AScan_block == n_send_vec_block)) | ((n_AScan_block == n_rec_vec_block) & (1 == n_send_vec_block))) & ((n_AScan_block == n_speed_vec_block) | (1 == n_speed_vec_block)) ))
          {out     = mxCreateDoubleMatrix(0, 0, mxCOMPLEX); //out  == empty
          mexPrintf("Error: Blocked sizes parameter mismatch. Size(AScan,2) has to be size(rec_vec,2) and/or size(send_vec,2), n_rec_vec_block:%i n_send_vec_block:%i n_speed_vec_block:%i break\n",n_rec_vec_block, n_send_vec_block, n_speed_vec_block);
          return;}


         if (!(((mxGetPi(AScan)==NULL) & (mxGetPi(IMAGE_SUM)==NULL)) | ((mxGetPi(AScan)!=NULL) & (mxGetPi(IMAGE_SUM)!=NULL))))
         {mexPrintf("Error: Mismatch complex/real AScan vs complex/real sum image, break\n");
          return;

          //if here blocked mode successful initalized!
          addsig2vol_mode = 1;
        }
         else
             //normal version
              addsig2vol_mode = 0;
    }
		setBufferDim[0] = (*BUFFER_ptr) * INTERP_RATIO;
		setBufferDim[1] = 1;   //z.b: 400000x1
        
        #ifdef addsig2vol_debug  
        mexPrintf("(*BUFFER_ptr) :  %i\n", (*BUFFER_ptr));
		mexPrintf("mxGetNumberOfElements(IMAGE_XYZ):  %i\n",mxGetNumberOfElements(IMAGE_XYZ));
		#endif
        
		//DEBUG
		//if (n_AScan<(2*wid+2)) mexErrMsgTxt("1. Array not 2 times +2 greater value 3 ");
		/*mexPrintf("n_AScan :  %i\n\n", n_AScan );
		mexPrintf("n_AScan_block :  %i\n\n", n_AScan_block );
		mexPrintf("n_Z:  %i\n\n", n_Z);
		mexPrintf("n_Y :  %i\n\n", n_Y);
		mexPrintf("n_X :  %i\n\n", n_X);
		mexPrintf("n_IMAGE :  %i\n\n", n_IMAGE);
		mexPrintf("ewcpos:  %f\n\n", (float*)mxGetPr(rec_pos));
		mexPrintf("sendpos:  %f\n\n", (float*)mxGetPr(send_pos));
		mexPrintf("pix_vectz :  %f\n\n",(float*)mxGetPr(pix_vect));
		mexCallMATLAB(0, NULL, 1, &send_pos, "disp");*/
//         ascan_output_buffer = mxCreateDoubleMatrix(3000,1,mxREAL);
//         recpos_output_buffer = mxCreateNumericMatrix(3,1,mxSINGLE_CLASS,mxREAL);

		//Check if complex
		AScan_pi = mxGetPi(AScan);

		//anlegen puffer fuer xsum & image
		if (AScan_pi != NULL) 
		    { //complex
			out     = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);//out        = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);             //out     = mxCreateDoubleMatrix(n_Index,1,mxCOMPLEX);
			if (out==NULL) {mexPrintf("Image_SUM alloc failed!\n");return;} 
            mxSetDimensions(out, setImageDim, mxGetNumberOfElements(IMAGE_XYZ)); // UNSAFE if SCHLAUCHINPUT!!! mxGetNumberOfDimensions(IMAGE_SUM)); //bsp. 3000x1  -> (3000,1) ,2                      //bsp. 3000x1  -> (3000,1) ,2
			pr = mxMalloc(n_IMAGE * sizeof(double));
            if (pr==NULL) {mexPrintf("Image_SUM alloc failed!\n");return;} 
			pi = mxMalloc(n_IMAGE * sizeof(double));
            if (pi==NULL) {mexPrintf("Image_SUM alloc failed!\n");return;}
			mxSetPr(out, pr);
			mxSetPi(out, pi);
			
			buffer = mxCreateDoubleMatrix(0, 0, mxCOMPLEX);                         //Sum buffer laenge ascan
			mxSetDimensions(buffer, setBufferDim, mxGetNumberOfDimensions(AScan));  //bsp. 3000x1  -> (3000,1) ,2
			pr = mxMalloc(INTERP_RATIO * n_AScan * sizeof(double));
			pi = mxMalloc(INTERP_RATIO * n_AScan * sizeof(double));
			mxSetPr(buffer, pr);
			mxSetPi(buffer, pi);
		} else { //real
			out        = mxCreateDoubleMatrix(0, 0, mxREAL);        //out     = mxCreateDoubleMatrix(n_Index,1,mxCOMPLEX);
			if (out==NULL) {mexPrintf("Image_SUM alloc failed!\n");return;} 
            mxSetDimensions(out, setImageDim,mxGetNumberOfElements(IMAGE_XYZ)); // UNSAFE if SCHLAUCHINPUT!!!mxGetNumberOfDimensions(IMAGE_SUM)); //bsp. 3000x1  -> (3000,1) ,2
			pr = mxMalloc(n_IMAGE * sizeof(double));
            if (pr==NULL) {mexPrintf("Image_SUM alloc failed!\n");return;}
			mxSetPr(out, pr);
			mxSetPi(out, NULL);

			buffer = mxCreateDoubleMatrix(0, 0, mxREAL);              //Sum buffer laenge ascan
            if (buffer==NULL) {mexPrintf("buffer alloc failed!\n");return;}
			mxSetDimensions(buffer, setBufferDim, mxGetNumberOfDimensions(AScan));  //bsp. 3000x1  -> (3000,1) ,2
			pr = mxMalloc(INTERP_RATIO * n_AScan * sizeof(double));
            if (pr==NULL) {mexPrintf("buffer alloc failed!\n");return;}
			mxSetPr(buffer, pr);
			mxSetPi(buffer, NULL);

		}
	    #ifdef addsig2vol_debug  
        mexPrintf("mxGetNumberOfElements(IMAGE_SUM):  %i\n", mxGetNumberOfElements(out));
		mexPrintf("mxGetNumberOfElements(IMAGE_XYZ):  %i\n",mxGetNumberOfElements(IMAGE_XYZ));
		mexPrintf(" mxGetNumberOfElements(out):  %i\n", mxGetNumberOfElements(out));
		mexPrintf(" mxGetNumberOfElements(buffer):  %i\n", mxGetNumberOfElements(buffer));
        mexPrintf("  setImageDim[0]:  %i\n", setImageDim[0]);
        mexPrintf("  setImageDim[1]:  %i\n", setImageDim[1]);
        mexPrintf("  setImageDim[2]:  %i\n", setImageDim[2]);
		#endif
               
		////first Ascan
		// combined REAL & COMPLEX VERSION
    
		//no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
		as2v_MT(mxGetPr(out), mxGetPr(AScan), n_AScan, mxGetPr(buffer),
			     (float*)mxGetPr(pix_vect), n_X, (float*)mxGetPr(rec_pos), (float*)mxGetPr(send_pos),
			     (float*)mxGetPr(speed),
			     (float*)mxGetPr(res), (float*)mxGetPr(timeint), mxGetPi(AScan), mxGetPi(
				     buffer), mxGetPi(
				     out), n_Y, n_Z, mxGetPr(IMAGE_SUM), mxGetPi(
				     IMAGE_SUM));
		
		//loop over ascans > 1
		for (i = 2; i <= n_AScan_block; i++) {
			//check for complex ascan only increase if available because NULL-Pointer +something -> not anymore a nullpointer!
			if (AScan_pi != NULL)	AScan_pi = mxGetPi(AScan) + (n_AScan * (i - 1)); //set to next value
            if (1<n_rec_vec_block)  rec_vec_ptr  = (float*)mxGetPr(rec_pos)  + (3 * sizeof(float) * (i - 1)); else  rec_vec_ptr  = (float*)mxGetPr(rec_pos);
            if (1<n_send_vec_block) send_vec_ptr = (float*)mxGetPr(send_pos) + (3 * sizeof(float) * (i - 1)); else send_vec_ptr  = (float*)mxGetPr(send_pos);
            if (1<n_speed_vec_block) speed_vec_ptr = (float*)mxGetPr(speed) + (1 * sizeof(float) * (i - 1)); else speed_vec_ptr  = (float*)mxGetPr(speed);

// 			mexPrintf("delta_address :  %i\n\n", n_AScan*(i-1));
// 			mexPrintf("n_AScan_address :  %p\n\n", mxGetPr(AScan)+(n_AScan)*(i-1));
// 			mexPrintf("mxGetPr(out):  %p\n\n", mxGetPr(out));
// 			mexPrintf("mxGetPi(out):  %p\n\n", mxGetPi(out));
// 			mexPrintf("mxGetPr(Ascan):  %p\n\n", mxGetPr(AScan));
// 			mexPrintf("mxGetPi(Ascan):  %p\n\n", mxGetPi(AScan));
// 			mexPrintf("i:  %i\n\n", i);

			// combined REAL & COMPLEX VERSION
			//no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
			//matlab seems to do some nasty errors on ptr to float...adding 8bytes instead of 4!!!! workaround implemented
		as2v_MT(mxGetPr(out), mxGetPr(
					     AScan) + (n_AScan * (i - 1)), n_AScan, mxGetPr(
					     buffer), (float*)mxGetPr(
					     pix_vect), n_X, (float*)rec_vec_ptr,
				     (float*)send_vec_ptr, (float*)speed_vec_ptr, (float*)mxGetPr(res), (float*)mxGetPr(
					     timeint), AScan_pi, mxGetPi(buffer),
				     mxGetPi(out), n_Y, n_Z, mxGetPr(out), mxGetPi(
					     out));
		
		
		/*	as2v_MT(mxGetPr(out), mxGetPr(
					     AScan) + (n_AScan * (i - 1)), n_AScan, mxGetPr(
					     buffer), mxGetPr(
					     pix_vect), n_X, (char*)mxGetPr(
					     rec_pos) + (3 * sizeof(float) * (i - 1)),
				     mxGetPr(send_pos), mxGetPr(
					     speed), mxGetPr(res), mxGetPr(
					     timeint), AScan_pi, mxGetPi(buffer),
				     mxGetPi(out), n_Y, n_Z, mxGetPr(out), mxGetPi(
					     out));*/
			//  test
			//as2v_complex(mxGetPr(out),mxGetPr(AScan)+(n_AScan*(i-1)),n_AScan,mxGetPr(buffer),mxGetPr(pix_vect),n_X,(char*)mxGetPr(rec_pos)+(3*sizeof(float)*(i-1)),mxGetPr(send_pos),mxGetPr(speed),mxGetPr(res),mxGetPr(timeint),AScan_pi,NULL,NULL,n_Y,n_Z,mxGetPr(out),NULL,qwordbuffer);

//             mxSetPr(ascan_output_buffer,mxGetPr(AScan)+(n_AScan*(i-1)));
//             mxSetPr(recpos_output_buffer,(char*)mxGetPr(recpos_org)+(3*sizeof(float)*(i-1)));
//
//             mexCallMATLAB(0, NULL, 0, NULL, "figure");
//             mexCallMATLAB(0, NULL, 1, &ascan_output_buffer, "plot");
//             mexCallMATLAB(0, NULL, 1, &recpos_output_buffer, "disp");
		
            }
		   out2 = buffer;
       //mxDestroyArray(buffer);
		   break;	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	} //end switch

	return;
}



//////////////////////////////////////////////////////////////////////////////////

void as2v_MT(double*outz, double*AScanz, unsigned int n_AScanz, double*bufferz, float*pix_vectz,
		    unsigned int n_Xz, float*rec_posz, float*send_posz, float*speedz, float*resz,
		    float*timeintz, double*AScan_complexz,
		    double*buffer_complexz, double*out_complexz, unsigned int n_Yz, unsigned int n_Zz,
		    double *IMAGE_SUMz, double *IMAGE_SUM_complexz)

{

  //pthread variables
#ifdef p_threads
    pthread_t mythread[NUMCORES]; //numCPU -1
    int rc = 0; //return-value from thread functions
#endif
    Addsig2vol_param threadArg[NUMCORES]; //numCPU -1 BUT for NUM 1 
    float pix_vecz_buffer[NUMCORES][3]= {0,0,0};
    unsigned int n_Zz_start = 0;
    unsigned int n_Zz_num = 0;
    //unsigned int nCores = NUMCORES;
    unsigned int i = 0;
      
      
    ////////Z-Layer multithreading, use multi-threading if enough Z-layers imaged
    if (n_Zz>1)
      {
      //limit to usefull number
      if (n_Zz < nCores)  nCores = n_Zz;
    
      #ifdef addsig2vol_debug
      mexPrintf("Z-Dim multithreading\n");
      #endif
      
      //Generate parameter structs for Z-multithreading
		  for (i=0;i<=nCores-1;i++) //WITH mainthread! (=)!
        {
        n_Zz_num  = (unsigned int)(floor(n_Zz/nCores));
        n_Zz_start = n_Zz_num*i*n_Xz*n_Yz;
        
         //set picture startpoint       
        pix_vecz_buffer[i][0]= *pix_vectz;
        pix_vecz_buffer[i][1]= *(pix_vectz+1);
        pix_vecz_buffer[i][2]= (*(pix_vectz+2))+n_Zz_num*i* *resz;
        
        if (i==nCores-1)
        {      //compensated number for FLOOR/CEIL rounding probs
        n_Zz_num = (unsigned int)(n_Zz-floor(n_Zz/nCores)*(nCores-1));}
        
        /* mexPrintf("delta_address :  %p\n\n", (pix_vectz));
        mexPrintf("z :  %f\n\n", *(pix_vectz+2));
        mexPrintf("z_num :  %i\n\n", n_Zz_num);
        mexPrintf("z_start :  %i\n\n", n_Zz_start);
        mexPrintf("x :  %f\n\n", pix_vecz_buffer[i][0]);
        mexPrintf("y :  %f\n\n", pix_vecz_buffer[i][1]);
        mexPrintf("z :  %f\n\n", pix_vecz_buffer[i][2]);
        mexPrintf("p_x :  %p\n\n", &(pix_vecz_buffer[i][0]));*/
      
        //fill parameter struct
        threadArg[i].outz=outz+n_Zz_start;//
        threadArg[i].AScanz=AScanz;
        threadArg[i].n_AScanz=n_AScanz;
        threadArg[i].bufferz=bufferz;
        threadArg[i].pix_vectz=&(pix_vecz_buffer[i][0]);/*mxGetPr(pix_vect)*/
        threadArg[i].n_Xz=n_Xz;
        threadArg[i].rec_posz=rec_posz;
        threadArg[i].send_posz=send_posz;
        threadArg[i].speedz=speedz;
        threadArg[i].resz=resz;
        threadArg[i].timeintz=timeintz;
        threadArg[i].AScan_complexz=AScan_complexz;
        threadArg[i].buffer_complexz=buffer_complexz;
        threadArg[i].out_complexz=out_complexz+n_Zz_start;//
        threadArg[i].n_Yz=n_Yz;
        threadArg[i].n_Zz=n_Zz_num;/*n_Zz*/
        threadArg[i].IMAGE_SUMz=IMAGE_SUMz+n_Zz_start;//
        threadArg[i].IMAGE_SUM_complexz=IMAGE_SUM_complexz+n_Zz_start; //
        }
      }
      else  {
      
     //Y-DIM multithreading because Z-layer = 1, use multithreading if enough Y-layers imaged
      if (n_Yz>1)
      {
      //limit to usefull number
      if (n_Yz < nCores)  nCores = n_Yz;
    
      #ifdef addsig2vol_debug
      mexPrintf("Y-Dim multithreading\n");
      #endif
    
      //Generate parameter structs for Y-multithreading
	  	for (i=0;i<=nCores-1;i++) //WITH mainthread! (=)!
        {
        n_Zz_num  = (unsigned int)(floor(n_Yz/nCores));
        n_Zz_start = n_Zz_num*i*n_Xz;
        
         //set picture startpoint       
        pix_vecz_buffer[i][0]= *pix_vectz;
        pix_vecz_buffer[i][1]= (*(pix_vectz+1))+n_Zz_num*i* *resz;
        pix_vecz_buffer[i][2]= *(pix_vectz+2);
        
        if (i==nCores-1)
        {      //compensated number for FLOOR/CEIL rounding probs
        n_Zz_num = (unsigned int)(n_Yz-floor(n_Yz/nCores)*(nCores-1));}
        
        //fill parameter struct
        threadArg[i].outz=outz+n_Zz_start;//
        threadArg[i].AScanz=AScanz;
        threadArg[i].n_AScanz=n_AScanz;
        threadArg[i].bufferz=bufferz;
        threadArg[i].pix_vectz=&(pix_vecz_buffer[i][0]);/*mxGetPr(pix_vect)*/
        threadArg[i].n_Xz=n_Xz;
        threadArg[i].rec_posz=rec_posz;
        threadArg[i].send_posz=send_posz;
        threadArg[i].speedz=speedz;
        threadArg[i].resz=resz;
        threadArg[i].timeintz=timeintz;
        threadArg[i].AScan_complexz=AScan_complexz;
        threadArg[i].buffer_complexz=buffer_complexz;
        threadArg[i].out_complexz=out_complexz+n_Zz_start;//
        threadArg[i].n_Yz=n_Zz_num; /*n_Yz*/
        threadArg[i].n_Zz=n_Zz;
        threadArg[i].IMAGE_SUMz=IMAGE_SUMz+n_Zz_start;//
        threadArg[i].IMAGE_SUM_complexz=IMAGE_SUM_complexz+n_Zz_start; //
        }
      
      }else  
        { //no multi-threading
        nCores = 1;
        
        #ifdef addsig2vol_debug
        mexPrintf("No multithreading\n");
        #endif
                  
        //fill parameter struct
        threadArg[nCores-1].outz=outz;//
        threadArg[nCores-1].AScanz=AScanz;
        threadArg[nCores-1].n_AScanz=n_AScanz;
        threadArg[nCores-1].bufferz=bufferz;
        threadArg[nCores-1].pix_vectz=pix_vectz;
        threadArg[nCores-1].n_Xz=n_Xz;
        threadArg[nCores-1].rec_posz=rec_posz;
        threadArg[nCores-1].send_posz=send_posz;
        threadArg[nCores-1].speedz=speedz;
        threadArg[nCores-1].resz=resz;
        threadArg[nCores-1].timeintz=timeintz;
        threadArg[nCores-1].AScan_complexz=AScan_complexz;
        threadArg[nCores-1].buffer_complexz=buffer_complexz;
        threadArg[nCores-1].out_complexz=out_complexz;//
        threadArg[nCores-1].n_Yz=n_Yz;
        threadArg[nCores-1].n_Zz=n_Zz;
        threadArg[nCores-1].IMAGE_SUMz=IMAGE_SUMz;//
        threadArg[nCores-1].IMAGE_SUM_complexz=IMAGE_SUM_complexz; //
        }
      }    
          
      ////////////////////imaging

      //interpol & X-SUM (in the case of NUMCORE=1 only call)  
      #ifdef C_CODE
      xsum_c(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]); 
      #else
      xsum_complex(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);    
      #endif
     
      ////release threads 
      for (i=0;i<nCores-1;i++)
	  {
      #ifdef p_threads 
      rc = pthread_create( &(mythread[i]), NULL, thread_function, (void*)&(threadArg[i]));
      if (rc) { mexPrintf("ERROR: return code from pthread_create() is %d\n", rc); return;}
      #endif
      }
    	
	//imaging-call of last part (in the case of NUMCORE=1 only call)  
      #ifdef C_CODE
      as2v_c(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      #else
      if (addsig2vol_mode==0) as2v_complex(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      if (addsig2vol_mode==2) as2v_complex_sm(&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1],&threadArg[nCores-1]);
      
      #endif
			//catches threads again	     
			for (i=0;i<nCores-1;i++)
		  {
              #ifdef p_threads
              rc = pthread_join ( mythread[i], NULL );
              if (rc) { mexPrintf("ERROR: return code from pthread_join() is %d\n", rc); return;}
              #endif
        }	     

      //set because because potentially reduced by imagesize
      //if (nCores_bench >0) nCores = nCores_bench;      
      
      
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void as2v_c(Addsig2vol_param* tt, Addsig2vol_param* t1, Addsig2vol_param* t2, Addsig2vol_param* t3)
{
    
      ///TODO: COMPLEX part!!!!

float dist_sv[3] = {0,0,0};
float dist_rv[3] = {0,0,0};
float factor = 0;

unsigned int index, image_index = 0;
unsigned int z,y,x,sampl = 0;

//decompose variables from struct
  	    double*outz = tt->outz;
        double*AScanz = tt->AScanz;
        double*out_complexz=tt->out_complexz;
        double*bufferz = tt->bufferz;
        float*pix_vectz = tt->pix_vectz;  
		    double*buffer_complexz=tt->buffer_complexz;
        float*rec_posz = tt->rec_posz; 
        float*send_posz = tt->send_posz;
        float*speedz = tt->speedz;
        float*resz = tt->resz;
        float*timeintz = tt->timeintz;
        double*AScan_complexz = tt->AScan_complexz;
		    double *IMAGE_SUMz = tt->IMAGE_SUMz;
        double *IMAGE_SUM_complexz = tt->IMAGE_SUM_complexz;
        unsigned int n_Yz = tt->n_Yz;
        unsigned int n_Zz = tt->n_Zz;
        unsigned int n_AScanz = tt->n_AScanz;
        unsigned int n_Xz = tt->n_Xz;
		    
///bildgebung

factor =  INTERP_RATIO / (*speedz * *timeintz);

for (z=1; z<=n_Zz; z++)
	{
		//dist_sv[2] = pow(send_posz[2] - (((float)z* *resz)+pix_vectz[2]) ,2);
        dist_sv[2]  = send_posz[2] - (((float)z* *resz)+pix_vectz[2]);
        dist_sv[2] *= dist_sv[2];
		//dist_rv[2] = pow(rec_posz[2] - (((float)z* *resz)+pix_vectz[2]) ,2);
        dist_rv[2]  = rec_posz[2] - (((float)z* *resz)+pix_vectz[2]);
        dist_rv[2] *= dist_rv[2];

	for (y=1; y<=n_Yz; y++)
		{
		//dist_sv[1] = pow(send_posz[1] - ((y* *resz)+pix_vectz[1]) ,2);
        dist_sv[1]  = send_posz[1] - (((float)y* *resz)+pix_vectz[1]);
        dist_sv[1] *= dist_sv[1];
		//dist_rv[1] =  pow(rec_posz[1] - ((y* *resz)+pix_vectz[1]) ,2);
        dist_rv[1]  = rec_posz[1] - (((float)y* *resz)+pix_vectz[1]);
        dist_rv[1] *= dist_rv[1];

		for (x=1; x<=n_Xz; x++)
			{
				//dist_sv[0] = pow(send_posz[0] - (x* *resz),2);
				//dist_rv[0] =  pow(rec_posz[0] - (x* *resz),2);
				dist_sv[0] = send_posz[0] - ((x* *resz)+pix_vectz[0]);
				dist_rv[0] =  rec_posz[0] - ((x* *resz)+pix_vectz[0]);

				    //dist = (sqrt(dist_sv[1]+ dist_sv[2] + (dist_sv[0]*dist_sv[0])) + sqrt(dist_rv[1]+ dist_rv[2]+ (dist_rv[0]*dist_rv[0])) );
				index = (unsigned int) floor( ( sqrt(dist_sv[1]+ dist_sv[2] + (dist_sv[0]*dist_sv[0])) + sqrt(dist_rv[1]+ dist_rv[2]+ (dist_rv[0]*dist_rv[0])) ) * factor);
				    //mexPrintf("index:  %i\n\n", index);
				    //outz[image_index] = (double)index;
				if ((index >= n_AScanz*INTERP_RATIO) | (index < 0))
					outz[image_index] = IMAGE_SUMz[image_index]; //nix addiert
				else
					outz[image_index] = IMAGE_SUMz[image_index] + bufferz[index];//AScanz[index];
					
				image_index++;
			}
		}
	}

}

///////////////////////////////////////////////

void xsum_c(Addsig2vol_param* tt, Addsig2vol_param* t1, Addsig2vol_param* t2, Addsig2vol_param* t3)
{  ///TODO: COMPLEX part!!!!

float dist_sv[3] = {0,0,0};
float dist_rv[3] = {0,0,0};
float factor = 0;

unsigned int image_index = 0;
unsigned int i,j,sampl = 0;

double i_buffer=0.0;

double* sec_buffer;

//decompose variables from struct
  /* as2v_c(tt->outz, tt->AScanz, tt->n_AScanz, tt->bufferz, tt->pix_vectz,
		    tt->n_Xz, tt->rec_posz, tt->send_posz, tt->speedz, tt->resz,
		    tt->timeintz, tt->AScan_complexz,
		    tt->buffer_complexz, tt->out_complexz, tt->n_Yz, tt->n_Zz,
		    tt->IMAGE_SUMz, tt->IMAGE_SUM_complexz);*/
		    
		double*outz = tt->outz;
        double*AScanz = tt->AScanz;
        double*out_complexz=tt->out_complexz;
        double*bufferz = tt->bufferz;
        float*pix_vectz = tt->pix_vectz;  
		    double*buffer_complexz=tt->buffer_complexz;
        float*rec_posz = tt->rec_posz; 
        float*send_posz = tt->send_posz;
        float*speedz = tt->speedz;
        float*resz = tt->resz;
        float*timeintz = tt->timeintz;
        double*AScan_complexz = tt->AScan_complexz;
		    double *IMAGE_SUMz = tt->IMAGE_SUMz;
        double *IMAGE_SUM_complexz = tt->IMAGE_SUM_complexz;
        unsigned int n_Yz = tt->n_Yz;
        unsigned int n_Zz = tt->n_Zz;
        unsigned int n_AScanz = tt->n_AScanz;
        unsigned int n_Xz = tt->n_Xz;
		 
#ifdef addsig2vol_debug
for (i=0;i<n_AScanz*INTERP_RATIO;i++)
{	bufferz[i] = i; //set marking for NON-set or initalized values
}  
#endif       
        
//C-Version needs second buffer
sec_buffer = mxMalloc(INTERP_RATIO * n_AScanz * sizeof(double));


//interp for first Samples
for (i=0;i<(unsigned int) floor(INTERP_RATIO/2);i++)
{	sec_buffer[i] = AScanz[0] * (unsigned int) ((i+1)/floor(INTERP_RATIO/2));
}

//almost all samples
for (i=0;i<n_AScanz-1;i++)
{    	for (j=0;j<INTERP_RATIO;j++)
	{	sec_buffer[(unsigned int) floor(INTERP_RATIO/2)+(i* INTERP_RATIO)+j] = AScanz[i+1] * j/INTERP_RATIO + AScanz[i] * (INTERP_RATIO-j)/INTERP_RATIO; //
	}
}

//interp for last Samples
for (i=0;i<(unsigned int) floor(INTERP_RATIO/2)+1;i++)
{	sec_buffer[n_AScanz*INTERP_RATIO-(unsigned int) floor(INTERP_RATIO/2)-1+i] = AScanz[n_AScanz-1] * ((floor(INTERP_RATIO/2)+1-i) / (floor(INTERP_RATIO/2)+1));
}
///end interp



/////xsum
sampl = (unsigned int)(ceil((float)1.7*(( *resz / *speedz)/ (*timeintz/INTERP_RATIO)) /2)); //halbe breite

i_buffer = 0;
for (i=0;i<sampl;i++)
{	i_buffer = i_buffer + sec_buffer[i]/(2*sampl);
}

for (i=0;i<sampl;i++)
{	if (i+sampl<n_AScanz*INTERP_RATIO){
    i_buffer = i_buffer +sec_buffer[i+sampl]/(2*sampl);
	bufferz[i] = i_buffer;}
}

for (i=sampl;i<(n_AScanz*INTERP_RATIO)-sampl;i++)  
{ if (i+sampl<n_AScanz*INTERP_RATIO){
    i_buffer = i_buffer + sec_buffer[i+sampl]/(2*sampl) - sec_buffer[i-sampl]/(2*sampl); 
	bufferz[i] = i_buffer;}
}

for (i=n_AScanz*INTERP_RATIO-sampl;i<n_AScanz*INTERP_RATIO;i++)
{	if (i-sampl>=0){
    i_buffer = i_buffer - sec_buffer[i-sampl]/(2*sampl);
	bufferz[i] = i_buffer / (sampl-(n_AScanz*INTERP_RATIO)-i);}
}
/////end xsum

//free sec_buffer (only buffer needed now)
mxFree(sec_buffer);

}

///////////////////////////////////////////////

 void *thread_function(void *arg) 
{   
    Addsig2vol_param* tt = (Addsig2vol_param*) arg;

    //imaging call with four times pointer to struct
    #ifdef C_CODE
    as2v_c(arg,arg,arg,arg); //compatible win64 & linxu64 function-call
    #else
    if (addsig2vol_mode==0) as2v_complex(arg,arg,arg,arg);
    if (addsig2vol_mode==2) as2v_complex_sm(arg,arg,arg,arg);
      
    #endif
   
   //decomposing for old function 
   /* as2v_c(tt->outz, tt->AScanz, tt->n_AScanz, tt->bufferz, tt->pix_vectz,
		    tt->n_Xz, tt->rec_posz, tt->send_posz, tt->speedz, tt->resz,
		    tt->timeintz, tt->AScan_complexz,
		    tt->buffer_complexz, tt->out_complexz, tt->n_Yz, tt->n_Zz,
		    tt->IMAGE_SUMz, tt->IMAGE_SUM_complexz);*/
     	
     return NULL;
}

////////////////////////////////////////////////////////////////////

void as2v_bench(uint64_t throughput[], uint64_t latency[])
{      
    #pragma fenv_access (on)  
    uint32_t i,j,l,k;
	uint32_t n_AScan;
	uint32_t n_AScan_block;
	int n_X;
	int n_Y;
	int n_Z;
	int n_IMAGE;
	mwSize setImageDim[3];
	mwSize setBufferDim[2];
	mwSize setAScanDim[2];
	float pix_vec_bench[3]  = {2,1,77};
	float send_vec_bench[3] = {1,2,2};
	float rec_vec_bench[3]  = {5,6,9};
    float float_bench = 1;

	double* pr;
   /* double*pr1;
    double* pr2;
    double*pr3;
    double*pr4;*/
    
    #define MAX_AVERAGE 128
    uint64_t average_buffer[MAX_AVERAGE]={0};
    uint64_t throughput_sort_buffer[NUMCORES]={0};
    uint64_t counter, counter2, stdabw, mittelw;
    uint64_t bench_ref_time=0; //ns
    
    uint64_t minBenchTime = 100000000; //ns
   
    uint32_t minAverage = 7;
    uint32_t average = 8;
    uint32_t nVoxel_throughput = 256;   
    
    mxArray *AScan_bench;
	mxArray *buffer_bench;
	mxArray *out_bench; 
	mxArray *image_sum_bench;
	mxArray *time_bench;    
    
    //autodetect free memory
   /*  do { nVoxel_throughput=(uint32_t) nVoxel_throughput/2;
       pr = mxMalloc(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput * sizeof(double));
     if (pr!=NULL)  mxFree(pr);
       nVoxel_throughput=(uint32_t) nVoxel_throughput/2; //nocnmal halbieren f�r 2Bilder
   } while (pr==NULL);*/
    
  mexPrintf("Benchmarking System:\nLatency with %d Byte Image\nThroughput with %d kByte Image.\n",(int)(MIN_VOXEL*NUMCORES*sizeof(double)),(int)(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput*sizeof(double))/1024);    
  mexPrintf("Threads | Throughput in MVoxel/s                     | Latency in nsec                              | Median-Size\n");    
  mexPrintf("        | Median | Mean   | Min    | Max    | Std    | Median    | Min       | Max       | Std      | \n");    
  mexPrintf("-----------------------------------------------------------------------------------------------------------\n");    
 
  //check for minimum amount
  if (nVoxel_throughput<NUMCORES*MIN_VOXEL) nVoxel_throughput=NUMCORES*MIN_VOXEL;
  /////kill eventually running tic-toc!!! (not needed)
   //mexCallMATLAB(0, NULL, 0, NULL, "toc");
	
	////fix variables
	time_bench = mxCreateDoubleMatrix(1,1,mxREAL);
	n_AScan        = 3000;        //gesamtanzahl elemente IN EINEM ASCAN!!!
	n_AScan_block  = 1;    //2 dim %number of parallel
  
    //buffer
    setBufferDim[0] = n_AScan * INTERP_RATIO;
	setBufferDim[1] = 1;   //z.b: 400000x1
	
	//ascan
	setAScanDim[0] = n_AScan;
	setAScanDim[1] = 1;   //z.b: 400000x1
		
	buffer_bench = mxCreateDoubleMatrix(0, 0, mxREAL);   //Sum buffer laenge ascan
	mxSetDimensions(buffer_bench, setBufferDim, 2);  //bsp. 3000x1  -> (3000,1) ,2
	pr = mxMalloc(INTERP_RATIO * n_AScan * sizeof(double));
	mxSetPr(buffer_bench, pr);
	mxSetPi(buffer_bench, NULL);
    //pr1=pr;
			
	AScan_bench = mxCreateDoubleMatrix(0, 0, mxREAL);   //Sum buffer laenge ascan
	mxSetDimensions(AScan_bench, setAScanDim, 2);  //bsp. 3000x1  -> (3000,1) ,2
	pr = mxMalloc(n_AScan * sizeof(double));
	mxSetPr(AScan_bench, pr);
	mxSetPi(AScan_bench, NULL);
    //pr2=pr;
	
    ///////////////throughput calc (assumption for Latency = Zero because working in saturation)
	n_X             = nVoxel_throughput;
	n_Y             = nVoxel_throughput;
    n_Z             = nVoxel_throughput;
		    
    //number of voxel
		n_IMAGE         = n_X * n_Y * n_Z;
					
		//////NOT Schlauchimage Code
		setImageDim[0]  = n_X;
		setImageDim[1]  = n_Y;        
		setImageDim[2]  = n_Z;    
   
    		out_bench        = mxCreateDoubleMatrix(0, 0, mxREAL);        //out     = mxCreateDoubleMatrix(n_Index,1,mxCOMPLEX);
			mxSetDimensions(out_bench, setImageDim, 3);                   //bsp. 3000x1  -> (3000,1) ,2
			pr = mxMalloc(n_IMAGE * sizeof(double));
			mxSetPr(out_bench, pr);
			mxSetPi(out_bench, NULL);
            //pr3=pr;
			
			image_sum_bench  = mxCreateDoubleMatrix(0, 0, mxREAL);        //out     = mxCreateDoubleMatrix(n_Index,1,mxCOMPLEX);
			mxSetDimensions(image_sum_bench, setImageDim, 3);                   //bsp. 3000x1  -> (3000,1) ,2
			pr = mxMalloc(n_IMAGE * sizeof(double));
			mxSetPr(image_sum_bench, pr);
			mxSetPi(image_sum_bench, NULL);
            //pr4=pr;
            
            
       //first call to get REAL memory from system (WARMUP later on faster)
      //set nCore!!!
	  nCores = 1;
      do{
      counter = TimeCounter();
      as2v_MT(mxGetPr(out_bench), mxGetPr(AScan_bench), n_AScan, mxGetPr(buffer_bench),
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, mxGetPi(AScan_bench), mxGetPi(
				     buffer_bench), mxGetPi(
				     out_bench), n_Y, n_Z, mxGetPr(image_sum_bench), mxGetPi(
				     image_sum_bench));
       counter2 = TimeCounter();
       } while (counter2<counter); //retry on error like used wrong core
                 
        //get througput time time for setup average
        bench_ref_time= counter2-counter; 
        //mexPrintf("benchreftime %llu c2:%llu c:%llu\n",bench_ref_time, counter, counter2);
        average = (uint32_t) ((minBenchTime/bench_ref_time)+1); //integer ceil
        
        if (average < minAverage) average = minAverage;  
        if (average > MAX_AVERAGE) average = MAX_AVERAGE;  
    
    
	for (i=1;i<=NUMCORES;i++) 
	{ 
	  //set nCore!!!
		nCores = i;
     
	    //benchmark throughput
        //mexCallMATLAB(0, NULL, 0, NULL, "tic");
        
		for (j=0;j<average;j++)
         {
         do{ counter=TimeCounter();
          //no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
		  as2v_MT(mxGetPr(out_bench), mxGetPr(AScan_bench), n_AScan, mxGetPr(buffer_bench),
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, mxGetPi(AScan_bench), mxGetPi(
				     buffer_bench), mxGetPi(
				     out_bench), n_Y, n_Z, mxGetPr(image_sum_bench), mxGetPi(
				     image_sum_bench));
          /*as2v_MT(pr3, pr2, n_AScan, pr1,
			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, NULL,NULL, NULL, n_Y, n_Z, pr4,NULL);*/
		  //mexCallMATLAB(1, &time_bench, 0, NULL, "toc");
          //mexPrintf("%llu ms\n",((TimeCounter()-counter)/1000000)); 
          counter2=TimeCounter();
         } while(counter2<counter); //retry on error like used wrong core
          average_buffer[j]= counter2-counter; 
        }        
                
        //bubblesort (small time top) for median
		  for (k=average-1;k>0;k--)
          {  for (l=average-1;l>0;l--){
                 if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;} 
             }
          }
        
         mittelw=0;
         for  (k=0;k<average;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/average;
         
         stdabw=0;
         for  (k=0;k<average;k++)        
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);} 
         stdabw = (uint64_t) sqrt( (int64_t)stdabw/(int64_t)average);  
           
         #ifdef addsig2vol_debug 
          mexPrintf("minTime: %f, timeCounter: %llu counter: %llu\n",(float)(uint64_t)bench_ref_time,counter2,counter);
           for (k=0;k<average;k++)               
           { mexPrintf("%llu average\n",average_buffer[k]);
           }
         #endif           
        
        //minimum selection
        bench_ref_time = (uint64_t) average_buffer[0];
        //median selection
        bench_ref_time = (uint64_t) average_buffer[(uint32_t) ceil(average/2)];
            
        //throughput[i-1] = (uint64_t) ((double)n_IMAGE/ ((double) ((uint32_t) bench_ref_time)/1000)); //gekuerzt ns und MVoxel ->Mvoxel/s // index 0 = core1; index1 = core 2 etc...(!)
        throughput[i-1] =  ((1000*(uint64_t)n_IMAGE)/ bench_ref_time); //gekuerzt ns und MVoxel ->Mvoxel/s // index 0 = core1; index1 = core 2 etc...(!)
  	
        mexPrintf("%7i |%7llu |%7llu |%7llu |%7llu |%7llu", i, (int64_t)throughput[i-1],((1000*(int64_t)n_IMAGE)/(int64_t)mittelw), (1000*(int64_t)n_IMAGE)/((int64_t)average_buffer[average-1]), (1000*(int64_t)n_IMAGE)/((int64_t)average_buffer[0]), (1000*(int64_t)n_IMAGE)/((int64_t)stdabw));    
  	
        
       //get throughput time time for setup average
       average = (uint32_t) ((minBenchTime/bench_ref_time)+1); //in ns; +1=integer-ceil
       if (average < minAverage) average = minAverage;  
       if (average > MAX_AVERAGE) average = MAX_AVERAGE;  
       
		 // printf("%e MVoxel/s for %i threads, %e sec through\n",(double)(4*n_IMAGE)/(*mxGetPr(time_bench)),n_IMAGE,*mxGetPr(time_bench));
      //  	  printf("%f",(float)(((4*2000000)/(*(double*)mxGetPr(time_bench)))));
	  
//         //first call to get REAL memory form system (later on faster)
//       as2v_MT(mxGetPr(out_bench), mxGetPr(AScan_bench), n_AScan, mxGetPr(buffer_bench),
// 			     &pix_vec_bench[0], n_X, &rec_vec_bench[0], &send_vec_bench[0],
// 			     &float_bench,
// 			     &float_bench, &float_bench, mxGetPi(AScan_bench), mxGetPi(
// 				     buffer_bench), mxGetPi(
// 				     out_bench), n_Y, n_Z, mxGetPr(image_sum_bench), mxGetPi(
// 				     image_sum_bench));
	
	   ///////////////////////////////////////////////////////benchmark latency
	//	 mexCallMATLAB(0, NULL, 0, NULL, "tic");
       	  
		for (j=0;j<minAverage;j++)
         {
          do{
          counter=TimeCounter();
			//no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
		  as2v_MT(mxGetPr(out_bench), mxGetPr(AScan_bench), n_AScan, mxGetPr(buffer_bench),
			     &pix_vec_bench[0],(unsigned int) MIN_VOXEL, &rec_vec_bench[0], &send_vec_bench[0],
			     &float_bench,
			     &float_bench, &float_bench, mxGetPi(AScan_bench), mxGetPi(
				     buffer_bench), mxGetPi(
				     out_bench), (unsigned int) 1,(unsigned int) NUMCORES, mxGetPr(image_sum_bench), mxGetPi(
				     image_sum_bench));
				     				 
		//  mexCallMATLAB(1, &time_bench, 0, NULL, "toc");
        //latency[i] =(*mxGetPr(time_bench)/minAverage)/1;	 //assumption for 1 PIXEL!
		 counter2 = TimeCounter(); } while(counter2<counter); //retry on error like used wrong core
         average_buffer[j]= counter2-counter; 
        }     
        //bubblesort (small time top)
		 for (k=minAverage-1;k>0;k--)
         {  for (l=minAverage-1;l>0;l--){
               if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;} 
            }
         }    
         #ifdef addsig2vol_debug 
            for (k=0;k<minAverage;k++)               
           { mexPrintf("%llu average\n",average_buffer[k]);
           }
        #endif  
        
         mittelw=0;
         for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/minAverage;
         
         stdabw=0;
         for  (k=0;k<minAverage;k++)        
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);} 
         stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);  
       
        //minimum selection
        latency[i-1] = (uint64_t) average_buffer[0];	 //assumption for 1 PIXEL! // index 0 = core1; index1 = core 2 etc...(!)
        //median selection
        latency[i-1] = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
	    mexPrintf("|%10llu |%10llu |%10llu |%10llu |%8i\n", latency[i-1],average_buffer[0],average_buffer[minAverage-1],stdabw,average);    
  		}         
 
		for (i=NUMCORES;i>0;i--)
        { throughput_sort_buffer[i-1]=throughput[i-1];
		}
		    
      //bubblesort biggest on top (bottom up)
      for (i=NUMCORES-1;i>=1;i--)
      { if (throughput_sort_buffer[i]>(throughput_sort_buffer[i-1])) throughput_sort_buffer[i-1] = throughput_sort_buffer[i]; } 

      //find core-number according to perf.-value
      for (i=0;i<NUMCORES;i++)
      { if (throughput[i] == throughput_sort_buffer[0]) break; } //break on first core-number which fits 
     
      //set up used Cores
      nCores = i+1;      
      //backup in nCores_bench
      nCores_bench = i+1; 
                  
      switch (i+1)
      {
      case 1:
          mexPrintf("Detected Single-core System, 1 thread prefered\n"); 
          break;
      case 2:
          mexPrintf("Detected Dual-core or Hyperthreading System, 2 threads prefered\n");
          break;
      case 3:
          mexPrintf("Detected Triple-core or Hyperthreading System, 3 threads prefered\n");
          break;
      case 4:
          mexPrintf("Detected Quad-Core system, 4 threads prefered\n");
          break;
      case 8:
          mexPrintf("Detected Octa-Core system (or Quadcore with HT), 8 threads prefered\n");
          break;
      case 16:
          mexPrintf("Detected 16 Core system, 16 threads prefered\n");
          break;
      default: 
          mexPrintf("Detected %i-core system, %i threads prefered (?)\n",i+1,i+1);
      }      
          
      
      //benchmark perf-per size
      mexPrintf("\nPerformance for various imagesize in Voxel (with potentially %i Cores)\n",nCores);
      mexPrintf("     Voxel | Throughput in kVoxel/s         | Time in mikros  | Malloc time (mikro-sec)\n");
      mexPrintf("           | Median   | Mean     | Std      |           | Median | mean  | Std  | min  | max  \n"); 
      mexPrintf("------------------------------------------------------------------------------------------------\n");
      for (i=MIN_VOXEL;i<=(nVoxel_throughput*nVoxel_throughput*nVoxel_throughput);i=i*2)
      {          
          for (j=0;j<minAverage;j++)
          {
              do {
                  counter=TimeCounter();
                  //no sizeof(double) needed because compilers assumes already double as datatype for pointer!!!
                  as2v_MT(mxGetPr(out_bench), mxGetPr(AScan_bench), n_AScan, mxGetPr(buffer_bench),
                  &pix_vec_bench[0], (uint32_t) MIN_VOXEL, &rec_vec_bench[0], &send_vec_bench[0],
                  &float_bench,
                  &float_bench, &float_bench, mxGetPi(AScan_bench), mxGetPi(
                  buffer_bench), mxGetPi(
                  out_bench), (uint32_t) 1, (uint32_t) (i/MIN_VOXEL),  mxGetPr(image_sum_bench), mxGetPi(
                  image_sum_bench));
                  counter2 = TimeCounter(); } while(counter2<counter); //retry on error like used wrong core
              average_buffer[j]= counter2-counter;
          }
          //bubblesort (small time top)
          for (k=minAverage-1;k>0;k--)
          {  for (l=minAverage-1;l>0;l--){
                     if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;} 
                 }
          }          
         mittelw=0;
         for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
         mittelw=(uint64_t) mittelw/minAverage;
         
         stdabw=0;
         for  (k=0;k<minAverage;k++)        
         {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);} 
         stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);            
          
          #ifdef addsig2vol_debug 
           for (k=0;k<minAverage;k++)               
           { mexPrintf("%llu values-range \n",average_buffer[k]);}
          #endif                   
          
          //minimum selection
          counter2 = (uint64_t) average_buffer[0];	 //assumption for 1 PIXEL! // index 0 = core1; index1 = core 2 etc...(!)
          //median selection
          counter2 = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
          // mexPrintf("%10luu , %10llu, %10llu \n",TimeCounter, counter, TimeCounter());
          //counter=pow(2,64)-counter;
          
          //mexCallMATLAB(1, &time_bench, 0, NULL, "toc");
          //mexPrintf("%llu ms\n",((TimeCounter()-counter)/1000000));
          mexPrintf("%10i | %8llu | %8llu | %8llu |%8llu ",i,(1000000*(uint64_t)i)/(counter2),(1000000*(uint64_t)i)/(mittelw),((uint64_t)i*1000000)/(stdabw),(uint64_t)(counter2/1000) );
          
          //fix  (UGLY!!!!)
          nCores=nCores_bench;          
          
          //benchmark mem-alloc
          minAverage=100;
          
              for (j=0;j<minAverage;j++)
              {
                  do {
                      mxDestroyArray(image_sum_bench);
                      
                      counter=TimeCounter();
                      setImageDim[0]  = i;
                      setImageDim[1]  = 1;
                      setImageDim[2]  = 1;                

                      image_sum_bench  = mxCreateDoubleMatrix(0, 0, mxREAL);        //out     = mxCreateDoubleMatrix(n_Index,1,mxCOMPLEX);
                      mxSetDimensions(image_sum_bench, setImageDim, 3);                   //bsp. 3000x1  -> (3000,1) ,2
                      pr = mxMalloc(n_IMAGE * sizeof(double));
                      mxSetPr(image_sum_bench, pr);
                      mxSetPi(image_sum_bench, NULL);
                     counter2 = TimeCounter(); } while(counter2<counter); //retry on error like used wrong core
                     average_buffer[j]= counter2-counter;
              }
              //bubblesort (small time top)
              for (k=minAverage-1;k>0;k--)
              {  for (l=minAverage-1;l>0;l--){
                  if (average_buffer[l]<average_buffer[l-1]) {counter2=average_buffer[l-1]; average_buffer[l-1] = average_buffer[l]; average_buffer[l]=counter2;}
              }
              }
              
              mittelw=0;
             for  (k=0;k<minAverage;k++)  {mittelw=average_buffer[k]+mittelw;}
             mittelw=(uint64_t) mittelw/minAverage;

             stdabw=0;
             for  (k=0;k<minAverage;k++)        
             {if (mittelw>average_buffer[k]) stdabw=stdabw+(mittelw-average_buffer[k])*(mittelw-average_buffer[k]); else stdabw=stdabw+(average_buffer[k]-mittelw)*(average_buffer[k]-mittelw);} 
             stdabw = (uint64_t) sqrt((int64_t) stdabw/minAverage);            

             //minimum selection
             counter2 = (uint64_t) average_buffer[0];	 //assumption for 1 PIXEL! // index 0 = core1; index1 = core 2 etc...(!)
             //median selection
             counter2 = (uint64_t) average_buffer[(uint32_t) ceil(minAverage/2)];
             mexPrintf("| %8llu | %8llu | %8llu |%8llu |%8llu\n", (uint64_t)counter2/1000,(uint64_t)mittelw/1000,(uint64_t)stdabw/1000,(uint64_t) average_buffer[0]/1000,(uint64_t) average_buffer[minAverage-1]/1000);
              
          }     
                  
      
        //check status
        fpu_check();        
      
      
        //free memory
		  mxDestroyArray(out_bench);
		  mxDestroyArray(image_sum_bench);
        //free memory
		  mxDestroyArray(buffer_bench);
		  mxDestroyArray(AScan_bench);
		  mxDestroyArray(time_bench);        
}


uint64_t TimeCounter(void) {
uint64_t counter;  
    #ifndef __WIN32__
        #ifdef WIN32
        #define __WIN32__
        #endif
        #ifdef _WIN32
        #define __WIN32__
        #endif
        #ifdef __WIN32
        #define __WIN32__
        #endif
        #ifdef _WIN64
        #define __WIN32__
        #endif
        #ifdef WIN64
        #define __WIN32__
        #endif
        #ifdef _WINDOWS
        #define __WIN32__
        #endif
    #endif
            
    #ifdef __WIN32__ //fitting windows64 AND windows32    
    #include <windows.h>
    //register uint64_t temp; 
    uint64_t iFreq, iCount;
    QueryPerformanceFrequency((LARGE_INTEGER*)&iFreq);
    QueryPerformanceCounter((LARGE_INTEGER*)&iCount); 
    //counter = (uint64_t) (1000000000*((double)iCount/(double)iFreq)); //f�r nSekunden (balancing der multiplikation)
    counter = (uint64_t) ((1000000000*iCount)/iFreq); //f�r nSekunden (balancing der multiplikation)
    
//#ifdef addsig2vol_debug 
    //mexPrintf("iFreq:%llu, iCount:%llu, Counter:%llu\n",iFreq,iCount,counter); 
    //#endif
    #endif   
    
   #ifdef __linux__       
        
    #include <time.h>
    #include <unistd.h>

    struct timespec time_str;  
    struct timespec time_res;    
    clock_gettime(CLOCK_REALTIME, &time_str);
    counter = (uint64_t) time_str.tv_nsec + (uint64_t) (time_str.tv_sec*1000000000);
    
    /*
    clock_gettime(CLOCK_REALTIME, &time_str); clock_getres(CLOCK_REALTIME, &time_res);
    mexPrintf("CLOCK_REALTIME clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec); 
    clock_gettime(CLOCK_MONOTONIC, &time_str); clock_getres(CLOCK_MONOTONIC, &time_res);
    mexPrintf("CLOCK_MONOTONIC clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec); 
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_str); clock_getres(CLOCK_PROCESS_CPUTIME_ID, &time_res);
    mexPrintf("CLOCK_PROCESS_CPUTIME_ID clockRes:%f, nsec:%f, csec:%f\n",(double)time_res.tv_nsec,(double) time_str.tv_nsec,(double) time_str.tv_sec); 
    */
      //int clock_getres(CLOCK_REALTIME, struct timespec *res);
    //int clock_gettime(clockid_t clk_id, struct timespec *tp);
    #endif   
    return counter;
 }

void fpu_check()
{
// http://www.christian-seiler.de/projekte/fpmath/ 
    #ifdef __WIN32__
    //#include <float.h>
    uint32_t control_wordfp=0, status_wordfp=0;
    
    //read out
    control_wordfp = _controlfp(0, 0);//,&control_wordfp, 0);
    control_wordfp = getFPUStateX86();
    //http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz/
    //DAZ and FTZ in MXCRS
    //FTZ: The FTZ bit (bit 15)& The underflow exception (bit 11) 
    //DAZ: DAZ bit (bit 6) 
    //control_wordfp = _controlfp( _CW_DEFAULT, 0xfffff);
    //control_wordfp =_controlfp(_DN_FLUSH|32 ,_MCW_DN|32 ); //_controlfp(control_wordfp | 34816, 4294967295);
    //control_wordfp =_controlfp(_PC_64, _MCW_PC); //_PC_24, _PC_53 _PC_64
    //control_wordfp =_controlfp(_RC_NEAR, _MCW_RC); //_RC_UP _RC_CHOP _RC_DOWN _RC_NEAR
 
    status_wordfp = _statusfp(); 
    _clearfp();
    #else
    #include <fpu_control.h>
    fpu_control_t control_wordfp=0, fpu_cw=0;
    _FPU_GETCW(control_wordfp);
    
    //Status WORD under linux????

    #endif 
            
            
    mexPrintf("FPU/SSE Control-register(0x%.4x): ", control_wordfp );
    //RC round-control
     switch ((control_wordfp & __FPU_CW_ROUND_MASK__)>>0)//& 3072) 
    {case __FPU_CW_ROUND_NEAR__:
         mexPrintf("nearest rounding");
         break;
     case __FPU_CW_ROUND_UP__:
          mexPrintf("ceil-rounding");
         break;
     case __FPU_CW_ROUND_DOWN__:
          mexPrintf("floor-rounding");
         break;
     case __FPU_CW_ROUND_CHOP__:
          mexPrintf("truncation rounding");
         break;
    }
    // mexPrintf("%x %x %x",control_wordfp & _MCW_PC, _MCW_PC,_MCW_RC);

    //PC Precision-control
    switch ((control_wordfp & __FPU_CW_PREC_MASK__ )>>0 )//& 768) 
    {case __FPU_CW_PREC_SINGLE__:
         mexPrintf(", internal precision float (32bit)\n");
         break;
     case __FPU_CW_PREC_DOUBLE__:
         mexPrintf(", internal precision double (64bit)\n");
         break;
     case __FPU_CW_PREC_EXTENDED__:
         mexPrintf(", internal precision extended double (80bit)\n"); 
         break;
     default : 
         mexPrintf(", internal precision invalid \n"); 
    }       
    
#ifdef __WIN32__
    //FPU/SSE status
    mexPrintf("FPU/SSE status (0x%.4x): ", status_wordfp);
    if ((status_wordfp & 32)>>5)  mexPrintf("inexact ");
    if ((status_wordfp & 16)>>4)  mexPrintf("underflow ");
    if ((status_wordfp &  8)>>3)  mexPrintf("overflow ");
    if ((status_wordfp &  4)>>2)  mexPrintf("division-by-zero ");
    if ((status_wordfp &  2)>>1)  mexPrintf("denormal ");  
    if ((status_wordfp &  1)>>0)  mexPrintf("invalid operation mask ");
    if ( status_wordfp == 0 )  mexPrintf("OK");   
    mexPrintf("\n");
    #endif
}

