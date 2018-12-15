// Automatching.cpp : 
//

//Doing the rough matching

#include "stdafx.h"
#include "omp.h"
//#include <stdio.h>
#include <ctime> 
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm> //
//#include <CString>
#include "nabo/nabo.h"
#include "nabo/Eigen/Eigenvalues"
#include <windows.h>
#include <comdef.h>
using namespace std;
using namespace Nabo;
using namespace Eigen;


#define ITERMAX 200
#define ITERMIN 1
#define pi 3.1415926535

typedef struct {
		Matrix3d RT;      //3*3 rotation matrix
		Vector3d TR;      //3*1 translate matrix
		double* data_out;//target print out data
		double res;      //error
		int  iter;   //iter //for rough_matching it means the correct index of the eight PCA axis direction
	}ICPOUT;

__int64 CountLines(char *filename)  ;
int ReadTXT(char* fileName, double* matData);
//int Readparam(char* fileName, float* trig,int* knum,float* anglethrold,float* scal) ;
//__int64 fextract(double* gpoint,__int64 pnum ,double*gfea,double trig,int knum,double anglethrold); //return feature num
void RoughMaching(double* target,double* model,__int64 numt,__int64 numm,ICPOUT* result_rough);
__int64 linefextract(double* gpoint,__int64 pnum ,double*gfea,double trig,int knum,double anglethrold);
double mulknn(double* target,double* model,__int64 trow,__int64 mrow,double* dataclose,int k=1);
void weldcorrect(double* weld,double* fea,double* data, __int64 wrow, __int64 frow, __int64 drow,double* newweld);
void icp_alg(double* model,__int64 mrow,double* target,__int64 trow,ICPOUT* result,double res,int maxiter,int miniter);
void rotationToQuaternion(double* RT,double* Q);
void quaternionToRotation(double* Q,double* R);
void extrapolateInTransformSpace(int t,double* dq,double* qs,double* err);
double extrapolate(double* v,double* d,double vmax);
void datafit(double *x,double *y,int n,double *a,int m,double *dt);
int Readparam(char* fileName, float* trig, int* knum, float* anglethrold, float* scal);

int main(int argc, char* argv[])
{
	//char* Tpath;
	//if(argv[1]!=NULL)
	//{
	//	_bstr_t buf(argv[1]);
	//    Tpath=buf;
	//}
	//else
	//{
	//    //cout<<"add the txt path"<<endl;
	//	//cin.get(Tpath,100);
	//	//cin>>Tpath;
	//}

	//**********   change file name,change other data**********//
	//**********// model is moved,data is static //**********//
	char* modelpath = "model\\CADmodel.txt";
	char* datapath = "data\\data.txt";
	char* CADweldpath = "model\\CADweld.txt";
	int i = 1;
	bool DOLFEA = false;
	bool DOICP = false;	
	bool USEWELDDOICP = false;
	while ( i < argc) {
		//cout << *argv[i] << endl;
		if (!strcmp(argv[i], "-m")) {
			modelpath = argv[++i];
		}
		else if (!strcmp(argv[i], "-d")) {
			datapath = argv[++i];
		}
		else if (!strcmp(argv[i], "-lf")) {
			DOLFEA = true; //Do line feature extraction 
		}
		else if (!strcmp(argv[i], "-icp")) {
			DOICP = true; //DO icp
			modelpath = argv[++i];
			CADweldpath = argv[++i];
		}
		else if (!strcmp(argv[i], "-wicp")) {
			USEWELDDOICP = true; //Use weld curve DO ICP
			modelpath = argv[++i];
			CADweldpath = argv[++i];
		}
		else if (argv[i][0] == '-') {
			cerr << "Unknown flag\n";
			return -1;
		}
		i++;
		
	}
	Eigen::initParallel();
	//Eigen::setNbThreads(6);

	clock_t start, stop;
	double icp_res=0.0001;
	__int64 maxiter=ITERMAX, miniter=ITERMIN;

	//char* curve="data\\textmodelcurve.txt";
	//********************************************//
	cout << endl << "reading model points num ..." << endl;
	__int64 num_m = CountLines(modelpath);//
	cout << " numModel: "<<num_m << endl;
	cout << endl << "reading data points num..." << endl;
	__int64 num_t = CountLines(datapath);//
	cout << " numData: "<<num_t << endl;
	
	double* model=(double*) calloc(3*num_m,sizeof(double));
	double* data=(double*) calloc(3*num_t,sizeof(double));
	double* T_rough=(double*) calloc(3* num_m,sizeof(double));

	 
	ICPOUT* result_rough=new ICPOUT;//C++ struct inition 
	result_rough->RT=Matrix3d::Identity();
	result_rough->TR=Vector3d::Zero();
	result_rough->data_out=T_rough;
	result_rough->iter=0;
	result_rough->res=0.0;
	cout << endl << "reading model points ..." << endl;
	ReadTXT(modelpath, model);
	cout << endl << "reading data points ..." << endl;
	ReadTXT(datapath, data);

	ICPOUT* result_ICP = new ICPOUT;//C++ struct inition 
	result_ICP->RT = Matrix3d::Identity();
	result_ICP->TR = Vector3d::Zero();
	//result_ICP->data_out = T_rough;
	result_ICP->iter = 0;
	result_ICP->res = 0.0;
	__int64 num_c = 0;
	double* CADweld;
	cout << endl << "reading CADweld points num ..." << endl;
	num_c = CountLines(CADweldpath);//
	cout << " numCADweld: " << num_c << endl;
	CADweld = (double*)calloc(3 * num_c, sizeof(double));
	cout << endl << "reading CADweld points ..." << endl;
	ReadTXT(CADweldpath, CADweld);
	if (USEWELDDOICP)
	{			
		result_ICP->data_out = (double*)calloc(3 * num_c, sizeof(double));
	}
	else
	{
		//CADweld = NULL;
		result_ICP->data_out = (double*)calloc(3 * num_m, sizeof(double));
	}
	double* gfea = (double*)calloc(3 * num_t, sizeof(double));
	__int64 feanum = 0;
	float scal = 1.0; float trig; int knum; float anglethrold;

	cout << "reading param for caculation" << endl;
	Readparam("param\\mparam.param", &trig, &knum, &anglethrold, &scal);

	//******* Adjust the point cloud ratio here ********//
	Map<MatrixXd> CADMODEL(model, 3, num_m);
	Map<MatrixXd> CADWELD(CADweld, 3, num_c);
	Vector3d m_cadmodel = CADMODEL.rowwise().mean();
	//Vector3d m_cadweld = CADWELD.rowwise().mean();
	CADMODEL = (CADMODEL-m_cadmodel*MatrixXd::Ones(1,num_m))*scal+ m_cadmodel*MatrixXd::Ones(1, num_m);
	CADWELD = (CADWELD - m_cadmodel*MatrixXd::Ones(1, num_c))*scal+ m_cadmodel*MatrixXd::Ones(1, num_c);
	//start matching
	start = clock();
	//******* rough matching doing here ********//
	
	cout << endl << "Doing the Rough Maching..." << endl;
	RoughMaching(model, data, num_m, num_t,result_rough);
	cout << endl << "Rough Maching Done!" << endl;

	//******* feature extraction doing here ********//


	if (DOLFEA)
	{		
		cout << endl << "Running feature extract from data pointcloud..." << endl;		
		feanum = linefextract(data, num_t, gfea, trig, knum, anglethrold);		
		cout << "feature point num: " << feanum << endl;
	}
	Map<MatrixXd> GFEA(gfea, 3, feanum);
	//********* ICP doing here *****************//
		
	if (DOICP || USEWELDDOICP)
	{
		if (USEWELDDOICP && DOLFEA)
		{			
			Map<MatrixXd> cadweld(CADweld, 3, num_c);
			cadweld = result_rough->RT*cadweld + result_rough->TR*MatrixXd::Ones(1, num_c);
			cout << endl << "Running ICP with CADweld and data feature extraction (point-to-point, no outliers)..." << endl;
			icp_alg(gfea, feanum, CADweld, num_c, result_ICP, icp_res, maxiter, miniter);
		}
		else if (USEWELDDOICP)
		{
			Map<MatrixXd> cadweld(CADweld, 3, num_c);
			cadweld = result_rough->RT*cadweld + result_rough->TR*MatrixXd::Ones(1, num_c);
			cout << endl << "Running ICP with CADweld and data (point-to-point, no outliers)..." << endl;
			icp_alg(data, num_t, CADweld, num_c, result_ICP, icp_res, maxiter, miniter);
		}
		else if (DOLFEA) 
		{
			
			cout << endl << "Running ICP with model and data feature extraction (point-to-point, no outliers)..." << endl;
			icp_alg(gfea, feanum, result_rough->data_out, num_m, result_ICP, icp_res, maxiter, miniter);
		}
		else
		{
			cout << endl << "Running ICP with model and data (point-to-point, no outliers)..." << endl;
			icp_alg(data, num_t, result_rough->data_out, num_m, result_ICP, icp_res, maxiter, miniter);
		}
	}
	//******************************************//

	stop=clock();
	printf("Run time = %g s \n",((double)(stop - start)) / CLOCKS_PER_SEC);
	

	char* resultRTpath="result\\result_RT.txt";

	Matrix4d resultRT= Matrix4d::Zero();
	resultRT.block<3, 3>(0, 0) = result_ICP->RT*result_rough->RT;
	resultRT.block<3, 1>(0, 3) = result_ICP->RT*result_rough->TR + result_ICP->TR;
	resultRT(3,3) = 1;
	ofstream outputTXT;
	cout << "Save resultRT to ( result\\result_RT.txt )" << endl;
	outputTXT.clear();
	outputTXT.open(resultRTpath, fstream::in | fstream::trunc);
	outputTXT<< resultRT<<endl;
	outputTXT.close();	
	cout << "result_RT: " << endl << resultRT << endl;
	if (DOLFEA)
	{
		cout << "Save the feature points to ( result\\featurePC.txt )" << endl;
		outputTXT.clear();
		outputTXT.open("result\\featurePC.txt", fstream::in | fstream::trunc);
		outputTXT << GFEA.transpose() << endl;
		outputTXT.close();
	}
	if (USEWELDDOICP)
	{
		//find weld
		//CADweld has changed before the ICP 
		Map<MatrixXd> cadweldmatched(result_ICP->data_out, 3, num_c);
		
		cout << "Save the CADweld points after matched to ( result\\CADweldmatched.txt )" << endl;
		outputTXT.clear();
		outputTXT.open("result\\CADweldmatched.txt", fstream::in | fstream::trunc);
		outputTXT << cadweldmatched.transpose() << endl;
		outputTXT.close();

		cout << "Save the weld curve points to ( result\\weldcurve.txt )" << endl;		
		double *weldcurve = (double*)calloc(3 * num_c, sizeof(double));
		//mulknn(result_ICP->data_out, gfea, num_c, feanum, weldcurve, 1);
		if(DOLFEA)
			weldcorrect(result_ICP->data_out, gfea, data, num_c, feanum, num_t, weldcurve);
		else
			mulknn(result_ICP->data_out, data, num_c, num_t, weldcurve, 1);
		Map<MatrixXd> WELDCURVE(weldcurve, 3, num_c); 
		outputTXT.clear();
		outputTXT.open("result\\weldcurve.txt", fstream::in | fstream::trunc);
		outputTXT << WELDCURVE.transpose() << endl;
		outputTXT.close();

		Map<MatrixXd> CADMODELMATCHED(model, 3, num_m);
		CADMODELMATCHED = resultRT.block<3, 3>(0, 0)*CADMODELMATCHED + resultRT.block<3, 1>(0, 3)*MatrixXd::Ones(1, num_m);
		cout << "Save the changed CADmodel points to ( result\\CADmodelmatched.txt )" << endl;
		outputTXT.clear();
		outputTXT.open("result\\CADmodelmatched.txt", fstream::in | fstream::trunc);
		outputTXT << CADMODELMATCHED.transpose() << endl;
		outputTXT.close();

		free(weldcurve);
	}
	else if (!USEWELDDOICP && DOICP)
	{
		Map<MatrixXd> CADMODELMATCHED(result_ICP->data_out, 3, num_m);
		
		cout << "Save the CADmodel points after ICP matched to ( result\\CADmodelmatched.txt )" << endl;
		outputTXT.clear();
		outputTXT.open("result\\CADmodelmatched.txt", fstream::in | fstream::trunc);
		outputTXT << CADMODELMATCHED.transpose() << endl;
		outputTXT.close();

		//CADweld has not changed before the ICP 
		Map<MatrixXd> cadweldmatched(CADweld, 3, num_c);
		cadweldmatched = resultRT.block<3, 3>(0, 0)*cadweldmatched + resultRT.block<3, 1>(0, 3)*MatrixXd::Ones(1, num_c);
		cout << "Save the changed CADweld points to ( result\\CADweldmatched.txt )" << endl;
		outputTXT.clear();
		outputTXT.open("result\\CADweldmatched.txt", fstream::in | fstream::trunc);
		outputTXT << cadweldmatched.transpose() << endl;
		outputTXT.close();

		cout << "Save the weld curve points  to ( result\\weldcurve.txt )" << endl;
		double *weldcurve = (double*)calloc(3 * num_c, sizeof(double));
		//mulknn(CADweld, gfea, num_c, feanum, weldcurve, 1);
		if (DOLFEA)
			weldcorrect(CADweld, gfea, data, num_c, feanum, num_t, weldcurve);
		else
			mulknn(CADweld, data, num_c, num_t, weldcurve, 1);
		Map<MatrixXd> WELDCURVE(weldcurve, 3, num_c);
		outputTXT.clear();
		outputTXT.open("result\\weldcurve.txt", fstream::in | fstream::trunc);
		outputTXT << WELDCURVE.transpose() << endl;
		outputTXT.close();

		free(weldcurve);
	}
	else
	{
		Map<MatrixXd> CADMODELMATCHED(result_rough->data_out, 3, num_m);

		cout << "Save the CADmodel points after rough matched to ( result\\CADmodelmatched.txt )" << endl;
		outputTXT.clear();
		outputTXT.open("result\\CADmodelmatched.txt", fstream::in | fstream::trunc);
		outputTXT << CADMODELMATCHED.transpose() << endl;
		outputTXT.close();
	}
	free(result_ICP->data_out);	
	free(model); free(data); free(T_rough); free(gfea); free(CADweld);
	//free(modelpath); free(datapath);
	delete result_rough;
	delete result_ICP;
	//system("pause");
	return 0;
}
__int64 CountLines(char *filename)  
{  
    ifstream ReadFile;  
    __int64 n=0;  
    string tmp;  
    ReadFile.open(filename,ios::in);//ios::in means file can read only  
    if(ReadFile.fail())//file openned fail, return -1 
    {  
		printf("Can not open the file!\n");
        return -1;  
    }  
    else//file exist  
    {  
        
		while(getline(ReadFile,tmp,'\n'))  
        {  
            n++;  
        }  
        ReadFile.close();  
        return n;  
    }  
}
int ReadTXT(char* fileName, double* matData)  
{  
	FILE *file;  
    __int64 LINES;    
    if(fopen_s(&file,fileName,"r"))
      {
           printf("Can not open the file!\n");
           return -1;
      } 
    else//file exist
    {  
        LINES=CountLines(fileName);  
                 
        int i;  
		for(i=0;i<LINES*3;i++) 
        {  
			fscanf_s(file,"%lf",&matData[i]);
			
        }  
		//cout<<matData[LINES*3-1]<<endl;
		fclose(file); //close file  
	return 0;
    }  
	
}
void RoughMaching(double* target,double* model,__int64 numt,__int64 numm,ICPOUT* result_rough)
{
	//model means the static pointcloud
	//target is the moved pointcloud
	Map<MatrixXd> MODEL(model,3,numm);
	Map<MatrixXd> TARGET(target,3,numt);
	Map<MatrixXd> RESULT(result_rough->data_out,3,numt);

	Vector3d m_model=MODEL.rowwise().mean();//(MODEL.rowwise().maxCoeff()+MODEL.rowwise().minCoeff())/2;
	Vector3d m_target=TARGET.rowwise().mean();//(TARGET.rowwise().maxCoeff()+TARGET.rowwise().minCoeff())/2;
	// caculate Covariance matrix
	MatrixXd model_anti=(MODEL-m_model*VectorXd::Ones(numm).transpose())*(MODEL-m_model*VectorXd::Ones(numm).transpose()).transpose();
	MatrixXd target_anti=(TARGET-m_target*VectorXd::Ones(numt).transpose())*(TARGET-m_target*VectorXd::Ones(numt).transpose()).transpose();
	//SVD: Different with matlab,but 
	JacobiSVD<MatrixXd> modelsvd(model_anti, ComputeFullU | ComputeFullV);//Different with matlab
	JacobiSVD<MatrixXd> targetsvd(target_anti, ComputeFullU | ComputeFullV);//Different with matlab
	
	MatrixXd vc_target=targetsvd.matrixU().transpose()*(TARGET-m_target*VectorXd::Ones(numt).transpose());
	MatrixXd vc_model=modelsvd.matrixU().transpose()*(MODEL-m_model*VectorXd::Ones(numm).transpose());
	
	MatrixXd quadbuf(8,3);
	quadbuf<<1,1,1,  //the eight possibilities of PCA axes directions
			1,1,-1,
			1,-1,1,
			1,-1,-1,
			-1,1,1,
			-1,1,-1,
			-1,-1,1,
			-1,-1,-1;
	//check the eight root mean square error
	VectorXd disbuf(8);
	MatrixXd dataclose(3,numt);
	for(int i=0;i<8;++i)
	{
		Matrix3d quadiag(quadbuf.row(i).asDiagonal());
		MatrixXd vc_buf=quadiag*vc_target;
		disbuf(i)=mulknn(vc_buf.data(),vc_model.data(),numt,numm,dataclose.data(),1);
		
	}
	cout << disbuf << endl;
	int minindex;
	disbuf.minCoeff(&minindex);
	Matrix3d realquadiag(quadbuf.row(minindex).asDiagonal());
	result_rough->RT=modelsvd.matrixU()*realquadiag*targetsvd.matrixU().transpose();
	result_rough->TR=m_model-result_rough->RT*m_target;
	RESULT=result_rough->RT*TARGET+result_rough->TR*VectorXd::Ones(numt).transpose();
	result_rough->res=disbuf(minindex);
	result_rough->iter=minindex;
	//cout<<disbuf<<endl;
	
}
__int64 linefextract(double* gpoint,__int64 pnum ,double*gfea,double trig,int knum,double anglethrold)
{
	//for line laser scaned data
	Map<MatrixXd> gp(gpoint,3,pnum);
	Map<MatrixXd> GFEA(gfea,3,pnum);
	MatrixXd GPOINT=gp.transpose();
	double ginit=GPOINT(0,0);
	__int64 gbegin=1,gend=1,k=0;
	MatrixXd gline(pnum,2);
	Vector2d veold,venew;
	for(__int64 n=0;n<pnum;++n)
	{
		if(abs(GPOINT(n,0)-ginit)>=trig)
		{
		    ginit=GPOINT(n,0);
			gend=n-1; gline(k,0)=gbegin; gline(k,1)=gend; gbegin=n;
			++k;
		}
	}
	gline(k,0)=gbegin; gline(k,1)=pnum-1;

    __int64 feanum=0;
	for(__int64 i=0;i<k+1;++i)
	{
		if(abs(gline(i,1)-gline(i,0))<knum*2) continue;

		for(__int64 j=(gline(i,0)+knum);j<(gline(i,1)-knum);++j)
		{
			veold=(gp.block(1,j,2,1)*MatrixXd::Ones(1,knum)-gp.block(1,j-knum,2,knum)).rowwise().mean();
			venew=(gp.block(1,j+1,2,knum)-gp.block(1,j,2,1)*MatrixXd::Ones(1,knum)).rowwise().mean();

			double nveold=veold.norm();
			double nvenew=venew.norm();
			//if((nveold+nvenew)/2 >knum ) continue;//Filter the outlier points 
			
			double gangle=acos( veold.dot(venew)/nvenew/nveold)*180/pi;
			if(gangle>anglethrold)
			{ 
				GFEA.col(feanum)=gp.col(j);
				++feanum;
			}
			
		}
	
	}

	return feanum;
}
double mulknn(double* target,double* model,__int64 trow,__int64 mrow,double* dataclose,int k)
{
	//target to model
	Map<MatrixXd> Mt(model,3,mrow);
	Map<MatrixXd> Tt(target,3,trow);
	MatrixXd M=Mt;  
	MatrixXd T=Tt;
	NNSearchD* nns = NNSearchD::createKDTreeLinearHeap(M);
	MatrixXi indices(k, T.cols());
	MatrixXd dists2(k, T.cols());
	nns->knn(T, indices, dists2, k, 0, NNSearchF::SORT_RESULTS);

	Map<MatrixXd> DATACLOSE(dataclose,3,trow);
	for(__int64 i=0;i<trow;++i)
		DATACLOSE.col(i)=Mt.col(indices(0,i));

	delete nns;
	return dists2.row(0).cwiseSqrt().mean();//return the average dis of nearest points
}
void weldcorrect(double* weld, double* fea, double* data, __int64 wrow, __int64 frow, __int64 drow, double* newweld)
{
	//int k = 1;
	Map<MatrixXd> WELD(weld,3,wrow);
	Map<MatrixXd> FEA(fea,3,frow);
	Map<MatrixXd> DATA(data, 3, drow);
	MatrixXd Weld = WELD;
	MatrixXd Fea = FEA;
	MatrixXd Data = DATA;
	MatrixXd ERRWELD(3, wrow);
	NNSearchD* nns = NNSearchD::createKDTreeLinearHeap(Fea);
	MatrixXi indices(1, wrow);
	MatrixXi ewindices(1, wrow);
	MatrixXd dists2(1, wrow);
	nns->knn(Weld, indices, dists2, 1, 0, NNSearchF::SORT_RESULTS);
	Map<MatrixXd> NEWWELD(newweld, 3, wrow);
	double meanerr = dists2.row(0).cwiseSqrt().mean();
	__int64 j = 0,jj = 0;
	for (__int64 i = 0; i < wrow; ++i)
	{
		if (sqrt(dists2(0, i)) > meanerr)
		{
			NEWWELD.col(i) = Vector3d::Zero();
			ewindices(0,j)= i;
			ERRWELD.col(j)= Weld.col(i);
			++j;
		}
		else 
		{
			NEWWELD.col(i) = FEA.col(indices(0, i));			
		}
	}
	MatrixXd Errweld = ERRWELD.block(0, 0, 3, j);
	MatrixXi Ewindices = ewindices.block(0, 0, 1, j);
	NNSearchD* nnsD = NNSearchD::createKDTreeLinearHeap(Data);
	MatrixXi dataindices(1, j);
	MatrixXd datadists2(1, j);
	nnsD->knn(Errweld, dataindices, datadists2, 1, 0, NNSearchF::SORT_RESULTS);
	
	for (__int64 i = 0;i<j;++i)
	{
		NEWWELD.col(Ewindices(0,i))=DATA.col(dataindices(0,i));
	}
	delete nns, nnsD;
}
void icp_alg(double* model,__int64 mrow,double* target,__int64 trow,ICPOUT* result,double res,int maxiter,int miniter)
{
	//model means the static pointcloud
	//target is the moved pointcloud

	Matrix3d RT=Matrix3d::Identity();
	Vector3d TR=Vector3d::Zero();
	
	double* dataclose=(double*) calloc(trow*3,sizeof(double));
	Map<MatrixXd> DATACLOSE(dataclose,3,trow);
	//Map<MatrixXd> MODEL(model,3,mrow);
	Map<MatrixXd> TARGET(target,3,trow);
	Map<MatrixXd> DATAOUT(result->data_out,3,trow);
	//MatrixXd Dataclose=DATACLOSE.transpose();
	MatrixXd TARGETbuff=TARGET; //TARGET backup
	MatrixXd Target=TARGET.transpose();
	
	//for extrapolateInTransformSpace
	double* qs=(double*) calloc(maxiter*7,sizeof(double));
	double* dq=(double*) calloc(maxiter*7,sizeof(double));
	double* err=(double*) calloc(maxiter,sizeof(double));
	Map<MatrixXd> QS(qs,7,maxiter);
	Map<MatrixXd> DQ(dq,7,maxiter);

	for(int t=1;t<=maxiter;++t)
	{
		result->iter=t;
		clock_t start=clock();
		mulknn(target,model,trow,mrow,dataclose); 

		clock_t stop=clock();
		cout<<"once libnabo during time: "<<((double)(stop - start)) / CLOCKS_PER_SEC<<"s."<<"iter :"<<result->iter<<endl;
	    
		Vector3d t_mean= Target.colwise().mean().transpose();
		Vector3d d_mean= DATACLOSE.rowwise().mean();

		// caculate Covariance matrix
		//mat_anti=[(DATACLOSE-d_mean)*(Target-t_mean)]'
		Matrix3d mat_anti=((DATACLOSE-d_mean*VectorXd::Ones(trow).transpose())*(Target-VectorXd::Ones(trow)*t_mean.transpose())).transpose();
		//cout<<mat_anti<<endl<<endl;
		//svd
		JacobiSVD<MatrixXd> svd(mat_anti, ComputeFullU | ComputeFullV);

		//caculate RT=v*u' in each iteration
		RT=svd.matrixV()*svd.matrixU().transpose();
		if(RT.determinant()<0)
		{
			Matrix3d V=svd.matrixV();
			V.col(2)=-V.col(2);
			RT=V*svd.matrixU().transpose();
		}
	
		//caculate TR
		TR=d_mean-RT*t_mean;

		//caculate new target point cloud
		DATAOUT=Target.transpose();
 		Target=(RT*DATAOUT+TR*VectorXd::Ones(trow).transpose()).transpose();
		//update result
		result->RT=RT*result->RT;
		result->TR=TR+RT*result->TR;
		//caculate err
		err[t-1]=(Target- DATAOUT.transpose()).rowwise().norm().mean();
		result->res=err[t-1];
		
		//extrapolateInTransformSpace
		//double* resultRT=result->RT.data(); //save index by column
		//result->RT.transpose().data() is the same with result->RT.data(),so need to add a buff

		Matrix3d rtbuff;
		rtbuff=result->RT.transpose();
		rotationToQuaternion(rtbuff.data(),&qs[7*(t-1)]);
		QS.block<3,1>(4,t-1)=result->TR;
		if(t>1)
		{
			DQ.col(t-1)=QS.col(t-1)-QS.col(t-2);
			if(t>3)
			{
				extrapolateInTransformSpace(t,dq,qs,err);   //With extrapolation, we might be able to converge faster
			    // Update transformation and data

				quaternionToRotation(&qs[7*(t-1)],rtbuff.data());//updata result->RT
				result->RT=rtbuff.transpose();
				result->TR=QS.block<3,1>(4,t-1);//updata result->TR

				
				DATAOUT=result->RT*TARGETbuff+result->TR*VectorXd::Ones(trow).transpose();//updata result->dataout
				Target=DATAOUT.transpose();//update Target
			}
		}
		TARGET=Target.transpose();

		/*cout<<result->RT<<endl<<endl<<result->TR<<endl<<endl;
		  cout<<"qs: "<<qs[7*(t-1)]<<" "<<qs[7*(t-1)+1]<<" "<<qs[7*(t-1)+2]<<" "<<qs[7*(t-1)+3]<<" "<<qs[7*(t-1)+4]<<" "<<qs[7*(t-1)+5]<<" "<<qs[7*(t-1)+6]<<" "<<endl;
		*/
		stop=clock();
		cout<<"once ICP during time: "<<((double)(stop - start)) / CLOCKS_PER_SEC<<"s."<<"iter :"<<result->iter<<endl
		    <<"result->res: "<<endl<<result->res<<endl;
		if(t>2)
		{	
			if(abs(err[t-1]-err[t-2])<res && t>=miniter)
			{
				cout<<"Done"<<endl;
				mulknn(target, model, trow, mrow, dataclose);
				cout << "Final matching err: " << (Target - DATACLOSE.transpose()).rowwise().norm().mean() << " mm" << endl;
				break;
			}
		}

	}

	free(dataclose);free(err);free(qs);free(dq);
}
void rotationToQuaternion(double* RT,double* Q)
{
	double* R=(double*) calloc(3*3,sizeof(double));
	Map<MatrixXd> RTBUF(RT,3,3);
	Map<MatrixXd> RBUF(R,3,3);
	RBUF=RTBUF.transpose();

	double t=R[0]+R[4]+R[8];
	double r,s;
	
	if(t>=0)
	{
		r=sqrt(1+t);s=0.5/r;
		Q[0]=0.5*r;Q[1]=(R[5]-R[7])*s;Q[2]=(R[6]-R[2])*s;Q[3]=(R[1]-R[3])*s;
	}
	else
	{
		double maxv=max(R[0],max(R[4],R[8]));
		if(maxv==R[0])
		{
			r =sqrt(1+R[0]-R[4]-R[8]);s=0.5/r;
			Q[0]=(R[5]-R[7])*s;Q[1]=0.5*r;Q[2]=(R[1]+R[3])*s;Q[3]=(R[6]+R[2])*s;
		}
		else if(maxv==R[4])
		{
		    r =sqrt(1+R[4]-R[0]-R[8]);s=0.5/r;
			Q[0]=(R[6]-R[2])*s;Q[1]=(R[1]+R[3])*s;Q[2]=0.5*r;Q[3]=(R[5]+R[7])*s;
		}
		else
		{
		    r =sqrt(1+R[8]-R[0]-R[4]);s=0.5/r;
			Q[0]=(R[1]-R[3])*s;Q[1]=(R[6]+R[2])*s;Q[2]=(R[5]+R[7])*s;Q[3]=0.5*r;
		}
	}
	free(R);
}
void quaternionToRotation(double* Q,double* R)
{
	R[0]=Q[0]*Q[0]+Q[1]*Q[1]-Q[2]*Q[2]-Q[3]*Q[3];
	R[1]=2*(Q[1]*Q[2]-Q[0]*Q[3]);
	R[2]=2*(Q[1]*Q[3]+Q[0]*Q[2]);
	R[3]=2*(Q[1]*Q[2]+Q[0]*Q[3]);
	R[4]=Q[0]*Q[0]+Q[2]*Q[2]-Q[1]*Q[1]-Q[3]*Q[3];
	R[5]=2*(Q[2]*Q[3]-Q[0]*Q[1]);
	R[6]=2*(Q[1]*Q[3]-Q[0]*Q[2]);
	R[7]=2*(Q[3]*Q[2]+Q[0]*Q[1]);
	R[8]=Q[0]*Q[0]+Q[3]*Q[3]-Q[1]*Q[1]-Q[2]*Q[2];
}
void extrapolateInTransformSpace(int t,double* dq,double* qs,double* err)
{
//Extrapolation in quaternion space. Details are found in:
//Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
// IEEE Transactions on pattern analysis and machine intelligence, 239-256.
	double dtheta1,dtheta2,n1,n2,n3;
	double angletheshold=10,scalefactor=25;
	Map<MatrixXd> QS(qs,7,Dynamic);
	Map<MatrixXd> DQ(dq,7,Dynamic);
	n1=DQ.col(t-3).norm();
	n2=DQ.col(t-2).norm();
	n3=DQ.col(t-1).norm();
	dtheta1=(180/pi)*acos(DQ.col(t-2).dot(DQ.col(t-3))/n1/n2);
	dtheta2=(180/pi)*acos(DQ.col(t-1).dot(DQ.col(t-2))/n2/n3);
	if(dtheta1<angletheshold && dtheta2<angletheshold )
	{
		double d[3]={err[t-1],err[t-2],err[t-3]};
		double v[3]={0,-n3,-n2-n3};
		double vmax=scalefactor*n3;
		double dv=extrapolate(v,d,vmax);
		if(dv!=0)
		{
			QS.col(t-1)=QS.col(t-1)+DQ.col(t-1)*dv/n3;
			QS.block<4,1>(0,t-1)=QS.block<4,1>(0,t-1)/QS.block<4,1>(0,t-1).norm();
		}
	}
}
double extrapolate(double* v,double* d,double vmax)
//Extrapolation in quaternion space. Details are found in:
//Besl, P., & McKay, N. (1992). A method for registration of 3-D shapes. 
// IEEE Transactions on pattern analysis and machine intelligence, 239-256.
{
	double dv=0,v1,v2;
	double* err=(double*) calloc(3,sizeof(double));
	double sco1[3],sco2[3],list[4]={0,0,0,vmax};
	datafit(v,d,3,sco1,2,err); //linear fit
    datafit(v,d,3,sco2,3,err); //parabolic fit
	v1=-sco1[0]/sco1[1];
	v2=-sco2[1]/sco2[2]/2;
	list[1]=v1;list[2]=v2;
	sort(list,list+4);
	if(list[0]==0)
	{
		dv=list[1];return dv;
	}
	else if(list[1]==0)
	{
	    dv=list[2];return dv;
	}
	else if(list[2]==0)
		dv=list[3];
	free(err);
	return dv;
}
void datafit(double *x,double *y,int n,double *a,int m,double *dt)
//  int n,m;   m=num(a);
//  double x[],y[],a[],dt[];
//y=a[0]+a[1]*(x-avg(x))+a[2]*(x-avg(x))^2+.....   
  { int i,j,k;
    double z,p,c,g,q,d1,d2,s[20],t[20],b[20];
    for (i=0; i<=m-1; i++) a[i]=0.0;
    if (m>n) m=n;
    if (m>20) m=20;
    z=0.0;
    for (i=0; i<=n-1; i++) z=z+x[i]/(1.0*n);
    b[0]=1.0; d1=1.0*n; p=0.0; c=0.0;
    for (i=0; i<=n-1; i++)
      { p=p+(x[i]-z); c=c+y[i];}
    c=c/d1; p=p/d1;
    a[0]=c*b[0];
    if (m>1)
      { t[1]=1.0; t[0]=-p;
        d2=0.0; c=0.0; g=0.0;
        for (i=0; i<=n-1; i++)
          { q=x[i]-z-p; d2=d2+q*q;
            c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
          }
        c=c/d2; p=g/d2; q=d2/d1;
        d1=d2;
        a[1]=c*t[1]; a[0]=c*t[0]+a[0];
      }
    for (j=2; j<=m-1; j++)
      { s[j]=t[j-1];
        s[j-1]=-p*t[j-1]+t[j-2];
        if (j>=3)
          for (k=j-2; k>=1; k--)
            s[k]=-p*t[k]+t[k-1]-q*b[k];
        s[0]=-p*t[0]-q*b[0];
        d2=0.0; c=0.0; g=0.0;
        for (i=0; i<=n-1; i++)
          { q=s[j];
            for (k=j-1; k>=0; k--)
              q=q*(x[i]-z)+s[k];
            d2=d2+q*q; c=c+y[i]*q;
            g=g+(x[i]-z)*q*q;
          }
        c=c/d2; p=g/d2; q=d2/d1;
        d1=d2;
        a[j]=c*s[j]; t[j]=s[j];
        for (k=j-1; k>=0; k--)
          { a[k]=c*s[k]+a[k];
            b[k]=t[k]; t[k]=s[k];
          }
      }
    dt[0]=0.0; dt[1]=0.0; dt[2]=0.0;
    for (i=0; i<=n-1; i++)
      { q=a[m-1];
        for (k=m-2; k>=0; k--)
          q=a[k]+q*(x[i]-z);
        p=q-y[i];
        if (fabs(p)>dt[2]) dt[2]=fabs(p);
        dt[0]=dt[0]+p*p;
        dt[1]=dt[1]+fabs(p);
      }
	//change scoeff
	double avg=(x[0]+x[1]+x[2])/3;
	if(m==3)  
	{
	    a[0]=a[0]-a[1]*avg+a[2]*avg*avg;
		a[1]=a[1]-2*a[2]*avg;
	}
  	if(m==2)
	{
		a[0]=a[0]-a[1]*avg;
	}
    return;
  }
int Readparam(char* fileName, float* trig, int* knum, float* anglethrold, float* scal)
{
	FILE *file;
	string name;
	if (fopen_s(&file, fileName, "r"))
	{
		printf("Can not open the mparam!\n");
		return -1;
	}
	else//文件存在  
	{
		int LINES = CountLines(fileName);
		for (int i = 0; i<LINES; i++)
		{
			fscanf_s(file, "%s", &name);
			if (!strcmp(name.data(), "trig")) { fscanf_s(file, "%f", trig); continue; }
			if (!strcmp(name.data(), "knum")) { fscanf_s(file, "%d", knum); continue; }
			if (!strcmp(name.data(), "anglethrold")) { fscanf_s(file, "%f", anglethrold); continue; }
			if (!strcmp(name.data(), "scale")) { fscanf_s(file, "%f", scal); continue; }
			if (!strcmp(name.data(), "datasavepath")) { fscanf_s(file, "%s", &name); continue; }

		}

		fclose(file); //关闭文件  
		return 0;
	}

}