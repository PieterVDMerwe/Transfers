#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <io.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <omp.h>
#include <sstream>  // for std::ostringstream
#include <string>   // for std::string
#include <algorithm>
#include <iterator>
#include <random>
#include <iomanip>
#include <map>
#include <bits/stdc++.h>
#include <chrono>

using namespace std;

const double G = 6.67259e-11;//(4.0*pow(M_PI,2))/m;//
const double me = 9.10938356e-31;//kg
const double ce = 1.602e-19;//coulomb
const double pie = 3.14159;
const double e0 = 8.85e-12;//c2/Nm2
const double M_sun = 1.989e30;
const double c = 2.997e8;
const double R = 2.0*G*M_sun/pow(c,2);
const double y_offset = 30e6;
const double x_offset = 30e6;//30e6;//160.0e3*R;
const double larmorcoeff = pow(ce,2)/(6.0*pie*e0*pow(c,3));
//const double y_offset = 10e7;
//const double x_offset = 10e7;
const int n_arr = 13;//29;//13 for best
const long int i_lim =  pow(2,n_arr);
const int width = 160;
const int height = 160;
const double seperation = 10e6;

// fftw_complex TotRadLos[i_lim];
// fftw_complex FFTResult[i_lim];
// fftw_complex *TotRadLos = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * i_lim);
// fftw_complex *FFTResult = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * i_lim);

double orbitCoeff[3][2];
double **pulseProp = NULL;




ofstream logfile;


double getRandomNumber(double upper,  double lower)
{
    std::random_device rd;
    std::mt19937 e2(rd());
    //std::mt19937 e2(25);
    std::uniform_real_distribution<double> dist(lower,  upper);
    return dist(e2);
}



void ReadFromFileOrbit()
{
    ifstream file("./fit.csv");
    for(int i = 0;i < 2;i++)
    {
        string line;
        getline(file, line);
        if(i>-1)
        {
            stringstream iss(line);
            for(int j = 0; j<3;j++)
            {
                string val;
                getline(iss,val,',');
                stringstream convertor(val);
                convertor >> orbitCoeff[j][i];
            }
        }
    }
}

void Orbit(double M1,double M2,double dt)
{
    ofstream myfile,myfile2;
    myfile.open("Particles/OrbitOutput.txt");
    myfile2.open("Particles/OrbitOutput2.txt");
    double dx = c*dt;
    double R1_0 = 40.0e6;//1e7;//140.0e3*R;
    double R2_0 = R1_0-seperation;//60.0e3*R;


    double* x1 = new double[i_lim];
    double* x2 = new double[i_lim];
    double* y1 = new double[i_lim];
    double* y2 = new double[i_lim];
    double* x = new double[i_lim];
    double* y = new double[i_lim];
    double* vx = new double[i_lim];
    double* vy = new double[i_lim];
    double* v_rel = new double[i_lim];
//    double a_sm[i_lim];
//    double e[i_lim];

    x1[0] = R1_0;
    x2[0] = R2_0;
    y1[0] = 0.0;
    y2[0] = 0.0;
    x[0] = x1[0]-x2[0];
    y[0] = y1[0]-y2[0];
    double m = M1*M2/(M1+M2);
    double G = 6.67259e-11;//(4.0*pow(M_PI,2))/m;//
    double Rs = 2.0*G*m/pow(c,2);
    vx[0] = 0.0;
    vy[0] = sqrt((G*m)/(R1_0-R2_0));// -0.5e6 ;//175e4;
    double* r = new double[i_lim];

    cout << "R: " << R << "\n";
    cout << "R1: " << R1_0 << "\n";
    cout << "R2: " << R2_0 << "\n";

    myfile2 << "x"<<','<<"y"<<'\n';
    myfile << "x1"<<','<<"y1"<<','<<"x2"<<','<<"y2" <<','<<"t" << '\n';
    for(int i=0;i<(i_lim-1);i=i+1)//1;i=i+1)//
    {

        r[i] = sqrt(pow(x[i],2) + pow(y[i],2));

        vx[i+1] = vx[i]-(G*m*x[i]*dt)/(r[i]*pow((r[i]-Rs),2));
        vy[i+1] = vy[i]-(G*m*y[i]*dt)/(r[i]*pow((r[i]-Rs),2));
        x[i+1] = x[i]+vx[i+1]*dt;
        y[i+1] = y[i]+vy[i+1]*dt;

        x1[i+1] = (M2*x[i+1])/(M1+M2)+x_offset;
        x2[i+1] = -(M1*x[i+1])/(M1+M2)+x_offset;
        y1[i+1] = (M2*y[i+1])/(M1+M2)+y_offset;
        y2[i+1] = -(M1*y[i+1])/(M1+M2)+y_offset;

        myfile2 << x[i+1]<<','<<y[i+1]<<'\n';
        myfile << x1[i+1]<<','<<y1[i+1]<<','<<x2[i+1]<<','<<y2[i+1] <<','<<static_cast<double>(i)*dt << '\n';
        //cout <<  x1[i+1]<<','<<y1[i+1]<<','<<x2[i+1]<<','<<y2[i+1] <<','<<static_cast<double>(i)*dt << '\n';
    }

    myfile2.close();
    myfile.close();
//
    delete[] x1;
    delete[] x2;
    delete[] y1;
    delete[] y2;
    delete[] x ;
    delete[] y ;
    delete[] vx;
    delete[] vy;
    delete[] v_rel;
    delete[] r;
}

class TestParticle
{
    public:


//     double ***localPotential = NULL;
    double curPos[3];
    double curVeloc[3];
    double curAcel[3];
    double prevPos[3];
    double prevVeloc[3];
    double prevAcel[3];
    double prevMomentum[3];
    double curMomentum[3];
    double time0;
//     double **orbitArr = NULL;
    double dt;
    double dx;
    double M1, M2;
    double Orbit1Coeff[3];
    double Orbit2Coeff[3];
    double RadLoss;
    double aperp;
    double apar;
    double velocDirection[3];
    double aPerpVec[3];
    int particleNum;
//     long double* TotalRadLos = new long double[i_lim];
    ofstream particleFile;
    double gamgam;
    

// double **tempOrbitArr,
    void initialize(double pos0[3], double vel0[3], double acel0[3], double time0temp,double dtTemp, double M1Temp, double M2Temp, int tempParticleNum)
    {
        particleNum=tempParticleNum;
        M1 = M1Temp;
        M2 = M2Temp;
        dt = dtTemp;
        time0 = time0temp;
        for(int q= 0;q<3;q++)
        {
            curPos[q] = pos0[q];
            curVeloc[q] = vel0[q];
            curAcel[q] = acel0[q];
        }
        curMomentum[0] = calcMomentumElement(curVeloc[0],curVeloc[1],curVeloc[2],curVeloc[0]);
        curMomentum[1] = calcMomentumElement(curVeloc[0],curVeloc[1],curVeloc[2],curVeloc[1]);
        curMomentum[2] = calcMomentumElement(curVeloc[0],curVeloc[1],curVeloc[2],curVeloc[2]);
        for(int m = 0;m<3;m++)
        {
            Orbit1Coeff[m] = orbitCoeff[m][0];
            Orbit2Coeff[m] = orbitCoeff[m][1];
        }
        dx = c*dt;

    }


    bool calcAccel(double t, double tprime1,  double tprime2, double dx)
    {
        double x1 = (curPos[0])-xOrbit(tprime1,Orbit1Coeff[0],Orbit1Coeff[1],Orbit1Coeff[2]);
        double y1 = (curPos[1])-yOrbit(tprime1,Orbit1Coeff[0],Orbit1Coeff[1],Orbit1Coeff[2]);
        double x2 = (curPos[0])-xOrbit(tprime2,Orbit1Coeff[0],Orbit2Coeff[1],Orbit2Coeff[2]);
        double y2 = (curPos[1])-yOrbit(tprime2,Orbit1Coeff[0],Orbit2Coeff[1],Orbit2Coeff[2]);
        double  z = curPos[2];
        double r1 = sqrt(pow(x1,2)+pow(y1,2)+pow(z,2));
        double r2 = sqrt(pow(x2,2)+pow(y2,2)+pow(z,2));

        double x1True = (curPos[0])-xOrbit(t,Orbit1Coeff[0],Orbit1Coeff[1],Orbit1Coeff[2]);
        double y1True = (curPos[1])-yOrbit(t,Orbit1Coeff[0],Orbit1Coeff[1],Orbit1Coeff[2]);
        double x2True = (curPos[0])-xOrbit(t,Orbit1Coeff[0],Orbit2Coeff[1],Orbit2Coeff[2]);
        double y2True = (curPos[1])-yOrbit(t,Orbit1Coeff[0],Orbit2Coeff[1],Orbit2Coeff[2]);
        double r1True = sqrt(pow(x1True,2)+pow(y1True,2)+pow(z,2));
        double r2True = sqrt(pow(x2True,2)+pow(y2True,2)+pow(z,2));


//         cout << '\n' << "x1True: " << x1True << "          r1True: " << r1True;

        double r_s1 = 2.0*G*M1/pow(c,2);
        double r_s2 = 2.0*G*M2/pow(c,2);

        if ((r1True<r_s1) || (r2True<r_s2))
        {
            return true;;
        }
        if ((r1True>9.0*seperation) && (r2True>9.0*seperation))
        {
            return true;;
        }

/*        curAcel[0] = acceli(r1, r_s1, x1, M1) + acceli(r2, r_s2, x2, M2);
        curAcel[1] = acceli(r1, r_s1, y1, M1) + acceli(r2, r_s2, y2, M2);
        curAcel[2] = acceli(r1, r_s1, z , M1) + acceli(r2, r_s2, z , M2);  */

//         curAcel[0] = acceli(r1True, r_s1, x1True, M1, 0) + acceli(r2True, r_s2, x2True, M2, 0);
//         curAcel[1] = acceli(r1True, r_s1, y1True, M1, 1) + acceli(r2True, r_s2, y2True, M2, 1);
//         curAcel[2] = acceli(r1True, r_s1, z , M1, 2)     + acceli(r2True, r_s2, z , M2, 2);
        for(int i =0;i<3;i++)
        {
            prevMomentum[i]=curMomentum[i];
            prevVeloc[i] = curVeloc[i];
        }

        double RargTemp = RArg(r1, r_s1, x1, M1, r2, r_s2, x2, M2, prevMomentum[0]);
        curMomentum[0] = RargTemp;
        RargTemp = RArg(r1, r_s1, y1, M1, r2, r_s2, y2, M2, prevMomentum[1]);
        curMomentum[1] = RargTemp;
        RargTemp = RArg(r1, r_s1, z, M1, r2, r_s2, z, M2, prevMomentum[2]);
        curMomentum[2] = RargTemp;

        curVeloc[0] = calcVelocityElement(curMomentum[0],curMomentum[1],curMomentum[2],curMomentum[0]);
        curVeloc[1] = calcVelocityElement(curMomentum[0],curMomentum[1],curMomentum[2],curMomentum[1]);
        curVeloc[2] = calcVelocityElement(curMomentum[0],curMomentum[1],curMomentum[2],curMomentum[2]);
        
        double velocNorm = sqrt(curVeloc[0]*curVeloc[0]+curVeloc[1]*curVeloc[1]+curVeloc[2]*curVeloc[2]);     
        velocDirection[0] = curVeloc[0]/velocNorm;
        velocDirection[1] = curVeloc[1]/velocNorm;
        velocDirection[2] = curVeloc[2]/velocNorm;

        for (int l = 0;l<3;l++ )
        {
            curAcel[l] = deriv(prevVeloc[l],curVeloc[l],dt);
        }
//         gamgam = (gamma(curVeloc[0],curVeloc[1],curVeloc[2]));
        gamgam = curMomentum[0]/(me*curVeloc[0]);
        apar = velocDirection[0]*curAcel[0]+velocDirection[1]*curAcel[1]+velocDirection[2]*curAcel[2];
        
        aPerpVec[0] = curAcel[0]-apar*velocDirection[0];
        aPerpVec[1] = curAcel[1]-apar*velocDirection[1];
        aPerpVec[2] = curAcel[2]-apar*velocDirection[2];
        
        aperp = sqrt(aPerpVec[0]*aPerpVec[0]+aPerpVec[1]*aPerpVec[1]+aPerpVec[2]*aPerpVec[2]);
//         RadLoss  = larmorcoeff*pow(gamgam,4)*((pow(curAcel[0],2)+pow(curAcel[1],2)+pow(curAcel[2],2))+(pow(gamgam,2))*pow((curAcel[0]*curVeloc[0]+curAcel[1]*curVeloc[1]+curAcel[2]*curVeloc[2])/c,2));
        if(RadLoss<0)
        {
            cout << "\n Negative radloss \n";
        }

//         int iTemp = t/dt;
//         TotalRadLos[iTemp] = RadLoss;






        return false;
    }

    double deriv(double vprev, double vcur, double dt)
    {
        return (vcur-vprev)/dt;
    }

    void calcPos()
    {
        for(int i=0;i<3;i++)
        {
            prevPos[i]=curPos[i];
            curPos[i]=prevPos[i]+curVeloc[i]*dt;
        }
    }

    double acceli(double r,  double rs, double xi, double Mi)
    {

        double GM = G*Mi;
        double rrs = pow((r-rs), -2.0);
        double xr = xi/r;

        return -GM*rrs*xr;
    }

    double calcMomentumElement(double veloc1, double veloc2, double veloc3, double velocElement)
    {
        return gamma(veloc1,veloc2,veloc3)*velocElement*me;
    }

    double calcVelocityElement(double momentum1,double momentum2, double momentum3, double momentumElement)
    {
        double pSqrd = pow(momentum1, 2)+pow(momentum2, 2)+pow(momentum3, 2);
        double coeff = pSqrd/pow((c*me) , 2);
        double tempRet = (momentumElement/me)*pow( (1.0-coeff*pow( (1.0+ coeff)  ,  -1)) , 0.5 );//*pow( 1.0/coeff , 0.5 );//*pow( (1.0-coeff*pow( (1.0+ coeff)  ,  -1)) , 0.5 );//
        return  tempRet;
    }

    double forceElement(double r,  double rs, double xi, double Mi)
    {
        double GM = G*Mi;
        double rrs = pow((r-rs), -2.0);
        double xr = xi/r;
        return -GM*rrs*xr;
    }

    double gamma(double vx, double vy, double vz)
    {
        double vsqr = pow(vx,2.0)+pow(vy,2.0)+pow(vz,2.0);
        double beta = vsqr/pow(c,2.0);
        double arg = 1.0-beta;
        return 1.0/pow(arg,0.5);
    }

    double RArg(double r1,  double rs1, double xi1, double M1,double r2,  double rs2, double xi2, double M2, double MomEl)
    {
        double temp = (forceElement(r1,rs1,xi1,M1)+forceElement(r2,rs2,xi2,M2))*dt*me+MomEl;
//         #pragma omp critcal
//         {
//             cout << temp  <<"  \n "                  ;
//         }
        return temp;
    }


//     void deleteOrbitArr()
//     {
//         for (int i=0;i<5;i++)
//         {
//             delete[] orbitArr[i];
//         }
//         delete[] orbitArr;
//     }

//     void createOrbitArr()
//     {
//         orbitArr = new double *[4];
//         for (int i=0;i<5;i++)
//         {
//             orbitArr[i] = new double [i_lim];
//         }
//     }

    double newtonCirc(double t,double a,double w,double x0,double xp,double y0,double yp,double zp,double tprime)
    {
       return  pow(c*(t-tprime),2) -pow(zp,2)-pow((xp-x0-a*cos(w*tprime)),2)-pow((yp-y0-a*sin(w*tprime)),2);
    }

    double derivNewtonCirc(double t,double a,double w,double x0,double xp,double y0,double yp,double zp,double tprime)
    {
       return  -2.0*pow(c,2)*(t-tprime)-2.0*(xp-x0)*a*w*sin(w*tprime)+2.0*(yp-y0)*a*w*cos(w*tprime);//-2.0*pow(a,2)*w*cos(w*tprime)*sin(w*tprime)+2.0*pow(a,2)*w*cos(w*tprime)*sin(w*tprime);
    }

    double retartedTime(double t,double a,double w,double x0,double xp,double y0,double yp,double zp)
    {
        double tprime0 = 0;
        double tprimediff = 1;
        while (tprimediff > 0.000000000001)
        {
            double tprime0old = tprime0;
            tprime0 = tprime0old-newtonCirc(t, a,w,x0,xp,y0,yp,zp,tprime0old)/derivNewtonCirc(t, a,w,x0,xp,y0,yp,zp,tprime0old);
            tprimediff = abs(tprime0-tprime0old);

        }

        return tprime0;

    }

    double xOrbit(double t,  double a,  double w, double x0)
    {
        return a*cos(w*t)+x0;

    }

    double yOrbit(double t,  double a,  double w, double y0)
    {
        return a*sin(w*t)+y0;

    }


    void particleTrajectory()
    {
        string str = std::to_string(particleNum);
        particleFile.open("Particles/particle_"+str+"_acceleration.txt");
//         particleFile << "Ax,Ay,Az,Vx,Vy,Vz,x,y,z,t,radlos,px,py,pz,gamma\n";
        particleFile << "t,gamma,apar,aperp\n";

        bool particlePlummet = false;
        for(int i = 0;i<i_lim;i++)
        {
//             particleFile << curAcel[0]<<","<<curAcel[1]<<","<<curAcel[2]<<","<< curVeloc[0]<<","<<curVeloc[1]<<","<<curVeloc[2]<<","<< curPos[0]<<","<<curPos[1]<<","<<curPos[2]<<","<<static_cast<double>(i)*dt<<","<<RadLoss<< " , " << curMomentum[0]<< " , "<< curMomentum[1]<< " , "<< curMomentum[2]<<" , " << gamgam << "\n";
            particleFile << static_cast<double>(i)*dt<<","<<gamgam<<","<< apar <<"," << aperp << "\n";

            for(int j=0;j<3;j++)
            {
                prevAcel[j]=curAcel[j];
            }
            double tprime1 = retartedTime(static_cast<double>(i)*dt,Orbit1Coeff[0],Orbit1Coeff[1],Orbit1Coeff[2],curPos[0],Orbit1Coeff[2],curPos[1],curPos[2]);
            double tprime2 = retartedTime(static_cast<double>(i)*dt,Orbit2Coeff[0],Orbit2Coeff[1],Orbit2Coeff[2],curPos[0],Orbit2Coeff[2],curPos[1],curPos[2]);
            particlePlummet = calcAccel(static_cast<double>(i)*dt, tprime1,  tprime2 ,dx);


            if (particlePlummet == true)
            {
                i = i_lim;
            }
            else
            {
                calcPos();
            }

        }

        particleFile.close();
    }

//     void killParticle()
//     {
//         delete[] TotalRadLos;
//     }

//     void destruct()
//     {
//         deleteOrbitArr();
//         deleteLocalPotential();
//     }
// 
// 
//     void createLocalPotential()
//     {
//             localPotential = new double **[3];
//             for (int i=0;i<3;i++)
//             {
//                 localPotential[i] = new double *[3];
//                 for (int j=0;j<3;j++)
//                 {
//                     localPotential[i][j] = new double [3];
// 
//                 }
//             }
//     }

//     void deleteLocalPotential()
//     {
//             for (int i=0;i<3;i++)
//             {
//                 for (int j=0;j<3;j++)
//                 {
// 
//                     delete[] localPotential[i][j];
//                 }
//                 delete[] localPotential[i];
//             }
//             delete[] localPotential;
//     }

};

int main()
{

    cout << i_lim << '\n';
    cout << "Starting..\n";
/////////////////////////////////////////////////////////////////////New code/////////////////////////////////////////////////////////////////////
    auto start = chrono::high_resolution_clock::now();

    // unsync the I/O of C and C++.
   // ios_base::sync_with_stdio(false);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//     for(int kk = 0; kk<i_lim;kk++)
//     {
//         TotalRadLos[kk] = 0.0;
//     }


            logfile.open("./log.txt");

    cout << "Time taken\n";
    double M1, M2;
    double dt = 1.5e1/c;
    double dx = c*dt;
    M1 = 30.0*M_sun;
    M2 = 35.0*M_sun;
//     pulseProp = new double *[4];
//     for (int i=0;i<5;i++)
//     {
//         pulseProp[i] = new double [i_lim];
//     }
    //Orbit(M1,M2,dt);


    //cout << "Starting python prog\n";
    //system("python ./curveFit.py");
    //cout << "python prog done\n";

    ReadFromFileOrbit();
    


    cout <<  "Starting Particle Trajectories...\n";
    int numParticles = 50;
    ofstream parmFile;
    parmFile.open("Parameter_file.csv");
    parmFile << "i_lim,dt,numPart\n";
    parmFile << i_lim << "," << dt << "," << numParticles << "\n";
    parmFile.close();
    int particlesDone = 0;
    cout << "Jip\n";
    TestParticle part1[numParticles];
    cout << "Nope\n";
    #pragma omp parallel for
    for (int i = 0; i < numParticles; i++)
    {
        
        
        double dummi1[3];
        double dummi2[3];
        double dummi3[3];

        for (int l = 0;l<3;l++)
        {
            if(l<2)
            {
                dummi1[l] = getRandomNumber(1.0e7, 4.5e7);
            }
            else
            {
                dummi1[l] = getRandomNumber(-4.5e7, 4.5e7);
            }
            dummi2[l] = getRandomNumber(-5000.0, 5000.0);
            dummi3[l] = 0.0;
        }
//         dummi1[0] = 4.31272e+07;
//         dummi1[1] = 2.01421e+07;
//         dummi1[2] =-2.45448e+06;
//         dummi2[0] =-4163.85;
//         dummi2[1] = 3264.69;
//         dummi2[2] = 3825.17;
//pulseProp,
        part1[i].initialize(dummi1,dummi2,dummi3,0.0,dt,M1,M2,i);
        part1[i].particleTrajectory();
//         string str1 = std::to_string(dt);
//         string str2 =std::to_string(i_lim);
//         string str3 = std::to_string(i);
//         system(("python ./fft.py "+str1+" "+str2+" "+str3).c_str());
        #pragma omp critical
        {
//             for(int b; b<i_lim; b++)
//             {
//                 TotRadLos[b][REAL] = TotRadLos[b][REAL] + part1[i].TotalRadLos[b];
//             }
            particlesDone++;
            cout <<  "Particles done----------------------------------" << particlesDone << "/" << numParticles << "\n";
        }
//         part1[i].killParticle();
// 
    }

 ///////////////////////////////FFT-stuff////////////////////////////////////////////
/*
    cout << "whoop whoop 0\n";
    for(int k = 0;k<i_lim;k++)
    {
        TotRadLos[k][IMAG] = 0.0;
    }
    cout << "whoop whoop 1\n";
    
    fftw_plan plan = fftw_plan_dft_1d(i_lim, TotRadLos, FFTResult, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_cleanup();
    
    
    ofstream FFTResultFile;
    cout << "whoop whoop 2\n";
    
    FFTResultFile.open("FFTResults.txt");
    FFTResultFile <<"Re , Im\n";
    for(int e = 0; e<i_lim; e++)
    {
        FFTResultFile << FFTResult[e][REAL] << " , " << FFTResult[e][IMAG] << "\n";
    }
    FFTResultFile.close();*/

 ///////////////////////////////////////////////////////////////////////////////////

//     for (int i = 0; i < numParticles; i++)
//     {
//         part1[i].killParticle();
//     }


//     for (int i=0;i<5;i++)
//     {
//         delete[] pulseProp[i];
//     }
//     delete[] pulseProp;



//     ofstream myfileRadlos;
//     myfileRadlos.open("Particles/Radloss.txt");
//     for(int i = 0; i<i_lim;i++)
//     {
//         myfileRadlos << TotalRadLos[i] << '\n' ;
//     }
//     myfileRadlos.close();


    


    logfile.close();
/////////////////////////////////////////////////////////////////////New code/////////////////////////////////////////////////////////////////////
    auto end = chrono::high_resolution_clock::now();

    // Calculating total time taken by the program.
    double time_taken =
      chrono::duration_cast<chrono::nanoseconds>(end - start).count();

    time_taken *= 1e-9;

    cout << "Time taken by program is : " << fixed
         << time_taken << setprecision(9);
    cout << " sec" << endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    return 0;
}
