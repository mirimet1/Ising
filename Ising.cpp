/*****************************************************************************/
/*** Ising_a                                                               ***/
/*** Mir Mohammad Ebrahimi at the (IASBS)                                               ***/
/*** Date: 13970501                                                        ***/
/*** Ising code Dr. Niri has used                                          ***/
/*****************************************************************************/

/* You can compile this code by the following command in shell:
   > g++ -o Ising.out Ising.cpp -Ofast
   , then execute it in shell by
   > ./Ising.out 
   and take the results. */

/* This is a sample code. It is not complete and it may contain some errors. 
   So, you can try it with your own risk! */

/* These code is based on the random number generator function

   double ran2(long *idum);
   
which is implemented in the Numerical recipes in C, chapter 7 (ran2)

``Long period (> 2 × 10^{18}) random number generator of L. Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). 

***!!! Call with idum a negative integer to initialize; !!!*** thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1."

Visit www.nr.com for the licence.*/

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>

using namespace std;

//BEGIN_FOLD - Random number generator section
//------------------------------------------------------------------------------

/* following routine is based on the random number generator

   double ran2(long *idum);
   
which is implemented in the Numerical recipes in C, chapter 7 (ran2)

``Long period (> 2 × 10^{18}) random number generator of L. Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). 

***!!! Call with idum a negative integer to initialize; !!!*** thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1."

Visit www.nr.com for the licence.*/

// This is a internal, 32 bit random number generator with uniform Distribution in range [0..1)
/* note #undef's at end of file */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum) {
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

//------------------------------------------------------------------------------

// My interface for the Numerical recipes PRNG routine

// Random number generator seed
long iseed = -36;

// You can initializing the random number generator seed by using the following routine.
void Randomize() {
  iseed = -time(NULL);  
}

// Random() returns a random number with uniform distribution in the range of [0..1), 
// where 1 is excluded.
inline double Random() { return ran2(&iseed); }

// Random(int N) returns an integer random number with uniform distribution between 0 
// to N-1 (both boundary are included).
inline int Random(int N) { return int(ran2(&iseed)*N); }

//------------------------------------------------------------------------------
//END_FOLD   - Random number generator section

const int L = 256;                // Lattice size; L = 2^n
int** s;                          // Spins lattice; s[i][j] shows the spin of cell at (i,j)

// All the following parameters are in the reduced unit system
float h = 0.0;                    // external field
float J = +1;                     // spin-spin coupling constant: +1 for ferro and -1 for antiferro
float T = 1.5;                    // tempreture
float m = 0.0;                    // magnetization

ofstream snapshot("snapshot.txt", ios::out | ios::trunc);            // An stream Keeps the model snapshots
ofstream snapshot1("magnetization.txt", ios::out | ios::trunc);          
                 

void Init();                      // Init() initializes the model, alocates the dynamic memory, ...
void Execute();                   // Execute() simulates the Ising model
void PostProcess();               // PostProcess() performs the post statistical analysis
void Done();                      // Done() deallocates  the dynamic memory, ...

// Initialize the model, alocate dynamic memory, ...
void Init() {
    cout << "Init ..." << endl;
    
    Randomize();                  // iseed is randomized
    
    // In the following lines the memory is allocated for s
    s = new int*[L];
    for (int i = 0; i<L; i++) {
      s[i] = new int[L];
      for (int j = 0; j<L; j++)
        s[i][j] = 1;
    }
    
    snapshot << L << endl;        // The lattice size is written to the snapshots' stream
}

// Deallocate dynamic memory, ...
void Done() {
    cout << "Done ..." << endl;
 
    snapshot.close();             // The snapshots stream is closed

    // In the following lines the memory of s is deallocated
    for (int i = 0; i<L; i++)
      delete[] s[i];
    delete[] s;
}

// Export single snapshot to the snapshots' stream
void Export_Single_Snapshot() {
    for (int i = 0; i<L; i++){
	
	
	 
      for (int j = 0; j<L; j++){
	  
         
      //  if ( s[i][j] > 0 )               	      
          if ( Random(L) % 2  )        
		  
		 
          
             snapshot << "1" << '\t';
         else
             snapshot<< "0" << '\t';
        }
    
         snapshot<<'\n';
   }  
} 

// Period(int i) rotate i periodically if it is cross the boundary
inline int Period(int i) {
   // const int LMask = L - 1;      // Mask is used to perfom periodic boundary condition
     return (i+L) % L;
   // return i & LMask;
     //or equivalently use: return (i+L) % L;
} 

// Calculate Delta E if s(i,j) flipped
float Delta_E_Flip(int i, int j) {
    float sum = 0.0;
    sum = s[i][Period(j+1)] + s[i][Period(j-1)] + s[Period(i-1)][j] + s[Period(i+1)][j];
    return 2*s[i][j]* (h + J * sum);
}
float magnet()
{
for (int i = 0; i < L; i++)	
for (int j = 0; j < L; j++)	
      m += s[i][j];
return m/(L*L);
}

// Perfom single Monte Carlo step
void Single_Monte_Carlo_Step() {
    // Loop on L^2 infinitesimate steps
    for (int c = 0; c < L*L; c++) {
        int i = Random(L);        // A random row between [0..L-1] is selected
        int j = Random(L);        // A random column between [0..L-1] is selected
        float dE = Delta_E_Flip(i, j);
        if ( dE < 0 ) {           // If the total energy decreases after spin-flip, accept the infinitesimal step
            s[i][j] = -s[i][j];
        } else if ( Random() < exp(-dE/T) ) {
            s[i][j] = -s[i][j];
        }
  
   
    }
     Export_Single_Snapshot();         
     
   
}
   

// Simulate Ising model
void Execute() {
    cout << "Execute the Metropolis algorithm ..." << endl;

    
   for (int i = 0; i < 500; i++ ){
   	
      Single_Monte_Carlo_Step();
      snapshot1<<magnet()<<endl;      
        m = 0;                           
   }
    
 
}
// Perform post statistical analysis
void PostProcess(){
    cout << "Perfom the post statistical analysis ..." << endl;
    // You should perfom the post statistical analysis here
    // .
    // .
    // .

}
// Main routine
int main (int argc, char *argv[]) {
    Init();   

    Execute();
    Done();
   // Zakhire va plot meghnatesh bar hasabe ghadamhaye monte_carlo 
   ofstream plot("plot.txt",ios::out | ios::trunc); 
                    

     plot << "set ylabel \"Magnetization\"\n"
               <<"set yrange [-1:1]\n"
               << "set xlabel \"MonteCarlo Step\"\n" 
               << "set autoscale" << endl;

    plot<<"plot(\"magnetization.txt\") with lp \n";
    plot <<"pause -1";
    plot.close();
    system("gnuplot plot.txt");
    
    cout << "\nFinish" << endl;
    
    return 0;
}
