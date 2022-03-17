#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <time.h>
#include <cstdlib>
#include <iomanip>

using namespace std;


/*
 * PSUEDOCODE FOR THE CUCKOO SEARCH ALGORITHM (CASE STUDY OF MINIMIZING)
 *
 * 1 - RANDOMLY INITIALIZE A MATRIX --->
 *
 *POPULATION SIZE ---
 *NO OF ITERATION
 *Pa = 0.25
 *
 * ITERATIONS LOOP{
 *                  PHASE 1
 *
 * 2 - CHOOSE THE BEST
 *
 * MAKE A LOOP FOR ALL THE NESTS
 *
 * 3 - CARRY OUT THE FOLLOWING OPERATIONS FOR THE NESTS
 *
 *      GENERATE
 *       U = randn * sigmau
 *       v = randn
 *       Xnew = randn * 0.01*s*(X(t) - Best)
 *       where s = u/|v|^1/β
 *
 *       check if Xnew is within bounds, if it is not, replace it with ub or lb depending on deviaton
 * 4 - PERFORM THE GREEDY SELECTION
 *       if f(Xnew) < f(X){replacement occurs}
 *       else{no replacement occurs}
 *
 * 5 - Then display the updated Nest Solutions after the First Phase
 *
 *
 *          PHASE 2
 *    LOOP THROUGH ALL THE NESTS
 *
 *    1 - r (random number btw 0 & 1) ---
 *              r -  It is a matrix, generated for every Nest, depending on Number of Factors
 *
 *          if {r < Pa}
     *          X is selected and modified { checking for x1 and x2 individual}
     *          Xnew = X + rand*(Xd1 - Xd2)
     *       else{Nothing is Done}
     *
 *    2 - Carry out greedy selection
 *
 *    GO TO THE NEXT ITERATION --->2
 *
 *
 *    Output the details of the Cuckoo only for Iteration 1
 *    After Every Iteration, Cout the BEST
 *   }
 * */

const int population {5};
const int maxiter {20};
const double Pa {0.25};

//Defining Function Prototypes

vector<vector<double>>randomInit (vector<vector<double>>);
double functions(vector<double>);
void iterations_min(vector<vector<double>>);
void iterations_max(vector<vector<double>>);
vector<vector<double>>phase1(vector<vector<double>>);
vector<vector<double>>phase2(vector<vector<double>>);
void printVec(vector<vector<double>>);


int main()
{
    vector<vector<double>> factor_bounds {{-5,5},{-5,5}};
    vector<vector<double>> randInit = randomInit(factor_bounds);
    printVec(randInit);
     return 0;
}

vector<vector<double>> randomInit(vector<vector<double>> bounds) {
    //recall that the formular is x = L + r(U - L)
    //the vector bounds has the lower and upper boundaries...
    vector<vector<double>> init;
    vector<double>randoms{}; //Vector to store the random numbers for the first wolf

    double x;
    srand((unsigned)time(NULL)); //Seeding the random numbers so that they vary always

    for (int k{ 0 }; k < population; k++) {
        vector<double>init_part{};
        double random{};

        for (int i{ 0 }; i < bounds.size(); i++) {
            //generate the random numbers here
            random = (float)rand() / RAND_MAX; //generating random numbers between 0 and 1

            //Obtaining the random numbers only for the first wolf.
            if(k==0){
              randoms.push_back(random);
            }

            x = bounds.at(i).at(0) + random * (bounds.at(i).at(1) - bounds.at(i).at(0));
            init_part.push_back(x);
        }
        double fofx = functions(init_part);
        init_part.push_back(fofx); // you add the objective functions result of all the parameters into the vector
        init.push_back(init_part); // you push back one wolf into the initial eqn
    }
    cout<<"The random numbers generated for the first wolf is: \n";

    //displaying the random numbers for the first wolf.
    for(auto c:randoms){
        cout<<c<<"   ";
    }
    cout<<endl;

    return init;
}

double functions(vector<double> facs){
    double val;

    //Defining the factors
    double x1 = facs.at(0);
    double x2 = facs.at(1);
    //Enter the Objective Function Here
    val = pow(x1,2) - x1*x2 + pow(x2,2) + 2*x1 + 4*x2 + 3; //Example from the video

    return val;
}

void printVec(vector<vector<double>> wolf) {
    int wolfSize = wolf.at(0).size();
    for (int i{ 0 }; i < wolf.size(); i++) {
        for (int j{ 0 }; j < wolfSize; j++) {
            if (j == (wolfSize - 1)) {
                cout << "  ";
            }
            cout << setw(9) << left << wolf.at(i).at(j) << " ";
        }
        cout << endl;
    }
}
