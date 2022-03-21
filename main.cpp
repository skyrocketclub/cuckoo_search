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

const int population {100};
const int maxiter {200};
const double Pa {0.25};
vector<vector<double>> factor_bounds {{-5,5},{-5,5}}; //Defining the boundaries as a universal variable

//Defining Function Prototypes

vector<vector<double>>randomInit (vector<vector<double>>);
double functions(vector<double>);
void iterations_min(vector<vector<double>>);
void iterations_max(vector<vector<double>>);
vector<vector<double>>phase1_max(vector<vector<double>>&,int);
vector<vector<double>>phase2_max(vector<vector<double>>&,int);
vector<vector<double>>phase1_min(vector<vector<double>>&,int);
vector<vector<double>>phase2_min(vector<vector<double>>&,int);
void printVec(vector<vector<double>>);
void printVec(vector<double>);
vector<double> find_best_max(vector<vector<double>>);
vector<double> find_best_min(vector<vector<double>>);


int main()
{
//    vector<vector<double>> factor_bounds {{-5,5},{-5,5}};
    vector<vector<double>> randInit = randomInit(factor_bounds);
    cout<<"\tRANDOMLY INITIALIZED MATRIX\n";
    printVec(randInit);
    iterations_min(randInit);
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

void printVec(vector<double>vec){
        for(auto c:vec){
            cout<<c<<" ";
        }
        cout<<endl;
}

vector<double> find_best_max(vector<vector<double>> facs){
    vector<double>best;
    best = facs.at(0);
    /*
     * 1 - Loop through the facs
     * 2 - Assign the first vec as the best vec
     * 3 - Compare the rest of the vec to the best vec and reassign if any of them beat the best vec
     * */
    size_t funcpos = facs.at(0).size() - 1;
    for(size_t i{0}; i<facs.size(); i++){
        if(facs.at(i).at(funcpos)>best.at(funcpos)){
            best = facs.at(i);
        }
    }
    return best;
}

vector<double> find_best_min(vector<vector<double>> facs){
    vector<double>best;
    best = facs.at(0);
    size_t funcpos = facs.at(0).size() - 1;
    for(size_t i{0}; i<facs.size(); i++){
        if(facs.at(i).at(funcpos)<best.at(funcpos)){
            best = facs.at(i);
        }
    }
    return best;
}

void iterations_max(vector<vector<double>>facs){

        for(int i{0}; i<maxiter; i++){
            facs = phase1_max(facs,i);
            if(i == 0){
                cout<<endl;
                printVec(facs);
            }
            facs = phase2_max(facs,i);

            if(i == 0){ cout<<endl; printVec(facs);}

            if(i>0){
                vector<double> best = find_best_max(facs);
                size_t bs = best.size();
                for(size_t j{0}; j<bs; j++){

                    cout<<best.at(j)<<",";
            }
               cout<<endl;
        }

    }
}

void iterations_min(vector<vector<double>>facs){
    for(int i{0}; i<maxiter; i++){
        facs = phase1_min(facs,i);
        if(i == 0){
            cout<<endl;
            printVec(facs);
        }
        facs = phase2_min(facs,i);

        if(i == 0){ cout<<endl; printVec(facs);}

        if(i>0){
            vector<double> best = find_best_min(facs);
            size_t bs = best.size();
            for(size_t j{0}; j<bs; j++){

                cout<<best.at(j)<<",";
        }
           cout<<endl;
    }

}
}

vector<vector<double>>phase1_max(vector<vector<double>>&facs, int num){
/*
 *  2 - CHOOSE THE BEST
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
 * */
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);
    std::uniform_real_distribution<double> distribution1(0,1);
    double sigmau = 0.6966;

    vector<double>best = find_best_max(facs);
    if(num == 0){
        cout<<endl;
        cout<<"Best: ";
        printVec(best);
    }

    //Looping through all the Nests...

    for(size_t i{0}; i<facs.size(); i++){

        vector<double>Xnew;

        //Looping through Nest one by one to find the Xnew Values
        for(size_t j{0}; j<facs.at(0).size()-1; j++){ //The Function is not included

            double U = distribution(generator)*sigmau;
            double v = distribution(generator);
            double s = U/pow(abs(v),1/1.5);

            if(num==0){
                cout<<endl;
                cout<<"For Nest "<<i+1<<" Factor "<<j<<endl;
                cout<<"U: "<<U<<" ; v: "<<v<<" ;s: "<<s<<endl;
            }

            double new_fac{};
            double randval = distribution(generator);

            if(num==0){
                cout<<endl;
                cout<<"For Xnew, randn: "<<randval<<endl;
            }

            new_fac = facs.at(i).at(j) + randval * 0.01 *s*(facs.at(i).at(j) - best.at(j));
            if(num == 0){
                cout<<"NewFac: "<<new_fac<<endl;
            }

           //Here you have to confirm that the newfac generated is within the bounds
            if(new_fac>=factor_bounds.at(j).at(0)&& new_fac<= factor_bounds.at(j).at(1)){

                if(num == 0){
                    cout<<"NewFac: "<<std::setprecision(7)<<std::showpoint<<new_fac<<" is within the bounds"<<endl;
                }
            }else if(new_fac<factor_bounds.at(j).at(0)){

                if(num == 0){
                    cout<<"NewFac: "<<std::setprecision(7)<<std::showpoint<<new_fac<<" is lower than the lower bounds"<<endl;
                }
                 new_fac = factor_bounds.at(j).at(0);
            }else{
                if(num == 0){
                    cout<<"NewFac: "<<std::setprecision(7)<<std::showpoint<<new_fac<<" is higher than the higher bounds"<<endl;
                }
                  new_fac = factor_bounds.at(j).at(1);
            }

            Xnew.push_back(new_fac); //All the factors will have their new value here
        }
           double func = functions(Xnew);
           Xnew.push_back(func);

           size_t pos_func = facs.at(i).size() -1;

           if(func>facs.at(i).at(pos_func)){
               facs.at(i)=Xnew;
               if(num == 0){
                    cout<<endl;
                   cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                   cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                   cout<<"Since F(Xnew) > F(X), Replacement will take place\n";
               }
           }
            else if (func<facs.at(i).at(pos_func)){

                   if(num == 0){
                       cout<<endl;
                       cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                       cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                       cout<<"Since F(Xnew) < F(X), Replacement will not take place\n";
                   }
               }
           else if(func==facs.at(i).at(facs.at(i).size()-1)) {

              if(num == 0){
                    cout<<endl;
                  cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                  cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                  cout<<"Since F(Xnew) = F(X), Replacement will not take place\n";
              }
           }
    }

    return facs;
}

vector<vector<double>>phase2_max(vector<vector<double>>&facs, int num){
/*
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
 * */

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);
    std::uniform_real_distribution<double> distribution1(0,1);


    srand(time(nullptr));
    if(num==0){cout<<"PHASE 2\n";}

    //Looping through all the Nests...
    for(size_t i{0}; i<facs.size(); i++){
     vector<double> Xnew;
     if(num==0){
     cout<<"NEST "<<i<<endl;
     }
        for(size_t j{0}; j<facs.at(i).size()-1; j++){
           if(num==0){cout<<"x"<<j+1<<endl;}

            double r = distribution1(generator);
            if(num == 0){cout<<"r: "<<r<<endl;}

            if(r<Pa){
               //This means that this factor will be manipulated
                int d1 = rand()%3 +1;
                int d2 = rand()%3 + 1;

                double random = distribution(generator);
                double new_fac{0};

                if(num==0){cout<<"r<Pa; d1:"<<d1<<" d2:"<<d2<<endl;}
                double current = facs.at(i).at(j);
                new_fac = current + random*(facs.at(d1).at(j) - facs.at(d2).at(j));

                if(num==0){cout<<"X"<<j+1<<" changed from "<<current<<" to "<< new_fac<<std::endl;}

                //Check if the value is within bounds
                //Here you have to confirm that the newfac generated is within the bounds
                 if(new_fac>=factor_bounds.at(j).at(0)&& new_fac<= factor_bounds.at(j).at(1)){

                     if(num == 0){
                         cout<<"NewFac: "<<std::setprecision(9)<<std::showpoint<<new_fac<<" is within the bounds"<<endl;
                     }
                 }else if(new_fac<factor_bounds.at(j).at(0)){

                     if(num == 0){
                         cout<<"NewFac: "<<std::setprecision(9)<<std::showpoint<<new_fac<<" is lower than the lower bounds"<<endl;
                     }
                      new_fac = factor_bounds.at(j).at(0);
                 }else{
                     if(num == 0){
                         cout<<"NewFac: "<<std::setprecision(9)<<std::showpoint<<new_fac<<" is higher than the higher bounds"<<endl;
                     }
                       new_fac = factor_bounds.at(j).at(1);
                 }


                Xnew.push_back(new_fac);
            }
            else{
                Xnew.push_back(facs.at(i).at(j));
            }
        }
        /*
         * 1 - Carry out Greedy Selection here
         * */
        double func = functions(Xnew);
        Xnew.push_back(func);
        size_t pos_func = facs.at(i).size() -1;

        if(func>facs.at(i).at(pos_func)){
            facs.at(i)=Xnew;
            if(num == 0){
                 cout<<endl;
                cout<<"F(Xnew): "<<std::setprecision(9)<<std::showpoint<<func<<endl;
                cout<<"F(X): "<<std::setprecision(9)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                cout<<"Since F(Xnew) > F(X), Replacement will take place\n";
            }
        }
         else if (func<facs.at(i).at(pos_func)){

                if(num == 0){
                    cout<<endl;
                    cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                    cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                    cout<<"Since F(Xnew) < F(X), Replacement will not take place\n";
                }
            }
        else if(func==facs.at(i).at(facs.at(i).size()-1)) {

           if(num == 0){
                cout<<endl;
               cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
               cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
               cout<<"Since F(Xnew) = F(X), Replacement will not take place\n";
           }
        }
    }
    return facs;
}

vector<vector<double>>phase1_min(vector<vector<double>>&facs, int num){

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);
    std::uniform_real_distribution<double> distribution1(0,1);
    double sigmau = 0.6966;

    vector<double>best = find_best_min(facs);
    if(num == 0){
        cout<<endl;
        cout<<"Best: ";
        printVec(best);
    }

    //Looping through all the Nests...

    for(size_t i{0}; i<facs.size(); i++){

        vector<double>Xnew;

        //Looping through Nest one by one to find the Xnew Values
        for(size_t j{0}; j<facs.at(0).size()-1; j++){ //The Function is not included

            double U = distribution(generator)*sigmau;
            double v = distribution(generator);
            double s = U/pow(abs(v),1/1.5);

            if(num==0){
                cout<<endl;
                cout<<"For Nest "<<i+1<<" Factor "<<j<<endl;
                cout<<"U: "<<U<<" ; v: "<<v<<" ;s: "<<s<<endl;
            }

            double new_fac{};
            double randval = distribution(generator);

            if(num==0){
                cout<<endl;
                cout<<"For Xnew, randn: "<<randval<<endl;
            }

            new_fac = facs.at(i).at(j) + randval * 0.01 *s*(facs.at(i).at(j) - best.at(j));
            if(num == 0){
                cout<<"NewFac: "<<new_fac<<endl;
            }

           //Here you have to confirm that the newfac generated is within the bounds
            if(new_fac>=factor_bounds.at(j).at(0)&& new_fac<= factor_bounds.at(j).at(1)){

                if(num == 0){
                    cout<<"NewFac: "<<std::setprecision(7)<<std::showpoint<<new_fac<<" is within the bounds"<<endl;
                }
            }else if(new_fac<factor_bounds.at(j).at(0)){

                if(num == 0){
                    cout<<"NewFac: "<<std::setprecision(7)<<std::showpoint<<new_fac<<" is lower than the lower bounds"<<endl;
                }
                 new_fac = factor_bounds.at(j).at(0);
            }else{
                if(num == 0){
                    cout<<"NewFac: "<<std::setprecision(7)<<std::showpoint<<new_fac<<" is higher than the higher bounds"<<endl;
                }
                  new_fac = factor_bounds.at(j).at(1);
            }

            Xnew.push_back(new_fac); //All the factors will have their new value here
        }
           double func = functions(Xnew);
           Xnew.push_back(func);

           size_t pos_func = facs.at(i).size() -1;

           if(func<facs.at(i).at(pos_func)){
               facs.at(i)=Xnew;
               if(num == 0){
                    cout<<endl;
                   cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                   cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                   cout<<"Since F(Xnew) < F(X), Replacement will take place\n";
               }
           }
            else if (func>facs.at(i).at(pos_func)){

                   if(num == 0){
                       cout<<endl;
                       cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                       cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                       cout<<"Since F(Xnew) > F(X), Replacement will not take place\n";
                   }
               }
           else if(func==facs.at(i).at(facs.at(i).size()-1)) {

              if(num == 0){
                    cout<<endl;
                  cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                  cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                  cout<<"Since F(Xnew) = F(X), Replacement will not take place\n";
              }
           }
    }

    return facs;
}

vector<vector<double>>phase2_min(vector<vector<double>>&facs, int num){
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,1);
    std::uniform_real_distribution<double> distribution1(0,1);


    srand(time(nullptr));
    if(num==0){cout<<"PHASE 2\n";}

    //Looping through all the Nests...
    for(size_t i{0}; i<facs.size(); i++){
     vector<double> Xnew;
     if(num==0){
     cout<<"NEST "<<i<<endl;
     }
        for(size_t j{0}; j<facs.at(i).size()-1; j++){
           if(num==0){cout<<"x"<<j+1<<endl;}

            double r = distribution1(generator);
            if(num == 0){cout<<"r: "<<r<<endl;}

            if(r<Pa){
               //This means that this factor will be manipulated
                int d1 = rand()%3 +1;
                int d2 = rand()%3 + 1;

                double random = distribution(generator);
                double new_fac{0};

                if(num==0){cout<<"r<Pa; d1:"<<d1<<" d2:"<<d2<<endl;}
                double current = facs.at(i).at(j);
                new_fac = current + random*(facs.at(d1).at(j) - facs.at(d2).at(j));

                if(num==0){cout<<"X"<<j+1<<" changed from "<<current<<" to "<< new_fac<<std::endl;}

                //Check if the value is within bounds
                //Here you have to confirm that the newfac generated is within the bounds
                 if(new_fac>=factor_bounds.at(j).at(0)&& new_fac<= factor_bounds.at(j).at(1)){

                     if(num == 0){
                         cout<<"NewFac: "<<std::setprecision(9)<<std::showpoint<<new_fac<<" is within the bounds"<<endl;
                     }
                 }else if(new_fac<factor_bounds.at(j).at(0)){

                     if(num == 0){
                         cout<<"NewFac: "<<std::setprecision(9)<<std::showpoint<<new_fac<<" is lower than the lower bounds"<<endl;
                     }
                      new_fac = factor_bounds.at(j).at(0);
                 }else{
                     if(num == 0){
                         cout<<"NewFac: "<<std::setprecision(9)<<std::showpoint<<new_fac<<" is higher than the higher bounds"<<endl;
                     }
                       new_fac = factor_bounds.at(j).at(1);
                 }


                Xnew.push_back(new_fac);
            }
            else{
                Xnew.push_back(facs.at(i).at(j));
            }
        }
        /*
         * 1 - Carry out Greedy Selection here
         * */
        double func = functions(Xnew);
        Xnew.push_back(func);
        size_t pos_func = facs.at(i).size() -1;

        if(func<facs.at(i).at(pos_func)){
            facs.at(i)=Xnew;
            if(num == 0){
                 cout<<endl;
                cout<<"F(Xnew): "<<std::setprecision(9)<<std::showpoint<<func<<endl;
                cout<<"F(X): "<<std::setprecision(9)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                cout<<"Since F(Xnew) < F(X), Replacement will take place\n";
            }
        }
         else if (func>facs.at(i).at(pos_func)){

                if(num == 0){
                    cout<<endl;
                    cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
                    cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
                    cout<<"Since F(Xnew) > F(X), Replacement will not take place\n";
                }
            }
        else if(func==facs.at(i).at(facs.at(i).size()-1)) {

           if(num == 0){
                cout<<endl;
               cout<<"F(Xnew): "<<std::setprecision(7)<<std::showpoint<<func<<endl;
               cout<<"F(X): "<<std::setprecision(7)<<std::showpoint<< facs.at(i).at(pos_func)<<endl;
               cout<<"Since F(Xnew) = F(X), Replacement will not take place\n";
           }
        }
    }
    return facs;
}
