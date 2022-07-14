#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <time.h>
#include <cstdlib>
#include <iomanip>

using namespace std;

const int population {10};
const int maxiter {400};
const double Pa {0.25};

//Putting down the bounds for LD, W and L
vector<vector<double>> factor_bounds {{5,20},{4.5,6.0},{60,90},{3,6},{15,30}}; //Defining the boundaries as a universal variable

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

    vector<vector<double>> randInit = randomInit(factor_bounds);
    cout<<"\tRANDOMLY INITIALIZED MATRIX\n";
    printVec(randInit);

    int choice{};
    cout<<"1 - MAXIMIZE\n2 - MINIMIZE\nCHOICE: ";
    cin>>choice;

    if(choice == 1){
        cout<<"Maximizing...\n";
        iterations_max(randInit);
    }else if(choice == 2){
        cout<<"Minimizing...\n";
        iterations_min(randInit);
    }else{
        cout<<"Kindly enter a valid input...";
    }

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
    double F = facs.at(0);
    double DC = facs.at(1);
    double PO = facs.at(2);
    double N = facs.at(3);
    double GS = facs.at(4);

    //Objective function for the Surface Roughness
    val =-1.53*F + 2.76*DC + 0.169*PO - 1.92*N - 0.565*GS + 0.00363*F*F + 0.222*DC*DC
+ 0.001081*PO*PO + 0.179*N*N + 0.01285*GS*GS + 0.086*F*DC + 0.00592*F*PO + 0.0630*F*N
+ 0.0080*F*GS - 0.0683*DC*PO;

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
            if(i == 1){
                cout<<endl<<"ITERATIONS 2 - "<<maxiter<<endl;
            }

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
        if(i == 1){
            cout<<endl<<"ITERATIONS 2 - "<<maxiter<<endl;
        }

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
    srand(time(nullptr));
    std::default_random_engine generator(rand());
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
    srand(time(nullptr));
    std::default_random_engine generator(rand());
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

                double random = distribution1(generator);
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

    srand(time(nullptr));
    std::default_random_engine generator(rand());
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
    srand(time(nullptr));
    std::default_random_engine generator(rand());
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
                int d1 = rand() % population ;
                int d2 = rand() % population ;

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
