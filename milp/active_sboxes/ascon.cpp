/*
Copyright 2022 Rusydi H. Makarim and Raghvendra Rohit

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include "main.h"
#include <sstream>
string itos(int i) {stringstream s; s << i; return s.str(); }
#define N 320
GRBEnv env = GRBEnv();

// Inequalities for Sbox
const int INEQ[50][11] = {{4, 9, 3, 9, 5, -2, 2, -1, 0, -6, 0},
{-2, 3, 0, 4, -2, 1, 1, 2, 0, 3, 0}, 
{13, 7, -5, -5, 8, -1, -2, 4, 4, 13, 0}, 
{-1, 3, 1, -3, -2, -1, 0, 0, -1, -2, 7}, 
{3, 5, 5, -2, 1, 3, 0, -2, -2, 3, 0}, 
{1, -1, 3, 5, 4, -2, 2, 6, -1, -2, 0}, 
{2, -2, -3, 2, 4, 0, -1, -4, 1, -1, 7}, 
{-1, -1, -1, -3, -2, 0, 2, 3, 3, 0, 5}, 
{-3, -2, -3, -4, -6, -1, -3, -3, -3, -1, 23}, 
{-3, 3, 2, 4, -2, 1, -2, -1, 0, 1, 4}, 
{-5, -4, -2, -3, -2, 1, 5, -7, -7, -3, 26}, 
{2, -1, -1, -2, -2, 0, 2, 0, 1, 0, 4}, 
{2, -2, -4, 3, 1, 4, -4, 1, -1, -4, 10}, 
{-3, -1, 2, -2, -3, -1, -2, 5, 5, -2, 9}, 
{-2, -1, -5, 5, -1, -4, -4, -3, 2, 6, 14}, 
{-2, 2, -2, -1, 1, 2, 0, 1, 1, 0, 3}, 
{4, 4, -3, 2, -2, 5, 1, 0, 3, 1, 0}, 
{-1, -2, -9, 7, 3, -7, 6, -4, 1, -5, 19},
{1, -2, 2, 1, 2, 0, 0, 1, 1, 1, 0}, 
{-3, -3, -4, -1, 2, 1, 4, -2, -2, 1, 11}, 
{-2, 9, 10, 8, -4, -3, 10, -1, -2, 12, 0}, 
{4, -3, -4, -1, 2, -1, -4, 0, -2, 1, 11}, 
{2, -3, 2, -1, 1, -1, -1, -3, 3, -3, 9}, 
{0, 2, -1, 2, 1, -1, -1, 1, 0, -1, 2}, 
{1, -1, -2, 3, -1, 4, 4, 2, 0, 4, 0}, 
{2, 2, -2, -3, -4, -4, 1, -1, -1, -1, 12}, 
{0, -1, 0, -1, 1, 0, 1, 1, -1, -1, 3}, 
{-1, -3, 2, 2, 1, -1, -1, 3, -1, 2, 4}, 
{-1, -3, -1, 2, -1, 3, -3, 1, -1, -2, 9}, 
{2, 1, 1, -3, 1, 2, 0, 2, 2, 2, 0}, 
{-2, 2, 3, 1, -1, -1, -3, -1, 0, -1, 6}, 
{-1, -2, 2, -2, 1, 0, -2, -1, -2, 2, 8}, 
{-2, -2, 1, 1, -1, 1, -1, 3, 2, -3, 6}, 
{-1, -1, 0, -1, 1, 0, 1, 0, 1, 1, 2}, 
{-1, 1, -1, -2, -2, 0, 0, 1, 1, -1, 5}, 
{1, -2, 2, -1, -1, -1, -1, 1, 2, 1, 4}, 
{-1, 1, -1, -1, 0, 1, 0, -1, -1, 0, 4}, 
{0, 0, 1, 1, -1, 1, 1, -1, 0, -1, 2}, 
{1, 0, 2, 1, -1, 2, -2, -2, 1, 2, 3},
{2, -4, -4, 2, 3, -1, -1, -1, -2, 1, 9}, 
{-1, -2, 0, 2, -1, -2, 2, 1, -1, -1, 6}, 
{2, -2, -1, -1, 1, 0, -1, 1, -1, -1, 5}, 
{-1, 2, 1, -2, -1, -1, 0, 0, 0, -1, 4}, 
{-1, 1, 2, 2, -2, -1, 2, -2, 1, 3, 3}, 
{2, 5, 1, -3, -5, -4, 1, -1, -4, -1, 13}, 
{2, -2, -2, 2, 3, 1, -1, -3, 2, -1, 5}, 
{0, 0, 1, 1, -1, -1, -1, -1, 0, -1, 4}, 
{0, -1, 1, -1, 1, 0, -1, -1, 1, -1, 4}, 
{-1, 1, -2, 1, 2, 1, -1, 2, 0, -1, 3}, 
{2, -1, -1, -1, 1, 0, 0, -2, -2, 2, 5}} ; 

int ascon_model(int threadNumber, int rounds);
int ascon_model(int threadNumber, int rounds) {
    int i, j, r ;
    try {
	    env.set(GRB_IntParam_LogToConsole, 1);
		env.set(GRB_IntParam_Threads, threadNumber);
		GRBModel model = GRBModel(env);

		// Variables 
		vector<vector<GRBVar>> X(rounds, vector<GRBVar>(320));
		vector<vector<GRBVar>> Y(rounds-1, vector<GRBVar>(320));
		vector<vector<GRBVar>> D(rounds, vector<GRBVar>(64)) ;
		vector<vector<GRBVar>> U(rounds-1, vector<GRBVar>(320));

		for(r = 0 ; r<rounds; r++){
			for(i = 0; i<320; i++){
				X[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "X_"+itos(r)+"_"+itos(i)); 
			}
			for(i = 0 ; i<64; i++){
				D[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "D_"+itos(r)+"_"+itos(i));
			}
		}

		for(r = 0 ; r<rounds-1; r++){
			for(i = 0; i<320; i++){
				Y[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "Y_"+itos(r)+"_"+itos(i));
				U[r][i] = model.addVar(0, 2, 0, GRB_INTEGER, "U_"+itos(r)+"_"+itos(i));
			}
		}
		
		for(r = 0; r<rounds-1; r++){
		    // Sbox constraints
		    for(i = 0 ; i<64; i++){
		        model.addConstr(X[r][i] + X[r][64 + i] + X[r][128 + i] + X[r][192 + i] + X[r][256 + i] >= D[r][i]);
				model.addConstr(X[r][i] + X[r][64 + i] + X[r][128 + i] + X[r][192 + i] + X[r][256 + i] <= 5*D[r][i]);
				model.addConstr(Y[r][i] + Y[r][64 + i] + Y[r][128 + i] + Y[r][192 + i] + Y[r][256 + i] >= D[r][i]);
				model.addConstr(Y[r][i] + Y[r][64 + i] + Y[r][128 + i] + Y[r][192 + i] + Y[r][256 + i] <= 5*D[r][i]);

				for(j = 0 ; j < 50; j++){
					GRBLinExpr tt = 0 ; 
					for(int k = 0 ; k<5; k++) tt += INEQ[j][k] * X[r][64*k+i];
					for(int k = 0 ; k<5; k++) tt += INEQ[j][5 + k] * Y[r][64*k+i];
					tt += INEQ[j][10] ; 
					model.addConstr(tt >= 0) ; 
				}
			}
		    // Linear layer constraints
			for(i = 0; i<64; i++){
			    model.addConstr(Y[r][0   + i] + Y[r][((64+i-19) % 64) + 0]   + Y[r][((64+i-28) % 64) + 0]   + X[r+1][i]      == 2*U[r][i]) ;
				model.addConstr(Y[r][64  + i] + Y[r][((64+i-61) % 64) + 64]  + Y[r][((64+i-39) % 64) + 64]  + X[r+1][64+i]   == 2*U[r][64+i]) ;
				model.addConstr(Y[r][128 + i] + Y[r][((64+i-1)  % 64) + 128] + Y[r][((64+i-6)  % 64) + 128] + X[r+1][128+i]  == 2*U[r][128+i]) ;
				model.addConstr(Y[r][192 + i] + Y[r][((64+i-10) % 64) + 192] + Y[r][((64+i-17) % 64) + 192] + X[r+1][192 +i] == 2*U[r][192+i]) ;
				model.addConstr(Y[r][256 + i] + Y[r][((64+i-7)  % 64) + 256] + Y[r][((64+i-41) % 64) + 256] + X[r+1][256 +i] == 2*U[r][256+i]) ;
			}
		}

		// Last round active Sboxes constraint
		for(i = 0 ; i<64; i++){
		    model.addConstr(X[r][i] + X[r][64 + i] + X[r][128 + i] + X[r][192 + i] + X[r][256 + i] >= D[r][i]);
		    model.addConstr(X[r][i] + X[r][64 + i] + X[r][128 + i] + X[r][192 + i] + X[r][256 + i] <= 5*D[r][i]);
		}

		vector<GRBLinExpr> KS(rounds);
		GRBLinExpr KSS;
		int AS[rounds];

		for(j = 0; j<rounds; j++){
		    for(i = 0 ; i<64; i++){
		        KS[j] += D[j][i];
		    }
		}
		//model.addConstr(D[0][0] == 1) ;         // Zero-th Sbox is fixed to be active to avoid rotated solutions.
		model.addConstr(KS[0] >= 1);
		for(i = 0 ; i<rounds; i++){
		    KSS += KS[i];
		}
		model.setObjective(KSS, GRB_MINIMIZE);

		// 5 round trail with 72 active Sboxes
		/*
         * Input State X[0], X[1] and X[2]
         * Extend from round 0 to round 4
         */

		char state72_0[321] = "1......................................1.1...................1..1......1...............................1.1...................1..1......1...............................1.1...................1.........1.................................1.............................1........................................................" ;
		char state72_1[321] = "....................................................................1..1........1.....................1..1....1.................................................................................................................................................1...1..1........1.....................11.1....1..............1..";
		char state72_2[321] = "................................................................1.............1...........................................1.....................................................................................................................................1..........1..11..1....................1.............1...1...1..";

		for(i = 0 ; i<320; i++){
            if(state72_0[i] == '1'){
                model.addConstr(X[0][i] == 1) ;
            }
            else{
                model.addConstr(X[0][i] == 0) ;
            }
		}
		for(i = 0 ; i<320; i++){
		    if(state72_1[i] == '1'){
		        model.addConstr(X[1][i] == 1) ;
		    }
		    else{
		        model.addConstr(X[1][i] == 0) ;
		    }
		}
		for(i = 0 ; i<320; i++){
		    if(state72_2[i] == '1'){
		        model.addConstr(X[2][i] == 1) ;
		    }
		    else{
		        model.addConstr(X[2][i] == 0) ;
		    }
		}



		// 5 round trail with 75 active Sboxes
		/*
         * Input State X[0], X[1] and X[2]
         * Extend from round 0 to round 4
         */
		/*
		char state75_0[321] = ".......1......1.1..1..1......1..1...1..................................1......1....1..1......1..1......1.....................1........................1.........1......1........................................1..1................1..................................1......1.1............1......1........................1..";
		char state75_1[321] = "................................................................................1.....1...1.....11..1..1......1......1..........................................................................................1.....1...1.....11..1..1......1......1.........................................................................." ;
		char state75_2[321] = "......................................................................................................................................................................................................1.........1.....1...................1.......1.............................................................................";

		for(i = 0 ; i<320; i++){
            if(state75_0[i] == '1'){
                model.addConstr(X[0][i] == 1) ;
            }
            else{
                model.addConstr(X[0][i] == 0) ;
            }
        }

		for(i = 0 ; i<320; i++){
		    if(state75_1[i] == '1'){
		        model.addConstr(X[1][i] == 1) ;
		    }
		    else{
		        model.addConstr(X[1][i] == 0) ;
		    }
		}

		for(i = 0 ; i<320; i++){
		    if(state75_2[i] == '1'){
		        model.addConstr(X[2][i] == 1) ;
		    }
		    else{
		        model.addConstr(X[2][i] == 0) ;
		    }
		}
		*/

		model.optimize();
		// Print the differential trail
		for(r= 0; r<rounds-1; r++){
		    for (i = 0; i < 5; i++) {
				for(j = 0; j<64; j++){
					if (round(X[r][64*i+j].get(GRB_DoubleAttr_Xn)) == 1) cout<<"1";
					else cout<<".";
					}
				cout << endl ;
			}
			cout << "----------------------------------------------------------------"<<"  p_S"<<endl ;
			for (i = 0; i < 5; i++) {
				for(j = 0; j<64; j++){
					if (round(Y[r][64*i+j].get(GRB_DoubleAttr_Xn)) == 1) cout<<"1";
					else cout<<".";
					}
				cout << endl ;
			}
			cout << "----------------------------------------------------------------"<<"  p_L" <<endl ;
		}
		for(i = 0; i<5; i++){
     		for(j = 0 ; j<64; j++) {
     			if (round(X[rounds-1][64*i+j].get(GRB_DoubleAttr_Xn)) == 1) cout<<"1";
					else cout<<".";
     		}
     		cout << endl; 
     	}
		cout << "----------------------------------------------------------------" <<endl ;

        for(i = 0 ; i<rounds; i++){
            AS[i] = round(KS[i].getValue());
        }
        cout << "Configuration: ";
        for(i = 0 ; i<rounds; i++){
            cout << AS[i] << " " ;
        }
        cout << endl ;

		int obj = round(model.getObjective().getValue());
		cout << "Sum: " << obj ;
		cout << endl ;

		if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			return -1;
		}
		else if ((model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
		    return 1 ;
		}
		else {
			return -2;
		}
	}
	catch (GRBException e) {
		cerr << "Error code = " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	
	catch (...) {
		cerr << "Exception during optimization" << endl;
	}
	return -1;
}

int min_sbox(int rounds, int threadNumber) {
    ascon_model(threadNumber, rounds);
    return 0;
 }
