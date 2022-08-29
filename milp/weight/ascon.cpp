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
#include <fstream>
#include <string>
string itos(int i) {stringstream s; s << i; return s.str(); }
#define N 320
GRBEnv env = GRBEnv();

// Constraints with weight 4
void sbox_w4(GRBModel& model, vector<GRBVar>& x, vector<GRBVar>& y, GRBVar z) {
    int i, j;
    for (i = 0; i < 37; i++) {
        GRBLinExpr tmp = 0;
        int sum = 0 ;
        for(j = 0; j<5; j++){
            tmp += W4[i][j]*x[j];
            if(W4[i][j] == -1){
                sum+=1;
            }
        }
        for(j = 5; j<10; j++){
            tmp += W4[i][j]*y[j-5];
            if(W4[i][j] == -1){
                sum+=1;
            }
        }
        model.addConstr(tmp - 10*z >= 1 - sum - 10 );
    }
}

// Constraints with weight 3
void sbox_w3(GRBModel& model, vector<GRBVar>& x, vector<GRBVar>& y, GRBVar z) {
    int i, j;
    for (i = 0; i < 50; i++) {
        GRBLinExpr tmp = 0;
        int sum = 0 ;
        for(j = 0; j<5; j++){
            tmp += W3[i][j]*x[j];
            if(W3[i][j] == -1){
                sum+=1;
            }
        }
        for(j = 5; j<10; j++){
            tmp += W3[i][j]*y[j-5];
            if(W3[i][j] == -1){
                sum+=1;
            }
        }
        model.addConstr(tmp - 10*z >= 1 - sum - 10 );
    }
}

// Constraints with weight 2
void sbox_w2(GRBModel& model, vector<GRBVar>& x, vector<GRBVar>& y, GRBVar z) {
    int i, j;
    for (i = 0; i < 22; i++) {
        GRBLinExpr tmp = 0;
        int sum = 0 ;
        for(j = 0; j<5; j++){
            tmp += W2[i][j]*x[j];
            if(W2[i][j] == -1){
                sum+=1;
            }
        }
        for(j = 5; j<10; j++){
            tmp += W2[i][j]*y[j-5];
            if(W2[i][j] == -1){
                sum+=1;
            }
        }
        model.addConstr(tmp - 10*z >= 1 - sum - 10 );
    }
}

int ascon_model(int threadNumber, int rounds, int *ID);
int ascon_model(int threadNumber, int rounds, int *ID) {
    int i, j, r ;
    try {
	    env.set(GRB_IntParam_LogToConsole, 1);
		env.set(GRB_IntParam_Threads, threadNumber);
		GRBModel model = GRBModel(env);

		// Variables 
		vector<vector<GRBVar>> X(rounds, vector<GRBVar>(320));
		vector<vector<GRBVar>> Y(rounds, vector<GRBVar>(320));
		vector<vector<GRBVar>> U(rounds-1, vector<GRBVar>(320));
		vector<vector<GRBVar>> D(rounds, vector<GRBVar>(64)) ;
		vector<vector<GRBVar>> w4(rounds, vector<GRBVar>(64)) ;
		vector<vector<GRBVar>> w3(rounds, vector<GRBVar>(64)) ;
		vector<vector<GRBVar>> w2(rounds, vector<GRBVar>(64)) ;
		vector<GRBVar> sx(5);
		vector<GRBVar> sy(5);

		for(r = 0 ; r<rounds; r++){
			for(i = 0; i<320; i++){
				X[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "X_"+itos(r)+"_"+itos(i));
				Y[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "Y_"+itos(r)+"_"+itos(i));
			}
			for(i = 0 ; i<64; i++){
				D[r][i]  = model.addVar(0, 1, 0, GRB_BINARY, "D_"+itos(r)+"_"+itos(i));
				w4[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "W4_"+itos(r)+"_"+itos(i));
				w3[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "W3_"+itos(r)+"_"+itos(i));
				w2[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "W2_"+itos(r)+"_"+itos(i));
			}
		}
		for(r = 0 ; r<rounds-1; r++){
			for(i = 0; i<320; i++){
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

				sx[0] = X[r][i] ;
				sx[1] = X[r][64 + i] ;
				sx[2] = X[r][128 + i] ;
				sx[3] = X[r][192 + i] ;
				sx[4] = X[r][256 + i] ;

				sy[0] = Y[r][i] ;
				sy[1] = Y[r][64 + i] ;
				sy[2] = Y[r][128 + i] ;
				sy[3] = Y[r][192 + i] ;
				sy[4] = Y[r][256 + i] ;

				model.addConstr(w4[r][i] + w3[r][i] + w2[r][i] == D[r][i]);
				sbox_w4(model, sx, sy, w4[r][i]);
				sbox_w3(model, sx, sy, w3[r][i]);
                sbox_w2(model, sx, sy, w2[r][i]);
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

        // Last round without linear layer
		for(i = 0 ; i<64; i++){
		    model.addConstr(X[r][i] + X[r][64 + i] + X[r][128 + i] + X[r][192 + i] + X[r][256 + i] >= D[r][i]);
		    model.addConstr(X[r][i] + X[r][64 + i] + X[r][128 + i] + X[r][192 + i] + X[r][256 + i] <= 5*D[r][i]);
		    model.addConstr(Y[r][i] + Y[r][64 + i] + Y[r][128 + i] + Y[r][192 + i] + Y[r][256 + i] >= D[r][i]);
		    model.addConstr(Y[r][i] + Y[r][64 + i] + Y[r][128 + i] + Y[r][192 + i] + Y[r][256 + i] <= 5*D[r][i]);

		    sx[0] = X[r][i] ;
		    sx[1] = X[r][64 + i] ;
		    sx[2] = X[r][128 + i] ;
		    sx[3] = X[r][192 + i] ;
		    sx[4] = X[r][256 + i] ;

		    sy[0] = Y[r][i] ;
		    sy[1] = Y[r][64 + i] ;
		    sy[2] = Y[r][128 + i] ;
		    sy[3] = Y[r][192 + i] ;
		    sy[4] = Y[r][256 + i] ;

		    model.addConstr(w4[r][i] + w3[r][i] + w2[r][i] == D[r][i]);
		    sbox_w4(model, sx, sy, w4[r][i]);
		    sbox_w3(model, sx, sy, w3[r][i]);
		    sbox_w2(model, sx, sy, w2[r][i]);
		}

		vector<GRBLinExpr> KS(rounds);
		GRBLinExpr KSS;
		int AS[rounds];
		for(j = 0; j<rounds; j++){
		    for(i = 0 ; i<64; i++){
		        KS[j] += D[j][i];
		    }
		}

		// Number of active sboxes constraint in each round
		for(i = 0 ; i<rounds; i++){
		    model.addConstr(KS[i] == ID[i]);
		}

		// Sum of active sboxes
		for(i = 0 ; i<rounds; i++){
		    KSS += KS[i];
		}

		vector<GRBLinExpr> WS(rounds);
		GRBLinExpr WSS;
		int W[rounds];

		// Weight in each round
        for(i = 0 ; i<rounds; i++){
            for(j = 0 ; j<64; j++){
                WS[i] += 4*w4[i][j] + 3*w3[i][j] + 2*w2[i][j] ;
            }
        }
        // Total weight
        for(i = 0 ; i<rounds; i++){
            WSS += WS[i];
        }

        model.setObjective(WSS, GRB_MINIMIZE);
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

		for(i = 0 ; i<rounds; i++){
		    W[i] = round(WS[i].getValue());
		}
		cout << "Weights: ";
		for(i = 0 ; i<rounds; i++){
		    cout << W[i] << " " ;
		}
		cout << endl ;

		int obj = round(model.getObjective().getValue());
		cout << "Total weight: " << obj ;
		cout << endl ;

		if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			return -1;
		}
		else if ((model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
		    return 1 ;
		}
		else {
			//cout << model.get(GRB_IntAttr_Status) << endl;
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

int minimum_weight(int rounds, int threadNumber) {
    int ID[3] = {1, 3, 11} ;
    ascon_model(threadNumber, rounds, ID);
    return 0;
 }
