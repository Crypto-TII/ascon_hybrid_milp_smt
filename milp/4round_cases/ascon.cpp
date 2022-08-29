#include "main.h"
#include <sstream>
#include <fstream>
#include <string>
string itos(int i) {stringstream s; s << i; return s.str(); }
GRBEnv env = GRBEnv();

int ascon_model(int threadNumber);
int ascon_model(int threadNumber) {
    int i, j, r ;
	
	try {
	    env.set(GRB_IntParam_LogToConsole, 1);
		env.set(GRB_IntParam_Threads, threadNumber);
		env.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_BESTBOUND);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 100000000);
		env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);

		GRBModel model = GRBModel(env);

		// Variables 
		GRBVar a, b, c, d;
		a = model.addVar(2, 23, 0, GRB_INTEGER, "A");
		b = model.addVar(2, 43, 0, GRB_INTEGER, "B");
		c = model.addVar(2, 43, 0, GRB_INTEGER, "C");
		d = model.addVar(2, 24, 0, GRB_INTEGER, "D");

		model.addConstr(a + b + c >= 16);
		model.addConstr(b + c + d >= 16);
		model.addConstr(a + b + c + d >= 24);
		model.addConstr(a + b + c + d <= 44);
		model.addConstr(a + b >= 6);
		model.addConstr(b + c >= 6);
		model.addConstr(c + d >= 6);
		model.addConstr(a + c >= 6);
		model.addConstr(b + d >= 6);
		model.addConstr(b <= 34);
		model.optimize();

		int solCount = model.get(GRB_IntAttr_SolCount);
		cout << "Number of solutions:" << solCount<<endl;
		if (solCount >= 2000000000) {
		    cerr << "Number of solutions is too large" << endl;
		    exit(0);
		}

		int aa, bb, cc, dd ;
		for(i = 0 ; i<solCount; i++){
		    model.set(GRB_IntParam_SolutionNumber, i);
		    aa = round(a.get(GRB_DoubleAttr_Xn));
		    bb = round(b.get(GRB_DoubleAttr_Xn));
		    cc = round(c.get(GRB_DoubleAttr_Xn));
		    dd = round(d.get(GRB_DoubleAttr_Xn));
		    cout << aa << " " << bb << " " << cc << " " << dd << endl ;

		}
        cout << solCount << endl ;
		if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			return -1;
		}
		else if ((model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
		    return 1 ;

		}
		else {
			cout << model.get(GRB_IntAttr_Status) << endl;
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

int possible_cases(int threadNumber){
    ascon_model(threadNumber);
    return 0;
 }


