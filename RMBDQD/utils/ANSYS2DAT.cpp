#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<vector>
#include<sstream>

template<class T>
void insertionsort(std::vector<T>& x){

    const int n = x.size();

    if ( n == 1 ) {
        return;
    }

    for ( int i=1; i<n; i++ ) {

        T val = x[i];
        
        for ( int j=i-1; j>=0; j-- ) {
            
            if ( val < x[j] ) {
                x[j+1] = x[j];
                x[j] = val;
            }
            else {
                break;
            }

        }
    }

}

void vprint(std::vector<auto>& x){

    const int n = x.size();
    std::cout << "[";

    for (int i=0; i<n-1; i++){
        std::cout << x[i] << ", ";
    }

    std::cout << x[n-1] << "]" << std::endl;

}


using namespace std;
typedef vector< vector<double> > matd;
typedef vector< vector<int> > mati;
typedef vector<double> vecd;
typedef vector<int> veci;

int main(int argc, const char *argv[])
{

    if (argc < 3) {
        cout << "Precisa de dois argumentos: input e output\n";
        exit(1);
    }


    int itmp;

    ifstream inputfile(argv[1]);

    matd nodes;
    mati conec;

    int nNodes = 0, nEls = 0;

    while (true) {
        inputfile >> itmp;
        if (itmp == -1) break;

        vecd tmp;
        double x;
        nNodes++;

        for (int i=0; i<3; i++ ) {
            inputfile >> x;
            tmp.push_back(x);
        }

        nodes.push_back(tmp);
    }

    cout << "nNodes: " << nNodes << endl;

    while (true) {

        inputfile >> itmp;
        if (itmp == -1) break;

        int xx;
        int x1, x2, x3, x4;

        nEls++;

        for (int i=0; i<10; i++ ) {
            inputfile >> xx;
        }

        veci iitmp;

        inputfile >> x1 >> x2;
        for (int i=0; i<2; i++) inputfile >> x3;
        for (int i=0; i<4; i++) inputfile >> x4;


        iitmp.push_back(x1);
        iitmp.push_back(x2);
        iitmp.push_back(x3);
        iitmp.push_back(x4);
        
        conec.push_back(iitmp);

    }

    cout << "nEls: " << nEls << endl;

    inputfile.close();

    // =========================================================================

    int nBC = 0;
    veci BC;

    for (int i=0; i<nNodes; i++){
        if (abs(nodes[i][2]) < 1e-12) {
            nBC++;
            BC.push_back(i+1);
        }
    }

    cout << "nBC: " << nBC << endl;
    //insertionsort(BC);
    //vprint(BC);
    //
    // ================= write dat ============================
    
    ofstream of(argv[2]);

    of << nEls << "\t" << nNodes << "\t" << 6<< endl;
    of << 2.16e6 << "\t" << 0.34 << "\t" << 609<< endl << endl;

    for (int i=0; i<nNodes; i++) {
        of << scientific << setprecision(15) <<
            setw(22) << nodes[i][0] << "\t" <<
            setw(22) << nodes[i][1] << "\t" <<
            setw(22) << nodes[i][2] << endl;
    }

    of << endl;

    for (int i=0; i<nEls; i++) {
        of << 
            setw(12) << conec[i][0] << "\t" <<
            setw(12) << conec[i][1] << "\t" <<
            setw(12) << conec[i][2] << "\t" <<
            setw(12) << conec[i][3] << endl;
    }

    of << endl;

    of << nBC << endl;

    for (int i=0; i<nBC; i++){
        of <<
            BC[i] << "\t0\t0\t0" << endl;
    }

    of.close();
    
    return 0;
}
