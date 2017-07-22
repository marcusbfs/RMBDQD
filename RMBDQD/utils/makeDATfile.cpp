#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

int main(int argc, const char *argv[])
{
    if (argc < 2 ) {
        cout << "Entre com o nome BASE da malha e o nome do output (opcional)" << endl;
        exit(1);
    }

    string basename, ELEfile, NODEfile, OPfile;
    OPfile = "converted_mesh.dat";

    basename = argv[1];
    ELEfile = basename + ".ele";
    NODEfile = basename + ".node";

    // Verificar se arquivos existem
    ifstream inputNODE(NODEfile), inputELE(ELEfile);

    if (! (inputNODE || inputELE)) {
        cout << "Input não encontrado" << endl;
        exit(1);
    }

    if (argc == 3) {
        OPfile = argv[2];
    }

    // ------------------------------------------------------------------------ //
    
    // Procurar nós com coordenada Z = 0; salva número do nó correspondente 
    
    int nNodes, nDim, itmp, nodeID, nBC;
    double  x, y, z;
    vector<int> NODES;

    inputNODE >> nNodes >> nDim >> itmp >> itmp;

    nBC = 0;
    for (int i=0; i<nNodes; i++ ){
        inputNODE >> nodeID >> x >> y >> z;
        if (z == 0 ) {
            nBC++;
            NODES.push_back(nodeID);
        }
    }
    
    inputNODE.close();


    // ------------------------------------------------------------------------ //

    // Gerar arquivo .DAT (apenas informações geométricas e a respeito da malha)

    
    ofstream output(OPfile);
    inputNODE.open(NODEfile);

    int nEls, nNPT;

    inputELE >> nEls >> nNPT >> itmp;

    inputNODE >> nNodes >> nDim >> itmp >> itmp;

    output << nEls << "\t" << nNodes << "\t" << "1\t" << endl;
    output << "1\t" << "1\t" << "1" << endl;
    output << endl;

    for (int i=0; i<nNodes; i++ ){
        inputNODE >> nodeID >> x >> y >> z;
        output << scientific << setprecision(15) << setw(25) << x <<
        setw(25) << y <<
        setw(25) << z << endl;
    }

    inputNODE.close();
    
    int n1, n2, n3, n4;
    for (int i=0; i<nEls; i++ ) {
        inputELE >> itmp >> n1 >> n2 >> n3 >> n4;
        output << setw(15) << n1 <<
            setw(15) << n2 <<
            setw(15) << n3 <<
            setw(15) << n4 << endl;
    }

    inputELE.close();

    output << nBC << endl;;
    for (int i=0; i<nBC; i++){
        output << NODES[i] << "    0   0   0" << endl;
    }

    output.close();

    cout << "Número de nós das condições de contorno: " << nBC << endl;
    return 0;
}
