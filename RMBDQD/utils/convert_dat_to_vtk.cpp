#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<vector>

void cvt_dat_to_vkt(std::string& DATnfile, std::string& OPnfile);

using namespace std;

int main(int argc, const char *argv[])
{
    if (argc<2) exit(1);

    string DATnfile, OPnfile;

    DATnfile = argv[1];
    OPnfile = "generated_vtk.vtk";

    if (argc == 3) OPnfile = argv[2];

    /*-------------------------------------------------------------------------------*/

    cvt_dat_to_vkt(DATnfile, OPnfile);
    
    return 0;
}

void cvt_dat_to_vkt(std::string& DATnfile, std::string& OPnfile){

    std::ifstream input;
    std::ofstream output;

    double tmp;
    int itmp;
    int nEls, nNodes, nBC;

    std::vector<std::vector<double> > nodes;
    std::vector<std::vector<int> > conect;
    std::vector<int> BC;

    input.open(DATnfile);

    input >> nEls >> nNodes >> tmp >> tmp >> tmp >> tmp;

    /* Lê os nós */

    for (int i=0; i<nNodes; i++){
        std::vector<double> VDtmp;
        for (int j=0; j<3; j++){
            input >> tmp;
            VDtmp.push_back(tmp);
        }
        nodes.push_back(VDtmp);
    }

    /* Lê conectividade */

    for (int i=0; i<nEls; i++){
        std::vector<int> VItmp;
        for (int j=0; j<4; j++){
            input >> itmp;
            VItmp.push_back(itmp);
        }
        conect.push_back(VItmp);
    }

    input >> nBC;

    for (int i=0; i<nBC; i++){
        input >> itmp >> tmp >> tmp >> tmp;
        BC.push_back(itmp);
    }

    input.close();

    /*-------------------------------------------------------------------------------*/

    output.open(OPnfile);

    output << "# vtk DataFile Version 2.0\n" << "Converted from DAT file\n";
    output << "ASCII\n";
    output << "DATASET UNSTRUCTURED_GRID\n";
    output << "POINTS " << nNodes << " float\n";

    for (int i=0; i<nNodes; i++){
        output << std::fixed << std::setprecision(12) <<
            std::setw(24) << nodes[i][0] <<
            std::setw(24) << nodes[i][1] <<
            std::setw(24) << nodes[i][2] << std::endl;
    }

    output << "CELLS " << nEls << " " << 5*nEls << "\n";

    for (int i=0; i<nEls; i++){
        output << "4 " << 
            std::setw(20) << conect[i][0] <<
            std::setw(20) << conect[i][1] <<
            std::setw(20) << conect[i][2] <<
            std::setw(20) << conect[i][3] << std::endl;
    }

    output << "CELL_TYPES " << nEls << " \n";

    for (int i=0; i<nEls; i++){
        output << "10\n";
    }

    output << "POINT_DATA " << nNodes << "\n";
    output << "SCALARS aaa float\n";
    output << "LOOKUP_TABLE default\n";

    /* Escreve as condições de contorno com outros valores,
       útil na visualização  */

    bool is_in_arr = false;

    for (int i=0; i<nNodes; i++){
        for (int j=0; j<nBC; j++){
            if (BC[j] == i+1){
                is_in_arr = true;
            }
        }

        if (is_in_arr) {
            output << "1\n";
        }
        else{
            output << "0\n";
        }

        is_in_arr = false;
    }

    output.close();
}
