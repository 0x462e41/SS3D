/*
    A bioinformatics tool, named Sequence Similarity 3D (SS3D)
    Copyright (C) 2020  Igor Daniel M. Lima <igor.sj13@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

//Include Files
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../Algorithm/ArrayList.hpp"
#include "../Algorithm/Math.hpp"
#include "../Class/Interaction.hpp"
#include "../Class/Matrices.hpp"
#include "ss3d_compare.hpp"

using namespace std;

bool ss3d::compare(ArrayList<Interaction> &protA, ArrayList<Interaction> &protB, Param_comp &param) {


    //Variable
    ofstream fout(param.path);
    float conscore;
    unsigned int cnt;

    //Calculates the conscore using the chosen matrix, and saves it into the file
    fout << "#IndexA IndexB Conscore\n";
    for(unsigned int i=0; i<protA.length(); i++){ //Loop through all 'A' interactions
        for(unsigned int j=0; j<protB.length(); j++){ //Loop through all 'B' interactions

            //Pairing ...
            if(protA[i].getIndexD()==protB[j].getIndexD() && protA[i].getIndexA()==protB[j].getIndexA()){

                conscore=0; //Conscore value;
                cnt=0; //Number of common residues around the contact

                //Lining up the interactions, and calculating the conscore
                for(unsigned int k=0; k<protA[i].molecules.length(); k++) //Loop through all residues around A
                    for(unsigned int l=0; l<protB[j].molecules.length(); l++) //Loop through all residues around B
                        if(protA[i].molecules[k].getIndex()==protB[j].molecules[l].getIndex()){
                            cnt++;
                            conscore+=Matrix::align(param.mat, protA[i].molecules[k], protB[j].molecules[l]);
                        }

                //Checking if the minimum number of common residues around this contact was attended
                if(cnt>=param.min){
                    fout << protA[i].getIndexD() << " " << protA[i].getIndexA() << " " << conscore << " " << "\n";
                    if(param.dup)
                        fout << protA[i].getIndexA() << " " << protA[i].getIndexD() << " " << conscore << " " << "\n";
                }

                break;

            }
        }
    }

    fout.close();
    cout << endl << "Done!" << endl;


    return true;
}
