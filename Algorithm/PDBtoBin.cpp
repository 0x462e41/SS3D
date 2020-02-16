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
#include "../Class/Structure.hpp"
#include "AATools.hpp"
#include "PDBtoBin.hpp"

//Const markers
#define ATOM "ATOM"
#define ENDF "END"

using namespace std;

Structure PDBtoBin::run(string pathPDB){

    //Variables
    ifstream fin(pathPDB);
    string str, subStr;
    string::size_type pos;
    unsigned short int index, lastIndex=0xFFFF;
    bool header = true;
    char code, type;
    AATools aaTools;
    Structure frame;

    //Creating binary
    while(getline(fin,str)){

        //Saving atom
        if(str.find(ATOM)!=string::npos){

            //Checking if the molecule is a residue
            subStr=str.substr(17,3);
            code = AATools::convertCode(subStr);

            if(code!=0x00){

                //Index
                subStr=str.substr(22,4);
                pos=subStr.find_last_of(" ");
                subStr=subStr.substr(pos+1);
                index = stoi(subStr);

                if(lastIndex!=index){
                    frame.molecules.add(index);
                    lastIndex=index;
                    frame.molecules.last().setType(code);
                    aaTools.setMolecule(code);
                }

                type = aaTools.atomIndex(str.substr(12,5));

                if(type != 100){

                    //Type and position
                    frame.molecules.last().atoms.add(static_cast<long>(stod(str.substr(31,7))*1000),
                            static_cast<long>(stod(str.substr(39,7))*1000),
                            static_cast<long>(stod(str.substr(47,7))*1000));
                    frame.molecules.last().atoms.last().atom=type;
                }

            }
        }

        //END frame
        else if(str.find(ENDF)!=string::npos){

            break;
        }

    }

    fin.close();

    return frame;
}
