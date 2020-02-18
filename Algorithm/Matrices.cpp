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
#include <string>
#include "../Class/Molecule.hpp"
#include "../Class/Matrices.hpp"

using namespace std;

int Matrix::align(string &mat, Molecule &molA, Molecule &molB){

    unsigned short int x, y;
    char type1 = molA.getType(), type2 = molB.getType();

    //PAM250 or BLOSUM62
    if(mat.find("pam250")!=string::npos || mat.find("blosum62")!=string::npos){

        if(type1=='A'){
            x=0;
        }else if(type1=='R'){
            x=1;
        }else if(type1=='N' || type1=='B'){
            x=2;
        }else if(type1=='D'){
            x=3;
        }else if(type1=='C' || type1=='X'){
            x=4;
        }else if(type1=='Q' || type1=='Z'){
            x=5;
        }else if(type1=='E'){
            x=6;
        }else if(type1=='G'){
            x=7;
        }else if(type1=='H'){
            x=8;
        }else if(type1=='I'){
            x=9;
        }else if(type1=='L'){
            x=10;
        }else if(type1=='K'){
            x=11;
        }else if(type1=='M'){
            x=12;
        }else if(type1=='F'){
            x=13;
        }else if(type1=='P'){
            x=14;
        }else if(type1=='S'){
            x=15;
        }else if(type1=='T'){
            x=16;
        }else if(type1=='W'){
            x=17;
        }else if(type1=='Y'){
            x=18;
        }else if(type1=='V'){
            x=19;
        }

        if(type2=='A'){
            y=0;
        }else if(type2=='R'){
            y=1;
        }else if(type2=='N' || type2=='B'){
            y=2;
        }else if(type2=='D'){
            y=3;
        }else if(type2=='C' || type2=='X'){
            y=4;
        }else if(type2=='Q' || type2=='Z'){
            y=5;
        }else if(type2=='E'){
            y=6;
        }else if(type2=='G'){
            y=7;
        }else if(type2=='H'){
            y=8;
        }else if(type2=='I'){
            y=9;
        }else if(type2=='L'){
            y=10;
        }else if(type2=='K'){
            y=11;
        }else if(type2=='M'){
            y=12;
        }else if(type2=='F'){
            y=13;
        }else if(type2=='P'){
            y=14;
        }else if(type2=='S'){
            y=15;
        }else if(type2=='T'){
            y=16;
        }else if(type2=='W'){
            y=17;
        }else if(type2=='Y'){
            y=18;
        }else if(type2=='V'){
            y=19;
        }

        if(mat.find("pam250")!=string::npos)
            return pam250[y][x];
        else
            return blosum62[y][x];
    }

}
