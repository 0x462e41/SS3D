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

#pragma once

//Include files
#include <iostream>
#include <string>
#include "Molecule.hpp"

using namespace std;

namespace Matrix{


    //Matrix PAM 250
    const static float pam250[20][20]=
    {
    {2,-2,0,0,-2,0,0,1,-1,-1,-2,-1,-1,-4,1,1,1,-6,-3,0}, //A
    {-2,6,0,-1,-4,1,-1,-3,2,-2,-3,3,0,-4,0,0,-1,2,-4,-2}, //R
    {0,0,2,2,-4,1,1,0,2,-2,-3,1,-2,-4,-1,1,0,-4,-2,-2}, //N
    {0,-1,2,4,-5,2,3,1,1,-2,-4,0,-3,-6,-1,0,0,-7,-4,-2}, //D
    {-2,-4,-4,-5,4,-5,-5,-3,-3,-2,-6,-5,-5,-4,-3,0,-2,-8,0,-2}, //C
    {0,1,1,2,-5,4,2,-1,3,-2,-2,1,-1,-5,0,-1,-1,-5,-4,-2}, //Q
    {0,-1,1,3,-5,2,4,0,1,-2,-3,0,-2,-5,-1,0,0,-7,-4,-2}, //E
    {1,-3,0,1,-3,-1,0,5,-2,-3,-4,-2,-3,-5,-1,1,0,-7,-5,-1}, //G
    {-1,2,2,1,-3,3,1,-2,6,-2,-2,0,-2,-2,0,-1,-1,-3,0,-2}, //H
    {-1,-2,-2,-2,-2,-2,-2,-3,-2,5,2,-2,2,1,-2,-1,0,-5,-1,4}, //I
    {-2,-3,-3,-4,-6,-2,-3,-4,-2,2,6,-3,4,2,-3,-3,-2,-2,-1,2}, //L
    {-1,3,1,0,-5,1,0,-2,0,-2,-3,5,0,-5,-1,0,0,-3,-4,-2}, //K
    {-1,0,-2,-3,-5,-1,-2,-3,-2,2,4,0,6,0,-2,-2,-1,-4,-2,2}, //M
    {-4,-4,-4,-6,-4,-5,-5,-5,-2,1,2,-5,0,9,-5,-3,-2,0,7,-1}, //F
    {1,0,-1,-1,-3,0,-1,-1,0,-2,-3,-1,-2,-5,6,1,0,-6,-5,-1}, //P
    {1,0,1,0,0,-1,0,1,-1,-1,-3,0,-2,-3,1,2,1,-2,-3,-1}, //S
    {1,-1,0,0,-2,-1,0,0,-1,0,-2,0,-1,-2,0,1,3,-5,-3,0}, //T
    {-6,2,-4,-7,-8,-5,-7,-7,-3,-5,-2,-3,-4,0,-6,-2,-5,17,0,-6}, //W
    {-3,-4,-2,-4,0,-4,-4,-5,0,-1,-1,-4,-2,7,-5,-3,-3,0,10,-2}, //Y
    {0,-2,-2,-2,-2,-2,-2,-1,-2,4,2,-2,2,-1,-1,-1,0,-6,-2,4}, //V
    };

    //Matrix BLOSUM62
    const static float blosum62[20][20]=
    {
    {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0}, //A
    {-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3}, //R
    {-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3}, //N
    {-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3}, //D
    {0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1}, //C
    {-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2}, //Q
    {-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2}, //E
    {0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3}, //G
    {-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3}, //H
    {-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3}, //I
    {-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1}, //L
    {-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2}, //K
    {-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1}, //M
    {-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1}, //F
    {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2}, //P
    {1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2}, //S
    {0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0}, //T
    {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3}, //W
    {-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1}, //Y
    {0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}, //V
    };

    float align(string &mat, Molecule &molA, Molecule &molB);

};