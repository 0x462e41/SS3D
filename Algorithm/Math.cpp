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
#include "Math.hpp"

using namespace std;

bool Math::is_a_number(string val){

    bool result = true;

    for(unsigned int i = 0; i < val.length(); i++){

        if(!isdigit(val[i])){
            result = false;
            break;
        }
    }

    if(!result)
        cout << "ERROR: Expecting a number for parameter. Aborting ..." << endl;

    return result;
}
