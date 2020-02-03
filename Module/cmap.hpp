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

//Include Files
#include <iostream>
#include "../Algorithm/ArrayList.hpp"
#include "../Class/Bond.hpp"
#include "../Class/Structure.hpp"

using namespace std;

namespace cmap{

    struct Param {
        unsigned int dist, skip;
    };

    ArrayList<Bond> extract(Structure &frame, Param &param);

};
