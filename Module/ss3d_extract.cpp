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
#include <vector>
#include <cmath>
#include <unordered_map>
#include "../Algorithm/ArrayList.hpp"
#include "../Algorithm/AATools.hpp"
#include "../Algorithm/Math.hpp"
#include "../Class/Structure.hpp"
#include "../Class/Bond.hpp"
#include "../Class/Interaction.hpp"
#include "ss3d_extract.hpp"

using namespace std;

ArrayList<Interaction> ss3d::extract(Structure &frame, ArrayList<Bond> &bonds, Param_ext &param) {


    //Variables
    AATools aaTools;
    ArrayList<Interaction> inter;
    unordered_map<unsigned int, vector<unsigned int> > aligned_calpha;

    //Reading frame
    for(unsigned int k = 0; k < frame.molecules.length(); k++){

        //Checking if the molecule is a residue
        if(frame.molecules[k].getType()=='O')
            break;

        unsigned int CA = AATools::findAtom('C',frame.molecules[k]);
        int x = frame.molecules[k].atoms[CA].x;
        int y = frame.molecules[k].atoms[CA].y;
        int z = frame.molecules[k].atoms[CA].z;

        double distanceD, distanceA;
        unsigned int CAd, CAa;


        for(unsigned int i = 0; i < bonds.length(); i++){

            //Getting CA index
            CAd = AATools::findAtom('C',bonds[i],true);
            CAa = AATools::findAtom('C',bonds[i],false);

            //Calculating distance CA-CA (donor)
            distanceD=pow((x - bonds[i].atomsD[CAd].x),2);
            distanceD+=pow((y - bonds[i].atomsD[CAd].y),2);
            distanceD+=pow((z - bonds[i].atomsD[CAd].z),2);
            distanceD=sqrt(distanceD);

            //Calculating distance CA-CA (acceptor)
            distanceA=pow((x - bonds[i].atomsA[CAa].x),2);
            distanceA+=pow((y - bonds[i].atomsA[CAa].y),2);
            distanceA+=pow((z - bonds[i].atomsA[CAa].z),2);
            distanceA=sqrt(distanceA);

            //If distance <= to the cutoff
            if(distanceA<=param.rad || distanceD<=param.rad) {

                aligned_calpha[i].push_back(k);
            }

        }
    }


    //Creating interaction object
    for(unsigned int i = 0; i < bonds.length(); i++){

        //Adicionando informações sobre a ligação
        inter.add(bonds[i].getIndexD(), bonds[i].getIndexA());
        inter.last().setTypeD(bonds[i].getTypeD());
        inter.last().setTypeA(bonds[i].getTypeA());
        inter.last().setFrequency(bonds[i].getFrequency());

        for(auto &mol : aligned_calpha[i]){

            inter.last().molecules.add(frame.molecules[mol].getIndex());
            inter.last().molecules.last().setType(frame.molecules[mol].getType());
        }
    }

    return inter;
}
