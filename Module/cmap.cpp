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
#include "cmap.hpp"


using namespace std;

inline unsigned Key(unsigned short i,unsigned short j) {return (unsigned) i << 16 | j;}

ArrayList<Bond> cmap::extract(Structure &frame, Param &param){

    //Variables
    AATools aaTools;
    ArrayList<Bond> bonds;
    unordered_map<unsigned int, unsigned int> map_index;
    double distance;
    unsigned int CA1, CA2;


    //Reading frame and finding contacts
    for(unsigned int i = 0; i < frame.molecules.length(); i++){

        //If the molecule isn't a residue
        if(frame.molecules[i].getType()=='O')
            break;

        CA1 = AATools::findAtom('C',frame.molecules[i]);

        for(unsigned int j = i+param.skip; j < frame.molecules.length()-param.skip; j++){

            if(frame.molecules[j].getType()=='O')
                break;

            //Getting CA 2 index
            CA2 = AATools::findAtom('C',frame.molecules[j]);

            //Calculating distance between them
            distance=pow((frame.molecules[i].atoms[CA1].x - frame.molecules[j].atoms[CA2].x),2);
            distance+=pow((frame.molecules[i].atoms[CA1].y - frame.molecules[j].atoms[CA2].y),2);
            distance+=pow((frame.molecules[i].atoms[CA1].z - frame.molecules[j].atoms[CA2].z),2);
            distance=sqrt(distance);

            //If it's <= to the cutoff
            if(distance<=param.dist) {

                //Check if bond exists already
                unsigned int key = Key(frame.molecules[i].getIndex(),frame.molecules[j].getIndex());

                //New bond
                if(map_index.find(key)==map_index.end()){

                    map_index[key]=bonds.length();

                    //Saving molecule into 'bonds'
                    bonds.add(frame.molecules[i].getIndex(),frame.molecules[j].getIndex());
                    bonds.last().setTypeD(frame.molecules[i].getType());
                    bonds.last().setTypeA(frame.molecules[j].getType());

                    //Saving atoms 1st molecule
                    frame.molecules[i].atoms.runLambda([&bonds](auto &at){
                        bonds.last().atomsD.add(at.atom);
                        bonds.last().atomsD.last().x=at.x;
                        bonds.last().atomsD.last().y=at.y;
                        bonds.last().atomsD.last().z=at.z;
                    });

                    //Saving atoms 2st molecule
                    frame.molecules[j].atoms.runLambda([&bonds](auto &at){
                        bonds.last().atomsA.add(at.atom);
                        bonds.last().atomsA.last().x=at.x;
                        bonds.last().atomsA.last().y=at.y;
                        bonds.last().atomsA.last().z=at.z;
                    });


                //Found bond
                }else{

                    bonds[map_index[key]].incFrequency();
                }

            }

        }

    }

    return bonds;
}
