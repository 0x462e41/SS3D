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
#include "Atom.hpp"

using namespace std;

class Bond{

    private:

    unsigned short int indexD, indexA, frequency = 0;
    char typeD, typeA;

    public:

    Bond(){}

    Bond(unsigned short int indexD, unsigned short int indexA) : indexD(indexD), indexA(indexA) {}

    virtual ~Bond(){}

    ArrayList<Atom> atomsD, atomsA;

    //Getters & Setters
    auto getIndexD(){
        return this->indexD;
    }

    auto getIndexA(){
        return this->indexA;
    }

    void setIndexD(unsigned short indexD){
        this->indexD=indexD;
    }

    void setIndexA(unsigned short indexA){
        this->indexA=indexA;
    }

    char getTypeD(){
        return this->typeD;
    }

    void setTypeD(char typeD){
        this->typeD=typeD;
    }

    char getTypeA(){
        return this->typeA;
    }

    void setTypeA(char typeA){
        this->typeA=typeA;
    }

    auto getFrequency(){
        return this->frequency;
    }

    void setFrequency(unsigned short frequency){
        this->frequency=frequency;
    }

    void incFrequency(){
        this->frequency++;
    }

};
