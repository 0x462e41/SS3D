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

#include <iostream>

template<class CLASS>
class ArrayList{

    private:
    CLASS **Objetos;
    unsigned int leng;
    unsigned int cap;
    bool _cleanable;

    public:

    ArrayList(){
        this->_cleanable=true;
        this->Objetos=0;
        this->leng=0;
        this->cap=1;
        this->Objetos = new CLASS*[cap];
    }

    ArrayList(unsigned int cap){
        this->_cleanable=true;
        this->Objetos=0;
        this->leng=cap;
        this->cap=cap;
        this->Objetos = new CLASS*[this->cap];
        for(unsigned int i=0;i<cap;i++)
            this->Objetos[i]=nullptr;
    }

    ~ArrayList(){
        if(this->_cleanable){
            for(int i=0;i<this->leng;i++)
                delete this->Objetos[i];
            delete [] this->Objetos;
        }else{
            delete [] this->Objetos;
        }
    }

    unsigned int length(){
        return this->leng;
    }

    unsigned int capacity(){
        return this->cap;
    }

    void add(CLASS *cls){
        if(this->leng==this->cap){
            this->cap++;
            this->cap*=2;
            CLASS **temp = new CLASS*[cap];
            for(int i=0;i<this->leng;i++)
                temp[i]=this->Objetos[i];
            delete [] this->Objetos;
            this->Objetos=temp;
        }

        this->Objetos[leng] = cls;
        this->leng++;
    }

    template <typename... TArgs>
    void add(TArgs&&... MArgs){
        if(this->leng==this->cap){
            this->cap++;
            this->cap*=2;
            CLASS **temp = new CLASS*[cap];
            for(int i=0;i<this->leng;i++)
                temp[i]=this->Objetos[i];
            delete [] this->Objetos;
            this->Objetos=temp;
        }

        this->Objetos[leng] = new CLASS(std::forward<TArgs>(MArgs)...);
        this->leng++;
    }

    void remove(unsigned int indice){
        if((this->leng>0)&&(indice<leng)){
            unsigned int cont=0;
            CLASS **temp = new CLASS*[this->cap];
            for(unsigned int i=0;i<this->leng;i++)
                if(i!=indice){
                    temp[cont]=this->Objetos[i];
                    cont++;
                }else if(this->_cleanable)
                    delete this->Objetos[i];

            delete [] this->Objetos;
            this->Objetos=temp;
            this->leng--;
        }
    }

    void remove(CLASS *cls){
        for(unsigned int i=0;i<this->leng;i++){
            if(this->Objetos[i] == cls){
                unsigned int cont=0;
                CLASS **temp = new CLASS*[this->cap];
                for(unsigned int j=0;j<this->leng;j++)
                    if(j!=i){
                        temp[cont]=this->Objetos[j];
                        cont++;
                    }else if(this->_cleanable)
                        delete this->Objetos[i];

                delete [] this->Objetos;
                this->Objetos=temp;
                this->leng--;
                break;
            }
        }
    }


    CLASS *getPointer(unsigned int i){
        if(i<this->leng && i>=0)
            return this->Objetos[i];
        else
            return nullptr;
    }

    CLASS &operator[](unsigned int i){
        if(i>=0 && i<this->leng)
            return *this->Objetos[i];
    }

    CLASS& last(){

        return *this->Objetos[this->leng-1];
    }

    CLASS *&getPointerReference(unsigned int i){
        if(i<this->leng && i>=0)
            return this->Objetos[i];
    }

    void clear(){

        if(this->_cleanable)
            for(int i=0;i<this->leng;i++)
                delete this->Objetos[i];
        delete [] this->Objetos;

        this->Objetos=nullptr;
        this->leng=0;
        this->cap=1;
        this->Objetos = new CLASS*[this->cap];
    }

    void destroy(){

        if(this->_cleanable)
            for(int i=0;i<this->leng;i++)
                delete this->Objetos[i];
        delete [] this->Objetos;

        this->Objetos=nullptr;
        this->leng=0;
        this->cap=0;
    }

    void setCleanable(bool clean){
        this->_cleanable=clean;
    }

    template <typename TArg>
    void runLambda(TArg functor){

        for(unsigned int i = 0; i < this->leng; i++){

            functor(*this->Objetos[i]);
        }
    }

};
