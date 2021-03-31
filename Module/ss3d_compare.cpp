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
#include <iomanip>
#include <algorithm>
#include <map>
#include "../Algorithm/ArrayList.hpp"
#include "../Algorithm/Math.hpp"
#include "../Class/Interaction.hpp"
#include "../Class/Matrices.hpp"
#include "../Algorithm/AATools.hpp"
#include "ss3d_compare.hpp"

//Define consts
#define ATOM "ATOM"
#define ENDF "END"


using namespace std;

//Calculates the score using the chosen matrix
void getScore(ArrayList<Interaction> &protA, ArrayList<Interaction> &protB,
        ss3d::Param_comp &param, vector<unsigned short int> &resA, vector<unsigned short int> &resB,
        vector<int> &score, bool flag, ofstream &fout){

    //Variable
    int tmp_s, tmp_sum;
    unsigned int cnt;
    vector<char> itA, itB;
    vector<int> scs;

    for(unsigned int i=0; i<protA.length(); i++){ //Loop through all 'A' interactions
        for(unsigned int j=0; j<protB.length(); j++){ //Loop through all 'B' interactions

            //Pairing ...
            if(protA[i].getIndexD()==protB[j].getIndexD() && protA[i].getIndexA()==protB[j].getIndexA()){

                tmp_sum=0; //score value;
                cnt=0; //Number of common residues around the contact

                //If raw flag is set, print raw informations about
                if(flag) {
                    fout << protA[i].getIndexD() << " " << protA[i].getIndexA() << "\n";
                    itA.clear();
                    itB.clear();
                    scs.clear();
                }

                //Lining up the interactions, and calculating the score
                for(unsigned int k=0; k<protA[i].molecules.length(); k++) //Loop through all residues around A
                    for(unsigned int l=0; l<protB[j].molecules.length(); l++) //Loop through all residues around B
                        if(protA[i].molecules[k].getIndex()==protB[j].molecules[l].getIndex()){
                            cnt++;
                            tmp_s=Matrix::align(param.mat, protA[i].molecules[k], protB[j].molecules[l]);
                            tmp_sum+=tmp_s;
                            if(flag) {
                                itA.push_back(protA[i].molecules[k].getType());
                                itB.push_back(protB[j].molecules[l].getType());
                                scs.push_back(tmp_s);
                            }
                        }

                //Checking if the minimum number of common residues around this contact was attended
                if(cnt>=param.min){
                    resA.push_back(protA[i].getIndexD());
                    resB.push_back(protA[i].getIndexA());
                    score.push_back(tmp_sum);
                    if(param.dup) {
                        resA.push_back(protA[i].getIndexA());
                        resB.push_back(protA[i].getIndexD());
                        score.push_back(tmp_sum);
                    }
                }

                if(flag){
                    for(auto const &x : itA) {
                        fout << x;
                    }
                    fout << "\n";
                    for(auto const &x : itB) {
                        fout << x;
                    }
                    fout << "\n|";
                    for(auto const &x : scs) {
                        fout << x << "|";
                    }
                    fout << "\n";
                    if(cnt>=param.min)
                        fout << "Sum: " << score.back() << "\n";
                    else
                        fout << "Sum: skipped" << "\n";
                    fout << "===============\n";
                }

                break;

            }
        }
    }
}

//Generating a matrix output
bool ss3d::compare(ArrayList<Interaction> &protA, ArrayList<Interaction> &protB, Param_comp &param) {

    //Variable
    ofstream fout(param.path);
    vector<unsigned short int> resA, resB;
    vector<int> score;


    //Calculates the score using the chosen matrix
    getScore(protA, protB, param, resA, resB, score, param.raw, fout);

    //Outputting result
    if(!param.raw){
        fout << "#IndexA IndexB Score\n";
        if(param.norm) {
            auto max = *max_element(score.begin(), score.end());
            for(int i=0; i < score.size(); i++) {
                fout << resA[i] << " " << resB[i] << " " << static_cast<float>(score[i])/static_cast<float>(max) << "\n";
            }
        } else {
            for(int i=0; i < score.size(); i++) {
                fout << resA[i] << " " << resB[i] << " " << score[i] << "\n";
            }
        }
    }

    fout.close();
    cout << endl << "Done!" << endl;


    return true;
}

//Mapping to bfactor
bool ss3d::map_bfactor(ArrayList<Interaction> &protA, ArrayList<Interaction> &protB, Param_comp &param) {

    //Variable
    ofstream fout(param.path);
    fout << setprecision(3);
    fout << fixed;
    vector<unsigned short int> resA, resB;
    vector<int> score;
    unsigned int nAtom = 1;
    unsigned short int lastIndex = 0xFFFF;
    AATools aaTools;
    map<unsigned short int, float> resid_score;
    map<unsigned short int, unsigned int> resid_cnt;
    pair<map<unsigned short int, float>::iterator,bool> ret;

    //Calculates the score using the chosen matrix
    param.dup=false;
    getScore(protA, protB, param, resA, resB, score, param.raw, fout);

    //Mapping sum
    for(unsigned int i = 0; i < score.size(); i++){
        ret = resid_score.insert(pair<unsigned short int, float>(resA[i],static_cast<float>(score[i])));
        if(!ret.second) {
            ret.first->second+=score[i];
            resid_cnt[resA[i]]+=1;
        } else {
            resid_cnt[resA[i]]=1;
        }

        ret = resid_score.insert(pair<unsigned short int, float>(resB[i],static_cast<float>(score[i])));
        if(!ret.second) {
            ret.first->second+=score[i];
            resid_cnt[resB[i]]+=1;
        } else {
            resid_cnt[resB[i]]=1;
        }
    }

    //Average
    for (auto &x : resid_score) {
        x.second/=static_cast<float>(resid_cnt[x.first]);
    }

    //Normalizes the result (if param exists)
    if(param.norm) {
        auto max_e = max_element(resid_score.begin(), resid_score.end(),
                   [](const pair<unsigned short int, float> &p1,
                   const pair<unsigned short int, float> &p2) {
                       return p1.second < p2.second;
                   }
        );

        auto min_e = min_element(resid_score.begin(), resid_score.end(),
                   [](const pair<unsigned short int, float> &p1,
                   const pair<unsigned short int, float> &p2) {
                       return p1.second < p2.second;
                   }
        );

        auto max = max_e->second;
	auto min = min_e->second;
        for (auto &x : resid_score) {
	    x.second=(x.second-min)/(max-min);
        }
    }

    //Outputting result

    //Printing header
    fout << "REMARK    GENERATED BY SS3D" << "\n";

    //Printing molecules
    param.frame->molecules.runLambda([&fout,&nAtom,&lastIndex,&aaTools,&resid_score](auto &mol){

    if(mol.getIndex()!=lastIndex){
        lastIndex=mol.getIndex();
        aaTools.setMolecule(mol.getType());
    }
    mol.atoms.runLambda([&mol,&fout,&nAtom,&aaTools,&resid_score](auto &at){

        if(resid_score.find(mol.getIndex()) != resid_score.end()) {
            fout << ATOM;
            fout << setw(7) << setfill(' ') << nAtom++;
            fout << " ";
            fout << aaTools.getAtom(at.atom);
            fout << AATools::convertCode(mol.getType());
            fout << setw(6) << setfill(' ') << mol.getIndex();
            fout << "     ";
            fout << setw(7) << setfill(' ') << at.x/1000.f << " ";
            fout << setw(7) << setfill(' ') << at.y/1000.f << " ";
            fout << setw(7) << setfill(' ') << at.z/1000.f << "  1.00 ";
            fout << resid_score[mol.getIndex()] << "\n";
        }
    });
    });

    //Printing footnote
    fout << "TER" << "\n";
    fout << ENDF << "\n";

    fout.close();
    cout << endl << "Done!" << endl;

    return true;
}
