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
#include <fstream>
#include "HandleInput.hpp"

#define vecStr vector<string>

using namespace std;

//Checks if the parameter is expected, and if it contains all obligatory parameters
bool HandleInput::checkInput(vecStr args, const string flag, const vecStr inputObligatory, const vecStr inputOptional, const unsigned int nOpt){

    if(args.size()>=inputObligatory.size()+nOpt){

        //Variables
        unsigned int cntO = 0;
        unsigned int cntV = 0;
        unsigned int pos = 0;
        bool check=true;
        bool check2;
        bool found;
        string optP;

        //Checks whether all mandatory parameters are present
        for(auto &obl : inputObligatory){

            //Checking if this parameter requires a sub-parameter (e.g path, number, etc)
            if(obl.find("?")!=string::npos){
                check2 = true;
                optP = obl.substr(0,obl.find("?"));
            }else{
                check2 = false;
                optP = obl;
            }
            pos = 0;
            for(auto &arg : args){
                pos++;
                if(arg==optP){
                    arg="\0";
                    cntO++;
                    break;
                }
            }
            if(check2 && args.size()>=pos+1){
                args[pos]="\0";
            }else{
                check = false;
            }

        }

        //Checks whether all optional parameters are present
        for(auto &opt : inputOptional){

            //Checking if this parameter requires a sub-parameter (e.g path, number, etc)
            if(opt.find("?")!=string::npos){
                check2 = true;
                optP = opt.substr(0,opt.find("?"));
            }else{
                check2 = false;
                optP = opt;
            }
            pos = 0;
            found = false;
            for(auto &arg : args){
                pos++;
                if(arg==optP){
                    arg="\0";
                    found = true;
                    cntV++;
                    break;
                }
            }
            if(check2 && found && args.size()>=pos+1){
                args[pos]="\0";
            }else if(check2 && found && args.size()<pos+1){
                check = false;
                break;
            }

        }

        //Checking for unexpected parameter
        for(auto &arg : args){
            if(arg!="\0"){
                cout << "Error in user input: '" << arg << "' is not expected." << endl;
                cout << "For more informations, try the 'ss3d '" << flag << " command." << endl;
                return false;
            }
        }

        //Showing error message
        if(!check || cntV<nOpt || cntO<inputObligatory.size()){
            cout << "Error in user input: Missing obligatory parameter." << endl;
            cout << "For more informations, try the 'ss3d " << flag << "' command." << endl;
            return false;
        }

        return true;

    //Program exec. without any parameter
    }else{

        cout << "Error in user input: Missing obligatory parameter." << endl;
        cout << "For more informations, try the 'ss3d " << flag << "' command." << endl;
        return false;
    }

}

//Returns a sub-parameter (path, number, etc)
string HandleInput::getParameter(vecStr args, string value){

    unsigned int pos = 0;
    bool found = false;
    for(unsigned int i = 0; i < args.size(); i++){

        if(args[i]==value){
            found = true;
            pos=i+1;
            break;
        }
    }

    if(found)
        return args[pos];
    else
        return args[0];
}

//Checks if a optional parameter exists
bool HandleInput::checkParameter(vecStr args, string value){

    for(unsigned int i = 0; i < args.size(); i++){

        if(args[i]==value){
            return true;
        }
    }

    return false;
}

//Checks if input file is accessible
bool HandleInput::checkIn(string path, string ext){

    //Checks for correct extension
    if(path.length()<=ext.length() || path.find(ext)==0 || path.find(ext)==string::npos){

        cout << "ERROR: Input file must be a " << ext << endl;
        return false;
    }

    ifstream fin(path.c_str());
    if(!fin.is_open()){

        cout << "ERROR: Could not open '" << path << "' file for reading. Is the path correct?" << endl;
        return false;
    }

    return true;
}

//Checks if output file is accessible
bool HandleInput::checkOut(string path, string ext){

    if(path.length()<=ext.length() || path.find(ext)==0 || path.find(ext)==string::npos){

        path+=ext;
        cout << "Output will be '" << path << "'." << endl << endl;
    }

    ifstream fin(path.c_str());

    //Checks if output is exists, and if yes, create a backup file
    if(fin.is_open()){
        fin.close();
        cout << "BACKUP: Output exist! Creating backup '" << path << "_backup'." << endl;
        string str="mv "+path+" "+path+"_backup";

        if(system(str.c_str())){
            cout << "ERROR: Could not move file." << endl;
            return false;
        }
    }else
        fin.close();

    ofstream fout(path.c_str());
    if(!fout.is_open()){
        cout << "ERROR: Could not create output '" << path << "' file. Is the path correct?" << endl;
        return false;
    }

    fout.close();

    return true;
}

