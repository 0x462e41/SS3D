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

/*
    Note: This tool has been extracted from larger software. Unfortunately,
    there are some small redundancies and small pieces of inaccessible code,
    which would be a challenge to remove. But in my tests, the compiler
    optimization at level 2 (-O2) fixed it and the compiled code works very
    satisfactorily.
*/

//Constant
#define SKIP 3 //minimum distance between residues in the primary structure to be considered
#define RADIUS 10 //maximum radius to search for interaction between alpha carbons, in Ångströms
#define DISTANCE 10 //maximum distance to be considered a contact between alpha carbons, in Ångströms
#define MIN 1 //minimum number of common residues around the contact
#define MAT "blosum62" // matrix to be used to evaluate the conscore (it can be "blosum62" or "pam250")


//Include Files
#include <iostream>
#include <vector>
#include "Module/cmap.hpp"
#include "Module/ss3d_extract.hpp"
#include "Module/ss3d_compare.hpp"
#include "Algorithm/HandleInput.hpp"
#include "Algorithm/PDBtoBin.hpp"
#include "Algorithm/Math.hpp"

using namespace std;

void kill_program() {
    cout << endl << "------------------------------------------------" << endl;
    cout << endl << "SS3D has exited with an error! Please read the error message above.";
    cout << endl << endl;
    exit(-1);
}

int main(int argc, char **argv){

    std::ios::sync_with_stdio(false);

    //Variables
    vector<string> args;
    cmap::Param param_cmap;
    ss3d::Param_ext param_ss3d_extract;
    ss3d::Param_comp param_ss3d_comp;

    //Converting input to vector and printing header
    cout << "\t\t\tSS3D - 2019.11.25" << endl << endl;
    cout << "Copyright (C) 2020 Igor Daniel M. Lima <igor.sj13@gmail.com>" << endl;
    cout << "SS3D is free software; you can redistribute it and/or modify it " << endl;
    cout << "under the terms of the GNU Lesser General Public License " << endl;
    cout << "as published by the Free Software Foundation; either version 3" << endl;
    cout << "of the License, or (at your option) any later version." << endl << endl;
    cout << endl << "------------------------------------------------" << endl << endl;

    cout << "Command line: ss3d";
    for(int i = 1; i < argc; i++){
        cout << " " << argv[i];
        args.push_back(argv[i]);
    }
    cout << endl << endl;

    //Help menu
    if(args.size() >= 1)
        if(args[0]=="help" || args[0]=="-h"){
        cout << "SYNOPSIS" << endl << endl;
        cout << "ss3d [-a [.pdb]] [-b [.pdb]] [-o [.xvg/.pdb]] [-dist <number>] [-rad <number>] [-skip <number>]" << endl;
        cout << "[-min <number>] [-mat <matrix type>] [-dup] [-norm] [-bfac] [-raw]" << endl;

        cout << endl << "OPTIONS:" << endl << endl;
        cout << "Obligatory options to specify input & output files:" << endl << endl;
        cout << " -a\n\tStructural file in .pdb format of the first protein." << endl;
        cout << " -b\n\tStructural file in .pdb format of the second protein." << endl;
        cout << " -o\n\tOutput text file in .xvg or .pdb format." << endl;
        cout << endl << "Other options:" << endl << endl;
        cout << " -dist\n\tThe maximum distance to be considered a contact between two alpha carbons, in Ångströms." << endl;
        cout << "\t(default distance is " << DISTANCE << " Ångströms)" << endl;
        cout << " -rad\n\tThe maximum search radius for identifying alpha carbons around a contact, in Ångströms." << endl;
        cout << "\t(default radius is " << RADIUS << " Ångströms)" << endl;
        cout << " -skip\n\tExclude this number of residues within the search radius upstream and downstream in primary " << endl;
        cout << "\tsequence of the alpha carbon involved in the contact. Used to diminish the contribution of " << endl;
        cout << "\tsequentially proximal residues in the scoring. (default skip is " << SKIP << ", which will exclude " << endl;
        cout << "\t" << SKIP*2 << " residues total: the " << SKIP << " previous and " << SKIP << " subsequent to the one ";
        cout << "involved in the contact)." << endl;
        cout << " -min\n\tMinimum (absolute) number of common residues between the two proteins around a particular ";
        cout << "contact" << endl;
        cout << "\tfor it to be accepted. (default is " << MIN << " residue[s])." << endl;
        cout << " -mat\n\tSelect the substitution matrix to be used to evaluate the score. The default matrix is " << MAT;
        cout << ".\n\tThe options are:" << endl;
        cout << "\tpam250 - The PAM250 matrix." << endl;
        cout << "\tblosum62 - The BLOSUM62 matrix." << endl;
        cout << " -dup\n\tDuplicates the matrix in a mirrored fashion." << endl;
        cout << " -norm\n\tNormalizes the result to be between 0 and 1." << endl;
        cout << " -bfac\n\tMaps the score value to the B-factor of the first protein PDB." << endl;
        cout << " -raw\n\tProduces a raw text file containing the local alignment for each contact, ignoring the" << endl;
        cout << "\tparameters \"-bfac\", \"-dup\" and \"-norm\", if present." << endl;

        return 0;
    }

    //Checking parameter
    if(HandleInput::checkInput(args, "-h", {"-a?","-b?","-o?"},
                {"-skip?","-dist?","-rad?","-min?","-mat?","-dup","-norm","-bfac", "-raw"}, 0)){

        //I/O Checking
        if(!HandleInput::checkIn(HandleInput::getParameter(args, "-a"), ".pdb"))
            kill_program();
        if(!HandleInput::checkIn(HandleInput::getParameter(args, "-b"), ".pdb"))
            kill_program();
        string outputPath=HandleInput::getParameter(args, "-o");

        //First protein

        //Reading PDB file
        auto frame = PDBtoBin::run(HandleInput::getParameter(args, "-a"));
        if(frame.molecules.length() == 0)
            kill_program();

        //Generating CA map

        //Checking optional parameter
        if(HandleInput::checkParameter(args,"-dist")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-dist"))) {
                param_cmap.dist = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-dist"))) *1000;
            } else
                kill_program();
        }else{
            param_cmap.dist = RADIUS*1000;
        }

        if(HandleInput::checkParameter(args,"-skip")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-skip"))) {
                param_cmap.skip = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-skip")));
            } else
                kill_program();
        }else{
            param_cmap.skip = SKIP;
        }

        auto bonds = cmap::extract(frame, param_cmap);

        if(bonds.length() == 0)
            kill_program();


        //Getting residues identities

        //Checking optional parameter
        if(HandleInput::checkParameter(args,"-rad")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-rad"))) {
                param_ss3d_extract.rad = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-rad"))) *1000;
            } else
                kill_program();
        }else{
            param_ss3d_extract.rad = RADIUS *1000;
        }

        auto identity_A = ss3d::extract(frame, bonds, param_ss3d_extract);

        if(identity_A.length() == 0)
            kill_program();


        //Second protein

        //Reading PDB file
        auto frame_b = PDBtoBin::run(HandleInput::getParameter(args, "-b"));
        if(frame_b.molecules.length() == 0)
            kill_program();

        //Generating CA map

        //Checking optional parameter
        if(HandleInput::checkParameter(args,"-dist")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-dist"))) {
                param_cmap.dist = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-dist"))) *1000;
            } else
                kill_program();
        }else{
            param_cmap.dist = RADIUS*1000;
        }

        if(HandleInput::checkParameter(args,"-skip")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-skip"))) {
                param_cmap.skip = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-skip")));
            } else
                kill_program();
        }else{
            param_cmap.skip = SKIP;
        }

        auto bonds_b = cmap::extract(frame_b, param_cmap);

        if(bonds_b.length() == 0)
            kill_program();

        //Getting residues identities

        //Checking optional parameter
        if(HandleInput::checkParameter(args,"-rad")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-rad"))) {
                param_ss3d_extract.rad = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-rad"))) *1000;
            } else
                kill_program();
        }else{
            param_ss3d_extract.rad = RADIUS *1000;
        }

        auto identity_B = ss3d::extract(frame_b, bonds_b, param_ss3d_extract);

        if(identity_A.length() == 0)
            kill_program();


        //Comparing the two protein SS3D

        //Checking optional parameter
        param_ss3d_comp.dup = HandleInput::checkParameter(args,"-dup");
        param_ss3d_comp.norm = HandleInput::checkParameter(args,"-norm");
        param_ss3d_comp.raw = HandleInput::checkParameter(args,"-raw");

        if(HandleInput::checkParameter(args,"-min")) {
            if(Math::is_a_number(HandleInput::getParameter(args, "-min"))) {
                param_ss3d_comp.min = static_cast<unsigned int>(stoul(HandleInput::getParameter(args,"-min")));
            } else
                kill_program();
        }else{
            param_ss3d_comp.min = MIN;
        }

        if(HandleInput::checkParameter(args,"-mat")) {
            param_ss3d_comp.mat = HandleInput::getParameter(args,"-mat");
            //Cheking if the matrix exists
            if(param_ss3d_comp.mat!="pam250" && param_ss3d_comp.mat!="blosum62"){
                cout << "ERROR: Matrix '" << param_ss3d_comp.mat
                    << "' does not exist, check the options in the help menu. Aborting ..." << endl;
                kill_program();
            }
        }else{
            param_ss3d_comp.mat = MAT;
        }

        bool result;
        if(HandleInput::checkParameter(args,"-bfac") && !param_ss3d_comp.raw){
            param_ss3d_comp.frame=&frame;
            if(!HandleInput::checkOut(outputPath, ".pdb"))
                kill_program();
            param_ss3d_comp.path = outputPath;
            result = ss3d::map_bfactor(identity_A, identity_B, param_ss3d_comp);
        } else if(!param_ss3d_comp.raw){
            if(!HandleInput::checkOut(outputPath, ".xvg"))
                kill_program();
            param_ss3d_comp.path = outputPath;
            result = ss3d::compare(identity_A, identity_B, param_ss3d_comp);
        } else {
            if(!HandleInput::checkOut(outputPath, ".txt"))
                kill_program();
            param_ss3d_comp.path = outputPath;
            result = ss3d::compare(identity_A, identity_B, param_ss3d_comp);
        }

        if(!result)
            kill_program();


        cout << endl << "------------------------------------------------" << endl;
        cout << endl << "SS3D has exited successfully!" << endl << endl;;

    }else{
        kill_program();
    }

    return 0;
}
