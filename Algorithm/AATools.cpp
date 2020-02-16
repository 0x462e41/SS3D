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
#include "../Class/Molecule.hpp"
#include "AATools.hpp"

using namespace std;

//Convert 3-letter code to 1-letter code
char AATools::convertCode(string res){

    if(res=="ALA"){
        return 'A';
    }else if(res=="ARG"){
        return 'R';
    }else if(res=="ASN"){
        return 'N';
    }else if(res=="ASP"){
        return 'D';
    }else if(res=="ASX"){
        return 'B';
    }else if(res=="CYS"){
        return 'C';
    }else if(res=="CYM"){
        return 'X';
    }else if(res=="GLU"){
        return 'E';
    }else if(res=="GLN"){
        return 'Q';
    }else if(res=="GLX"){
        return 'Z';
    }else if(res=="GLY"){
        return 'G';
    }else if(res=="HIS"){
        return 'H';
    }else if(res=="ILE"){
        return 'I';
    }else if(res=="LEU"){
        return 'L';
    }else if(res=="LYS"){
        return 'K';
    }else if(res=="MET"){
        return 'M';
    }else if(res=="PHE"){
        return 'F';
    }else if(res=="PRO"){
        return 'P';
    }else if(res=="SER"){
        return 'S';
    }else if(res=="THR"){
        return 'T';
    }else if(res=="TRP"){
        return 'W';
    }else if(res=="TYR"){
        return 'Y';
    }else if(res=="VAL"){
        return 'V';
    }else if(res=="SOL"){
        return 'O';
    }else{
        return 0x00;
    }

}

//Convert 1-letter code to 3-letter code
string AATools::convertCode(char code){

    switch(code){

        case 'A':
            return "ALA";
            break;
        case 'R':
            return "ARG";
            break;
        case 'N':
            return "ASN";
            break;
        case 'D':
            return "ASP";
            break;
        case 'B':
            return "ASX";
            break;
        case 'C':
            return "CYS";
            break;
        case 'X':
            return "CYM";
            break;
        case 'E':
            return "GLU";
            break;
        case 'Q':
            return "GLN";
            break;
        case 'Z':
            return "GLX";
            break;
        case 'G':
            return "GLY";
            break;
        case 'H':
            return "HIS";
            break;
        case 'I':
            return "ILE";
            break;
        case 'L':
            return "LEU";
            break;
        case 'K':
            return "LYS";
            break;
        case 'M':
            return "MET";
            break;
        case 'F':
            return "PHE";
            break;
        case 'P':
            return "PRO";
            break;
        case 'S':
            return "SER";
            break;
        case 'T':
            return "THR";
            break;
        case 'W':
            return "TRP";
            break;
        case 'Y':
            return "TYR";
            break;
        case 'V':
            return "VAL";
            break;
        case 'O':
            return "SOL";
            break;
        default:
            return "xxx";
    }

}

//Set the molecule
void AATools::setMolecule(char mol){

    this->molecule=mol;
}

//Returns the next atom of the residue
string AATools::getAtom(char nAtom){

    switch(this->molecule){

        //ALA
        case 'A':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " HB3 ";
                    break;
                case 8:
                    return " C   ";
                    break;
                case 9:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //ARG
        case 'R':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " HG1 ";
                    break;
                case 9:
                    return " HG2 ";
                    break;
                case 10:
                    return " CD  ";
                    break;
                case 11:
                    return " HD1 ";
                    break;
                case 12:
                    return " HD2 ";
                    break;
                case 13:
                    return " NE  ";
                    break;
                case 14:
                    return " HE  ";
                    break;
                case 15:
                    return " CZ  ";
                    break;
                case 16:
                    return " NH1 ";
                    break;
                case 17:
                    return "1HH1 ";
                    break;
                case 18:
                    return "2HH1 ";
                    break;
                case 19:
                    return " NH2 ";
                    break;
                case 20:
                    return "1HH2 ";
                    break;
                case 21:
                    return "2HH2 ";
                    break;
                case 22:
                    return " C   ";
                    break;
                case 23:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //ASN
        case 'N':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " OD1 ";
                    break;
                case 9:
                    return " ND2 ";
                    break;
                case 10:
                    return "1HD2 ";
                    break;
                case 11:
                    return "2HD2 ";
                    break;
                case 12:
                    return " C   ";
                    break;
                case 13:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //ASP
        case 'D':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " OD1 ";
                    break;
                case 9:
                    return " OD2 ";
                    break;
                case 10:
                    return " C   ";
                    break;
                case 11:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //CYS
        case 'C':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " SG  ";
                    break;
                case 8:
                    return " HG  ";
                    break;
                case 9:
                    return " C   ";
                    break;
                case 10:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //CYM
        case 'X':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " SG  ";
                    break;
                case 8:
                    return " C   ";
                    break;
                case 9:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //GLU
        case 'E':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " HG1 ";
                    break;
                case 9:
                    return " HG2 ";
                    break;
                case 10:
                    return " CD  ";
                    break;
                case 11:
                    return " OE1 ";
                    break;
                case 12:
                    return " OE2 ";
                    break;
                case 13:
                    return " C   ";
                    break;
                case 14:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //GLN
        case 'Q':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " HG1 ";
                    break;
                case 9:
                    return " HG2 ";
                    break;
                case 10:
                    return " CD  ";
                    break;
                case 11:
                    return " OE1 ";
                    break;
                case 12:
                    return " NE2 ";
                    break;
                case 13:
                    return "1HE2 ";
                    break;
                case 14:
                    return "2HE2 ";
                    break;
                case 15:
                    return " C   ";
                    break;
                case 16:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //GLY
        case 'G':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA1 ";
                    break;
                case 4:
                    return " HA2 ";
                    break;
                case 5:
                    return " C   ";
                    break;
                case 6:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //HIS
        case 'H':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " ND1 ";
                    break;
                case 9:
                    return " HD1 ";
                    break;
                case 10:
                    return " CE1 ";
                    break;
                case 11:
                    return " HE1 ";
                    break;
                case 12:
                    return " NE2 ";
                    break;
                case 13:
                    return " CD2 ";
                    break;
                case 14:
                    return " HD2 ";
                    break;
                case 15:
                    return " HE2 ";
                    break;
                case 16:
                    return " C   ";
                    break;
                case 17:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //ILE
        case 'I':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB  ";
                    break;
                case 6:
                    return " CG2 ";
                    break;
                case 7:
                    return "1HG2 ";
                    break;
                case 8:
                    return "2HG2 ";
                    break;
                case 9:
                    return "3HG2 ";
                    break;
                case 10:
                    return " CG1 ";
                    break;
                case 11:
                    return "1HG1 ";
                    break;
                case 12:
                    return "2HG1 ";
                    break;
                case 13:
                    return " CD  ";
                    break;
                case 14:
                    return " HD1 ";
                    break;
                case 15:
                    return " HD2 ";
                    break;
                case 16:
                    return " HD3 ";
                    break;
                case 17:
                    return " C   ";
                    break;
                case 18:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //LEU
        case 'L':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " HG  ";
                    break;
                case 9:
                    return " CD1 ";
                    break;
                case 10:
                    return "1HD1 ";
                    break;
                case 11:
                    return "2HD1 ";
                    break;
                case 12:
                    return "3HD1 ";
                    break;
                case 13:
                    return " CD2 ";
                    break;
                case 14:
                    return "1HD2 ";
                    break;
                case 15:
                    return "2HD2 ";
                    break;
                case 16:
                    return "3HD2 ";
                    break;
                case 17:
                    return " C   ";
                    break;
                case 18:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //LYS
        case 'K':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " HG1 ";
                    break;
                case 9:
                    return " HG2 ";
                    break;
                case 10:
                    return " CD  ";
                    break;
                case 11:
                    return " HD1 ";
                    break;
                case 12:
                    return " HD2 ";
                    break;
                case 13:
                    return " CE  ";
                    break;
                case 14:
                    return " HE1 ";
                    break;
                case 15:
                    return " HE2 ";
                    break;
                case 16:
                    return " NZ  ";
                    break;
                case 17:
                    return " HZ1 ";
                    break;
                case 18:
                    return " HZ2 ";
                    break;
                case 19:
                    return " HZ3 ";
                    break;
                case 20:
                    return " C   ";
                    break;
                case 21:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //MET
        case 'M':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " HG1 ";
                    break;
                case 9:
                    return " HG2 ";
                    break;
                case 10:
                    return " SD  ";
                    break;
                case 11:
                    return " CE  ";
                    break;
                case 12:
                    return " HE1 ";
                    break;
                case 13:
                    return " HE2 ";
                    break;
                case 14:
                    return " HE3 ";
                    break;
                case 15:
                    return " C   ";
                    break;
                case 16:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //PHE
        case 'F':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " CD1 ";
                    break;
                case 9:
                    return " HD1 ";
                    break;
                case 10:
                    return " CE1 ";
                    break;
                case 11:
                    return " HE1 ";
                    break;
                case 12:
                    return " CZ  ";
                    break;
                case 13:
                    return " HZ  ";
                    break;
                case 14:
                    return " CE2 ";
                    break;
                case 15:
                    return " HE2 ";
                    break;
                case 16:
                    return " CD2 ";
                    break;
                case 17:
                    return " HD2 ";
                    break;
                case 18:
                    return " C   ";
                    break;
                case 19:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //PRO
        case 'P':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " CD  ";
                    break;
                case 2:
                    return " HD1 ";
                    break;
                case 3:
                    return " HD2 ";
                    break;
                case 4:
                    return " CG  ";
                    break;
                case 5:
                    return " HG1 ";
                    break;
                case 6:
                    return " HG2 ";
                    break;
                case 7:
                    return " CB  ";
                    break;
                case 8:
                    return " HB1 ";
                    break;
                case 9:
                    return " HB2 ";
                    break;
                case 10:
                    return " CA  ";
                    break;
                case 11:
                    return " HA  ";
                    break;
                case 12:
                    return " C   ";
                    break;
                case 13:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //SER
        case 'S':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " OG  ";
                    break;
                case 8:
                    return " HG  ";
                    break;
                case 9:
                    return " C   ";
                    break;
                case 10:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //THR
        case 'T':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB  ";
                    break;
                case 6:
                    return " CG2 ";
                    break;
                case 7:
                    return "1HG2 ";
                    break;
                case 8:
                    return "2HG2 ";
                    break;
                case 9:
                    return "3HG2 ";
                    break;
                case 10:
                    return " OG1 ";
                    break;
                case 11:
                    return " HG1 ";
                    break;
                case 12:
                    return " C   ";
                    break;
                case 13:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //TRP
        case 'W':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " CD1 ";
                    break;
                case 9:
                    return " HD1 ";
                    break;
                case 10:
                    return " NE1 ";
                    break;
                case 11:
                    return " HE1 ";
                    break;
                case 12:
                    return " CE2 ";
                    break;
                case 13:
                    return " CZ2 ";
                    break;
                case 14:
                    return " HZ2 ";
                    break;
                case 15:
                    return " CH2 ";
                    break;
                case 16:
                    return " HH2 ";
                    break;
                case 17:
                    return " CZ3 ";
                    break;
                case 18:
                    return " HZ3 ";
                    break;
                case 19:
                    return " CE3 ";
                    break;
                case 20:
                    return " HE3 ";
                    break;
                case 21:
                    return " CD2 ";
                    break;
                case 22:
                    return " C   ";
                    break;
                case 23:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //TYR
        case 'Y':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB1 ";
                    break;
                case 6:
                    return " HB2 ";
                    break;
                case 7:
                    return " CG  ";
                    break;
                case 8:
                    return " CD1 ";
                    break;
                case 9:
                    return " HD1 ";
                    break;
                case 10:
                    return " CE1 ";
                    break;
                case 11:
                    return " HE1 ";
                    break;
                case 12:
                    return " CZ  ";
                    break;
                case 13:
                    return " OH  ";
                    break;
                case 14:
                    return " HH  ";
                    break;
                case 15:
                    return " CE2 ";
                    break;
                case 16:
                    return " HE2 ";
                    break;
                case 17:
                    return " CD2 ";
                    break;
                case 18:
                    return " HD2 ";
                    break;
                case 19:
                    return " C   ";
                    break;
                case 20:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //VAL
        case 'V':
            switch(nAtom){

                case 0:
                    return " N   ";
                    break;
                case 1:
                    return " H   ";
                    break;
                case 2:
                    return " CA  ";
                    break;
                case 3:
                    return " HA  ";
                    break;
                case 4:
                    return " CB  ";
                    break;
                case 5:
                    return " HB  ";
                    break;
                case 6:
                    return " CG1 ";
                    break;
                case 7:
                    return "1HG1 ";
                    break;
                case 8:
                    return "2HG1 ";
                    break;
                case 9:
                    return "3HG1 ";
                    break;
                case 10:
                    return " CG2 ";
                    break;
                case 11:
                    return "1HG2 ";
                    break;
                case 12:
                    return "2HG2 ";
                    break;
                case 13:
                    return "3HG2 ";
                    break;
                case 14:
                    return " C   ";
                    break;
                case 15:
                    return " O   ";
                    break;
                default:
                    return "\0";
            }
            break;

        //SOL
        case 'O':
            switch(nAtom){

                case 0:
                    return " OW  ";
                    break;
                case 1:
                    return " HW1 ";
                    break;
                case 2:
                    return " HW2 ";
                    break;
                default:
                    return "\0";
            }
            break;

        default:
            return "xxx";
    }

}

//Convert the atom to a index
char AATools::atomIndex(string res){

    switch(this->molecule){

        //ALA
        case 'A':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res==" CA  ")
                    return 2;
                else if(res==" HA  ")
                    return 3;
                else if(res==" CB  ")
                    return 4;
                else if(res==" HB1 ")
                    return 5;
                else if(res==" HB2 ")
                    return 6;
                else if(res==" HB3 ")
                    return 7;
                else if(res==" C   ")
                    return 8;
                else if(res==" O   ")
                    return 9;
                break;

        //ARG
        case 'R':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res==" CA  ")
                    return 2;
                else if(res==" HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB1 ")
                    return 5;
                else if(res== " HB2 ")
                    return 6;
                else if(res== " CG  ")
                    return 7;
                else if(res== " HG1 ")
                    return 8;
                else if(res== " HG2 ")
                    return 9;
                else if(res==" CD  ")
                    return 10;
                else if(res==" HD1 ")
                    return 11;
                else if(res==" HD2 ")
                    return 12;
                else if(res==" NE  ")
                    return 13;
                else if(res==" HE  ")
                    return 14;
                else if(res==" CZ  ")
                    return 15;
                else if(res==" NH1 ")
                    return 16;
                else if(res=="1HH1 ")
                    return 17;
                else if(res=="2HH1 ")
                    return 18;
                else if(res==" NH2 ")
                    return 19;
                else if(res=="1HH2 ")
                    return 20;
                else if(res=="2HH2 ")
                    return 21;
                else if(res==" C   ")
                    return 22;
                else if(res==" O   ")
                    return 23;
                break;

        //ASN
        case 'N':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " CG  ")
                    return 7;
               else if(res== " OD1 ")
                    return 8;
               else if(res== " ND2 ")
                    return 9;
                else if(res=="1HD2 ")
                    return 10;
                else if(res=="2HD2 ")
                    return 11;
                else if(res==" C   ")
                    return 12;
                else if(res==" O   ")
                    return 13;
                break;

        //ASP
        case 'D':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res== " CA  ")
                    return 2;
                else if(res== " HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB1 ")
                    return 5;
                else if(res== " HB2 ")
                    return 6;
                else if(res==" CG  ")
                    return 7;
                else if(res== " OD1 ")
                    return 8;
                else if(res== " OD2 ")
                    return 9;
                else if(res==" C   ")
                    return 10;
                else if(res==" O   ")
                    return 11;
                break;


        //CYS
        case 'C':
               if(res==" N   ")
                    return 0;
               else if(res== " H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " SG  ")
                    return 7;
               else if(res== " HG  ")
                    return 8;
               else if(res== " C   ")
                    return 9;
                else if(res==" O   ")
                    return 10;
                break;

        //CYM
        case 'X':
                if(res==" N   ")
                    return 0;
               else if(res== " H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " SG  ")
                    return 7;
               else if(res== " C   ")
                    return 8;
               else if(res== " O   ")
                    return 9;
                break;

        //GLU
        case 'E':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " CG  ")
                    return 7;
               else if(res== " HG1 ")
                    return 8;
               else if(res== " HG2 ")
                    return 9;
                else if(res==" CD  ")
                    return 10;
                else if(res==" OE1 ")
                    return 11;
                else if(res==" OE2 ")
                    return 12;
                else if(res==" C   ")
                    return 13;
                else if(res==" O   ")
                    return 14;
                break;

        //GLN
        case 'Q':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " CG  ")
                    return 7;
               else if(res== " HG1 ")
                    return 8;
               else if(res== " HG2 ")
                    return 9;
                else if(res==" CD  ")
                    return 10;
                else if(res==" OE1 ")
                    return 11;
                else if(res==" NE2 ")
                    return 12;
                else if(res=="1HE2 ")
                    return 13;
                else if(res=="2HE2 ")
                    return 14;
                else if(res==" C   ")
                    return 15;
                else if(res==" O   ")
                    return 16;
                break;

        //GLY
        case 'G':
               if(res==" N   ")
                    return 0;
               else if(res== " H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA1 ")
                    return 3;
               else if(res== " HA2 ")
                    return 4;
               else if(res== " C   ")
                    return 5;
               else if(res== " O   ")
                    return 6;
                break;

        //HIS
        case 'H':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res== " CA  ")
                    return 2;
                else if(res== " HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB1 ")
                    return 5;
                else if(res== " HB2 ")
                    return 6;
                else if(res== " CG  ")
                    return 7;
                else if(res== " ND1 ")
                    return 8;
                else if(res== " HD1 ")
                    return 9;
                else if(res==" CE1 ")
                    return 10;
                else if(res==" HE1 ")
                    return 11;
                else if(res==" NE2 ")
                    return 12;
                else if(res==" CD2 ")
                    return 13;
                else if(res==" HD2 ")
                    return 14;
                else if(res==" HE2 ")
                    return 15;
                else if(res==" C   ")
                    return 16;
                else if(res==" O   ")
                    return 17;
                break;

        //ILE
        case 'I':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res== " CA  ")
                    return 2;
                else if(res== " HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB  ")
                    return 5;
                else if(res== " CG2 ")
                    return 6;
                else if(res== "1HG2 ")
                    return 7;
                else if(res== "2HG2 ")
                    return 8;
                else if(res== "3HG2 ")
                    return 9;
                else if(res==" CG1 ")
                    return 10;
                else if(res=="1HG1 ")
                    return 11;
                else if(res=="2HG1 ")
                    return 12;
                else if(res==" CD  " || res==" CD1 ")
                    return 13;
                else if(res==" HD1 ")
                    return 14;
                else if(res==" HD2 ")
                    return 15;
                else if(res==" HD3 ")
                    return 16;
                else if(res==" C   ")
                    return 17;
                else if(res==" O   ")
                    return 18;
                break;

        //LEU
        case 'L':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res== " CA  ")
                    return 2;
                else if(res== " HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB1 ")
                    return 5;
                else if(res== " HB2 ")
                    return 6;
                else if(res== " CG  ")
                    return 7;
                else if(res== " HG  ")
                    return 8;
                else if(res== " CD1 ")
                    return 9;
                else if(res=="1HD1 ")
                    return 10;
                else if(res=="2HD1 ")
                    return 11;
                else if(res=="3HD1 ")
                    return 12;
                else if(res==" CD2 ")
                    return 13;
                else if(res=="1HD2 ")
                    return 14;
                else if(res=="2HD2 ")
                    return 15;
                else if(res=="3HD2 ")
                    return 16;
                else if(res==" C   ")
                    return 17;
                else if(res==" O   ")
                    return 18;
                break;

        //LYS
        case 'K':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res== " CA  ")
                    return 2;
                else if(res== " HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB1 ")
                    return 5;
                else if(res== " HB2 ")
                    return 6;
                else if(res== " CG  ")
                    return 7;
                else if(res== " HG1 ")
                    return 8;
                else if(res== " HG2 ")
                    return 9;
                else if(res==" CD  ")
                    return 10;
                else if(res==" HD1 ")
                    return 11;
                else if(res==" HD2 ")
                    return 12;
                else if(res==" CE  ")
                    return 13;
                else if(res==" HE1 ")
                    return 14;
                else if(res==" HE2 ")
                    return 15;
                else if(res==" NZ  ")
                    return 16;
                else if(res==" HZ1 ")
                    return 17;
                else if(res==" HZ2 ")
                    return 18;
                else if(res==" HZ3 ")
                    return 19;
                else if(res==" C   ")
                    return 20;
                else if(res==" O   ")
                    return 21;
                break;

        //MET
        case 'M':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
                else if(res== " CA  ")
                    return 2;
                else if(res== " HA  ")
                    return 3;
                else if(res== " CB  ")
                    return 4;
                else if(res== " HB1 ")
                    return 5;
                else if(res== " HB2 ")
                    return 6;
                else if(res== " CG  ")
                    return 7;
                else if(res== " HG1 ")
                    return 8;
                else if(res== " HG2 ")
                    return 9;
                else if(res==" SD  ")
                    return 10;
                else if(res==" CE  ")
                    return 11;
                else if(res==" HE1 ")
                    return 12;
                else if(res==" HE2 ")
                    return 13;
                else if(res==" HE3 ")
                    return 14;
                else if(res==" C   ")
                    return 15;
                else if(res==" O   ")
                    return 16;
                break;

        //PHE
        case 'F':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " CG  ")
                    return 7;
               else if(res== " CD1 ")
                    return 8;
               else if(res== " HD1 ")
                    return 9;
                else if(res==" CE1 ")
                    return 10;
                else if(res==" HE1 ")
                    return 11;
                else if(res==" CZ  ")
                    return 12;
                else if(res==" HZ  ")
                    return 13;
                else if(res==" CE2 ")
                    return 14;
                else if(res==" HE2 ")
                    return 15;
                else if(res==" CD2 ")
                    return 16;
                else if(res==" HD2 ")
                    return 17;
                else if(res==" C   ")
                    return 18;
                else if(res==" O   ")
                    return 19;
                break;

        //PRO
        case 'P':
                if(res==" N   ")
                    return 0;
                else if(res==" CD  ")
                    return 1;
               else if(res== " HD1 ")
                    return 2;
               else if(res== " HD2 ")
                    return 3;
               else if(res== " CG  ")
                    return 4;
               else if(res== " HG1 ")
                    return 5;
               else if(res== " HG2 ")
                    return 6;
               else if(res== " CB  ")
                    return 7;
               else if(res== " HB1 ")
                    return 8;
               else if(res== " HB2 ")
                    return 9;
                else if(res==" CA  ")
                    return 10;
                else if(res==" HA  ")
                    return 11;
                else if(res==" C   ")
                    return 12;
                else if(res==" O   ")
                    return 13;
                break;

        //SER
        case 'S':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " OG  ")
                    return 7;
               else if(res== " HG  ")
                    return 8;
               else if(res== " C   ")
                    return 9;
                else if(res==" O   ")
                    return 10;
                break;

        //THR
        case 'T':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB  ")
                    return 5;
               else if(res== " CG2 ")
                    return 6;
               else if(res== "1HG2 ")
                    return 7;
               else if(res== "2HG2 ")
                    return 8;
               else if(res== "3HG2 ")
                    return 9;
                else if(res==" OG1 ")
                    return 10;
                else if(res==" HG1 ")
                    return 11;
                else if(res==" C   ")
                    return 12;
                else if(res==" O   ")
                    return 13;
                break;

        //TRP
        case 'W':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " CG  ")
                    return 7;
               else if(res== " CD1 ")
                    return 8;
               else if(res== " HD1 ")
                    return 9;
                else if(res==" NE1 ")
                    return 10;
                else if(res==" HE1 ")
                    return 11;
                else if(res==" CE2 ")
                    return 12;
                else if(res==" CZ2 ")
                    return 13;
                else if(res==" HZ2 ")
                    return 14;
                else if(res==" CH2 ")
                    return 15;
                else if(res==" HH2 ")
                    return 16;
                else if(res==" CZ3 ")
                    return 17;
                else if(res==" HZ3 ")
                    return 18;
                else if(res==" CE3 ")
                    return 19;
                else if(res==" HE3 ")
                    return 20;
                else if(res==" CD2 ")
                    return 21;
                else if(res==" C   ")
                    return 22;
                else if(res==" O   ")
                    return 23;
                break;

        //TYR
        case 'Y':
               if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB1 ")
                    return 5;
               else if(res== " HB2 ")
                    return 6;
               else if(res== " CG  ")
                    return 7;
               else if(res== " CD1 ")
                    return 8;
               else if(res== " HD1 ")
                    return 9;
                else if(res==" CE1 ")
                    return 10;
                else if(res==" HE1 ")
                    return 11;
                else if(res==" CZ  ")
                    return 12;
                else if(res==" OH  ")
                    return 13;
                else if(res==" HH  ")
                    return 14;
                else if(res==" CE2 ")
                    return 15;
                else if(res==" HE2 ")
                    return 16;
                else if(res==" CD2 ")
                    return 17;
                else if(res==" HD2 ")
                    return 18;
                else if(res==" C   ")
                    return 19;
                else if(res==" O   ")
                    return 20;
                break;

        //VAL
        case 'V':
                if(res==" N   ")
                    return 0;
                else if(res==" H   ")
                    return 1;
               else if(res== " CA  ")
                    return 2;
               else if(res== " HA  ")
                    return 3;
               else if(res== " CB  ")
                    return 4;
               else if(res== " HB  ")
                    return 5;
               else if(res== " CG1 ")
                    return 6;
               else if(res== "1HG1 ")
                    return 7;
               else if(res== "2HG1 ")
                    return 8;
               else if(res== "3HG1 ")
                    return 9;
                else if(res==" CG2 ")
                    return 10;
                else if(res=="1HG2 ")
                    return 11;
                else if(res=="2HG2 ")
                    return 12;
                else if(res=="3HG2 ")
                    return 13;
                else if(res==" C   ")
                    return 14;
                else if(res==" O   ")
                    return 15;
                break;

        //SOL
        case 'O':
                if(res==" OW  ")
                    return 0;
                else if(res==" HW1 ")
                    return 1;
                else if(res== " HW2 ")
                    return 2;
                break;
    }
        return 100;
}

//Returns a specific atom from a residue
unsigned int AATools::findAtom(char code, Molecule &mol){

    AATools aaTools;
    aaTools.setMolecule(mol.getType());
    string codeConv;

    if(code=='O')
        codeConv = " O   ";
    else if(code=='N')
        codeConv = " N   ";
    else if(code=='H')
        codeConv = " H   ";
    else if(code=='C')
        codeConv = " CA  ";

    for(unsigned int i = 0; i < mol.atoms.length(); i++){

        if(aaTools.getAtom(mol.atoms[i].atom)==codeConv)
            return i;
    }

    return 0xFFFFFFFF;

}

//Return a specific atom from a interaction, true for A, false for B
unsigned int AATools::findAtom(char code, Bond &bond, bool type){

    AATools aaTools;
    if(type)
        aaTools.setMolecule(bond.getTypeD());
    else
        aaTools.setMolecule(bond.getTypeA());
    string codeConv;

    if(code=='O')
        codeConv = " O   ";
    else if(code=='N')
        codeConv = " N   ";
    else if(code=='H')
        codeConv = " H   ";
    else if(code=='C')
        codeConv = " CA  ";

    if(type){

        for(unsigned int i = 0; i < bond.atomsD.length(); i++){

            if(aaTools.getAtom(bond.atomsD[i].atom)==codeConv)
                return i;
        }
    }else{

        for(unsigned int i = 0; i < bond.atomsA.length(); i++){

            if(aaTools.getAtom(bond.atomsA[i].atom)==codeConv)
                return i;
        }
    }

    return 0xFFFFFFFF;

}

