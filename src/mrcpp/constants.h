
/*
 *
 */

#pragma once

const double MachinePrec = 1.0e-15L;
const double MachineZero = 1.0e-14L;
const int MaxOrder = 41; ///< Maximum scaling order
const int MaxDepth = 30; ///< Maximum depth of trees
const int MaxScale = 31; ///< Maximum scale of trees
const int MinScale = -31; ///< Minimum scale of trees
const int MaxSepRank = 1000;
//Max number of orbitals stored temporarily. Larger->more memory
//Also max size of orbitalvector that can be sent with send_OrbVec
const int workOrbVecSize = 10;


namespace Axis {
const int None = -1;
const int X = 0;
const int Y = 1;
const int Z = 2;
}


enum Spin { Paired, Alpha, Beta };
enum FuncType {	Legendre, Interpol };
enum SplitType { ExactSplit, NormalSplit, FastSplit };
enum CV_Transform { Forward, Backward };
enum MW_Transform { Compression, Reconstruction };
enum XC_Type { XC_undefined, XC_lda, XC_gga };
enum Traverse { TopDown, BottomUp };

//Math constants
const double pi = 3.1415926535897932384626433832795;
const double root_pi = 1.7724538509055160273;
const double C_x = -0.73855876638202240588; //Dirac exchange constant


