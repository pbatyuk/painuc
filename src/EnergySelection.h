// flv  - 25 March 2010
//
//  header file introduced to automatize the compilation of the 4 version of the
//  program:  positive / negatice and 160Mev / 68 MeV
//  The variable "PI70" is defined in the script ~/tosello/SOURCE/Build_all.sh
//  before launching the make
//

#ifndef PI70
#define P_MEAN 210.//214; 194 -- pri T kin =100
#define MFIELD 6.5 // magnetic fieald in order to have R the same as by T=100
#endif
#ifdef PI70
#define P_MEAN 156. // pri T kin=70
#define MFIELD 5.2 // magnetic fieald in order to have R the same as by T=100
#endif
