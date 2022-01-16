A repository containing the source code of the PAINUC experiment, has been created. The  codes in this  
repository are the old ones that complied with the old versions of the G++ compiler and the ROOT software. 
Our main goal is to update the codes so as to make them compatible with the modern G++ compiler and ROOT 
software.
   The codes are located in the src directory of this repository. The complete installation procedure can be 
started making use of compile.sh. 
   To verify correct execution, i.e. to be sure there are no segmentation violations and so on, one is referred to 
the negative/106_MeV directory (in this repository), where an example of an event to be measured is found.  

The automaton-dialog.cxx (put in the root of the repository) is the last version of program to be used for automatic measurements in interactive mode.
An emulation of batch mode can be achieved by commenting some obvious functions in the source code. The program is a 
standalone c++ program, so, presence of Makefile to be created is obligatory.
