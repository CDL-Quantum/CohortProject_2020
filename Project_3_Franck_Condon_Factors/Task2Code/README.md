# FC is a program to calculate Franck-Condon factors
#Copyright (C) 1999,2006, 2020 Pierre-Nicholas Roy

#############							##############
Note: This is not the full version, but a modified version by Matthew Schmidt
with permission from Pierre-Nicholas Roy for use with the CDL Cohort Project 2020.
The full version can be found on: https://github.com/pnroy/FC 
Modifications:
- added comments (including references to relevant papers)
- removed output files containing matrices and delta vector
- removed capability to input Force Constants instead of normal modes
#############	     	      	    	      	      	 	##############


#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 2
#of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


To compile:

1. Depending on what operating system you have, copy the corresponding makefile and
to "makefile".

Example: You are running linux
cp makefile-linux makefile

2. Compile the code:
make 

3. Run the code (example with input file V3)
./FCF_calc V3