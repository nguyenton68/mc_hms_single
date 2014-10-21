mc_hms_single
============================================

mc_hms_single is the standard Hall C HMS single arm Monte Carlo.

Compiling code
------------
* Be careful with Makefile, you have to change the CERNLIB correspond to CERNLIB in your computer
To compile type "make"

Main code is mc_hms_single.f

If include files are changed it is best to do "make Clean" and then "make"
since the Makefile is not smart enough to look for dependency
[] make clean
[] make
[] ./mc_hms_single
[] test1

We can use .sh for convinient
Create runit.sh:
./mc_hms_single<<EOF1
test1
EOF1


./mc_hms_single<<EOF2
test2
EOF2


./mc_hms_single<<EOF3
test3
EOF3

./mc_hms_single<<EOF4
test4
EOF4

./mc_hms_single<<EOF5
test5
EOF5


Then:
[]make clean
[]make
[]chmod +x runit.sh
[]./runit.sh


Running code
------------

mc_hms_single 
(ask for input file name ( assumed in infiles subdirectory with .inp extension)

* Input file : infile_name
* Output file is at outfiles/infile_name.out 
* The hbook file is at worksim/infile_name.rzdat 

* mc_hms_single.f (main code): uncommand radlength_tuna(z,...) because this one create trouble and
we don't use this kind of target.
When merge subrountine externals() to the main mc_hms_single will create some trouble. 
Because there are so many variables name.

nmc_org is a separate subrountine
But F1F2IN09 is a subrountine inside external() subrountine.
* EXTERNALS(rcebeam,rcep,rctheta,x_total)
rcebeam, rcep, rctheta have to be real*8 because in the main mc_hms_single they are real*8.
But x_total is real
Inside EXTERNALS sub we have to input input files by hand: /TARGT/; /IKK12/
And put them in "data"
Then pass value of E_beam, E', theta to E,E',theta inside EXTERNALS.
Set TARGET="HE3"

