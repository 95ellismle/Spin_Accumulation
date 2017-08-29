# Spin Accumulation

This where Matt's version of the spin accumulation code is kept.

### Explanation of files:
##### Include folder
The include folder is where all the peripheral code is kept i.e. the header files.

###### -- save.h:
&emsp; &emsp; This contains the functions used to save data.

###### --math_funcs.h
&emsp; &emsp; This contains all the maths related functions such as a numerical differentiator.

###### --parse_and_string_stuff.h
&emsp; &emsp; This contains any string related functions, and a function to parse integer command line arguments.

###### --solve.h
&emsp; &emsp; This contains the solving functions. The functions that implement the spin current equation and the explicit Euler scheme.

###### --test.h
&emsp; &emsp; This contains functions used to test the code.

###### --vec_fill.h
&emsp; &emsp; This only contains 1 function and it calculates diffusive parameters.

##### Other_Code
This just contains a couple python files I have been using to plot the data and the python simulation code. There's a separate summary of the files in that folder.

##### Useful_documents
Contains some useful documents. These are:
* Maths behind the Explcit Euler Scheme.
* A presentation with some initial results of the code.  

##### main.cpp
This is the code that links all the separate header files together. <b>Remember:</b><i>The filepath will need to be set here to make the code save in the correct place.</i>

##### compile_run.sh
This file can be ran and it will compile and run the main.cpp code. To run this you may have to use the command: <br> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; "chmod 755 compile_run.sh"
<br>
to make the file executable. To run that you can then use the command:
<br> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; "./compile_run.sh"
<br>
This should compile and run the code. It will use full optimisations and create a runnable file called `main'.
