File=main
if g++ -O3 -Wall -std=c++11 -I ./include/ -o $File $File.cpp;
then
	time ./$File;
	python3 ./Other_Code/Python/plotter_2.py
else
	echo "^^^^^^^^^^^^^^^^^^^^^^^^^
COULD NOT COMPILE
--------------------------";
fi
