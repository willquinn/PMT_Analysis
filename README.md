# PMT_Analysis

Set Up
------

This repo contains a Class file for reading any form of PMT waveform.

Requirements to use this code:
  - Must have python3
  - Must have pyroot
  
I include here, as it may be useful to others a quick tutorial in how to set up these two things. However, I will only do so for Mac OS as I have only figured it out for such.

First get Homebrew installed - https://brew.sh

Then install python3 and pip3:

  	$ brew install python3
  	$ brew install pip3
  
Now we have to install root. This is not simple as if you run brew install root6, it will work, but not for python3. Brew uses as default python2 and when you build root in this manner, importing it into python3 doesn't work. I you want to use python2 then you have no issue.

We have to build root from source. Just incase, you'll need:

  	$ brew install cmake
  	$ brew install make
  
Now go to root https://root.cern.ch/releases and get the latest release. Download this file to /usr/local and unzip. you willneed root priveleges.

Type the following:

  	$ sudo mkdir root6
  	$ cd root6
  
The next one is important you know which version of python you have downloaded (type python3 --version). Also change $Directory$ in the next line to the unzipped file you got from root-cern website:

  	$ sudo cmake ../$Directory$ -DPYTHON_EXECUTABLE=/usr/local/bin/python3 -DPYTHON_INCLUDE_DIR=/usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib -DPYTHON_LIBRARY=/usr/local/Cellar/python/3.7.6_1/Frameworks/Python.framework/Versions/3.7/lib/libpython3.7.dylib
  	$ sudo make
  
The last line should take a while and you will see a percentage number as the first character in each line.

Now in your .bash_profile add:

  	$ cd /usr/local/root6
  	$ source bin/thisroot.sh
  	$ cd $HOME

Once that is done you can now run python and import ROOT. And you shouldn't get any error.

Usage
-----
Now you should be able to run the test.py file. This will read a sample set of waveforms I have provided and plot some of these data to the screen and save a few files to a test root file. If this works then everything is running fine. The output file will be in the directory that you clone this repo and I advise you delete it.

The main file to focus on in this repo is the PMT_Classes.py file. The other files are examples of how I am using this code. The code is designed to be universal so that it can read any form of waveform.

With python3, pyroot installed properly and a good IDE you should be able to familiarise yourself with all the functionality of this trivial OOP for PMTs. I have tried to write a function for everything and have made the naming almost too clear so there should be no problem.
  
  
 
