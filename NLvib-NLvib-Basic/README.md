# NLvib
NLvib is a Matlab tool for nonlinear vibration analysis.

An overview of its capabilities, included examles, the monograph on Harmonic Balance, and some further presentation material can be found on https://www.ila.uni-stuttgart.de/nlvib/.

Please also see the manual in the SRC folder.

We always appreciate any kind of feedback you may have. If you encounter any problems, which you cannot solve with the manual or the book, please do not hesitate to contact the authors of this code (Johann Gross and Malte Krack; for contact details see headers or https://www.ila.uni-stuttgart.de/nlvib/).

The tool should work well with a wide range of Matlab releases. It mainly relies on the optimization toolbox.
To get it to run under OCTAVE, you need to change the line(s) that thes the solver options to something like:
   Solopt = optimset(optimset ("fsolve"),'Display','off',... 'Jacobian','on','MaxIter',50);
Also, some 'legend' calls (for figures) might require modification.

# Branches 
The main branch of NLvib is now split in two options:
1. `NLvib - Basic`: former main branch
2. `NLvib - PEACE`: adapted implementation for model refinement capabilities

All EXAMPLES included in `NLvib - Basic` can be run in `NLvib - PEACE`, but not vice versa.
`NLvib - Basic` is kept as a branch as it is more compact.
