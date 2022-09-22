This repo contains various, unrelated functions I wrote.  

`bandline.R` implements Stephen Few's BandLines;  

`myTransform.R` is a modification of `sp::spTransform` which does not throw an error when some points are unprojectable: it deletes them instead and cut the lines if needed;  

`pacman.py` implements Lazarus et al. 2012 pacman algorithm in python while `pacnsub.R` does so in R, plus implements various subsampling methods (SQS, UW, O2W, OW, and classic rarefaction) on the "pacman'd" dataset;  

`plot.digitizer.R` is a script to digitize data points from a publication figure;   

`sample.depth.R` computes the depth of a DSDP or ODP sample based on its name;

`nsb.jl` is a script in julia containing functions and examples to connect and uses the NSB database. Those functions are meant to produce the exact same results as their counterpart in package [NSBcompanion](http://github.com/plannapus/NSBcompanion).
