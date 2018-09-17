## create small test dataset containing Oxalis enneaphylla and O. lacinita to 
## develop the code for the niche equivalency test and the background similarity
## test

g.remove vect=samples
cd /Users/Stoffi/R_files/oxalis/diversity 
v.in.ascii in=data/samples.txt out=samples x=1 y=2 columns='x double precision, y double precision, spec varchar(30), sect varchar(20), orig varchar(5)' --o

db.connect 