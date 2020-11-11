# SPARKY_to_NMSTAR3.1
Converts SPARKY peaklists to NMRSTAR3.1

Converts SPARKY nhsqc, hnca, hncacb, hncoca, hnco, hbhaconh, chsqc, cch tocsy, and hcch tocsy peaklists to NMRSTAR3.1

While there is a built in plug in in SPARKY to do this already, that plug in does not tell you if you have mislabeled peaks (i.e. accidently labeling a CB as an HB in a cch tocsy). This script provides an option to set an standard deviation value that will print out any values above that standard deviation, and print-out what spectrum that value came from. 

E.G.
VAL 233 HA has a standard deviation of 5.4, when the threshold was set at 0.25
[chsqc,hbhaconh,cch,hcch]
[4.35,4.34,4.55,60.29]

Each value corresponds to the spectrum it came from, in the above example we can see the mislabeled peak is from hcch, which results in the high standard deviation. 

Additionally, this script can also compare your assigned chemical shifts to BMRB average chemical shift values, and will print out amino acids whos atoms are outside the BMRB range.  

This will enable easy checking of mislabeled, or even misassigned peaks within your peaklists. And if your peaklists are good, it has a generated NMRSTAR 3.1 file that can be immediately inputted into other programs (TALOS, CYANA once converted to XEASY https://bmrb.pdbj.org/software/conv/, etc.)

***To Use***

Simply put all peaklist files and sequence file in the directory with this script. ***Make sure the BMRB.csv file is also in the same directory***Go into the script and modify the various parameters with the names of your peaklist and sequence file. Leave it blank if you don't have a peaklist for a spectra

I.E.
```
nhsqc_file='nhsqc.list'
hnca_file='hnca.list'
hncacb_file='hncacb.list'
hncoca_file='' 
```
The standard deviation value is the threshold. Peaks below this value will be displayed (as in the example above). 

***For BMRB values, I have added a 5% variation from the range of bmrb values. Sometimes your values will only be slightly different from the range (1.34 for range 1.344-1.54), this 5% will account for that so only show values that are significantly outside the range are displayed).***
