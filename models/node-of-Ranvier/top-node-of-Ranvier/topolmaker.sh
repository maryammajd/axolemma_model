#to add include files in topology file
ls * > top.top
#sed:command. ^:prefix $:appendix \:for special chars(such as /). 
sed 's/^/#include "top\//; s/$/"/' top.top
