#!/bin/bash
# $1: open enhancers cnt
# $2: closed enhancers cnt
# $3: open enh random filelist "a@b@c"
# $4: closed enh random filelist 
# $5: rseed values
# $6: output file

echo -e "class\tfootprint_size" > $6
rseed_vals=`echo $5 | tr '@' ' '`
shuf_open_files=`echo $3 | tr '@' ' '` 
shuf_closed_files=`echo $4 | tr '@' ' '` 
echo $shuf_closed_files
echo $rseed_vals 
cnt=1 
for f in $shuf_open_files 
do 
    rseed=`echo $rseed_vals | awk '{print $v}' v=${cnt}`
    echo "${rseed}"
    awk '{print  "ropen_"v"\t"$NF}' v=${rseed} $f >> $6 
    cnt=`expr $cnt + 1`
done 

cnt=1 
for f in $shuf_closed_files
do 
    rseed=`echo $rseed_vals | awk '{print $v}' v=${cnt}`
    echo ${rseed}
    awk '{print  "rclose_"v"\t"$NF}' v=${rseed} $f >> $6 
    cnt=`expr $cnt + 1`
done 

awk '{print  "obs\t"$NF}' $1 >> $6
awk '{print  "close\t"$NF}' $2 >> $6
