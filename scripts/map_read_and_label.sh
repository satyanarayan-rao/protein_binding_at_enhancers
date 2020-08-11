while read id label
do 
    footprint=`grep -w "$id" $2`
    echo -e "$footprint\t$label"
done < $1  
