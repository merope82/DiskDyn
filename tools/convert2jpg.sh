out=`echo $1 | sed s/.fits//g`
fitspng -o ${out}.png -fl 0,0.0000007 -f sqr $1
