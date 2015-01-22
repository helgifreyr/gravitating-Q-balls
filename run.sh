# $1 = nr
# $2 = frequency
# $3 = alpha
# $4 = # of runs

echo "$1" > parameters
echo "$2" >> parameters
echo "$3" >> parameters
for i in `seq 1 $4`;
do
  ./v1 | tee output.txt
done    
math -script extract-global-data.m
python extract-global-data.py
mkdir -p m=$1/a=$3/w$5=$2
mv *txt m=$1/a=$3/w$5=$2/
cp *dat m=$1/a=$3/w$5=$2/
