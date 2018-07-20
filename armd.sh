#!/bin/sh

if [ -e datafiles ];
then
if [ -d datafiles ];
then
echo 'Using the already extant datafiles directory'
else
echo '***ERROR*** The file datafiles exists and is not a directory'
exit 1
fi
else
echo 'Creating the datafiles directory'
mkdir datafiles
fi

echo 'Running gp to make the data files'
gp -q -s 300000000 armd.gp > /dev/null

cd datafiles
echo 'Cleaning the data files'
for x in A*.txt m*.txt M*.txt
do
grep -v '^\?' $x | sed 's/ E/e/' > .tempfile.123 &&\
 echo 'END' >> .tempfile.123 && mv .tempfile.123 $x
done
cd ..

echo 'Turning the text data files into binaries'
for x in datafiles/*M.txt
do
NUM=`grep -c AT $x`
./sympow -txt2bin $NUM ${x//txt/bin} < $x
# rm -f $x
done
