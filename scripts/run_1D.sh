#!/bin/sh

if [ $# -lt 3 ] 
    then 
  
    echo """this script expects 2 arguments
            model label
            model index: 0=0+/0-, 1=0+/2m+, 
                         2=0+/0h+, 3=0+/1+, 
                         4=0+/1-,  5=0+/qq2m+
			 6=0+/2h+, 7=0+/2h-
			 8=0+/2b+
                        (default=-1, dummy, no effect)
            mu option: 0--> mu = 1 (default)
                       1--> prefitmu
         """
    exit
fi

cd cards_$1_1D/
cp ../../SignalSeparation/ex* .


combineCards.py hzz4l_2e2muS_7TeV_ALT.txt hzz4l_2e2muS_8TeV_ALT.txt hzz4l_4eS_7TeV_ALT.txt hzz4l_4eS_8TeV_ALT.txt hzz4l_4muS_7TeV_ALT.txt hzz4l_4muS_8TeV_ALT.txt > hzz4l_ALT.txt

#./execute_SignalSeparationCombine.sh ./ hzz4l_ALT.txt 2 $3
#root -l -n -b -q "extractSignificanceStats.C+ (true,\"$1\",\"$1\")"
