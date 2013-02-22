#! /bin/bash


DIR7TEV=$1
DIR8TEV=$2
MH=$3

cp ${DIR7TEV}/HCG/${MH}/hzz*7TeV*.txt ${DIR8TEV}/HCG/${MH}/
cp ${DIR7TEV}/HCG/${MH}/hzz*7TeV*.root ${DIR8TEV}/HCG/${MH}/

cd ${DIR8TEV}/HCG/${MH}/

combineCards.py hzz4l_2e2muS_7TeV.txt hzz4l_4eS_7TeV.txt hzz4l_4muS_7TeV.txt > hzz4l_4lS_7TeV.txt
combineCards.py hzz4l_2e2muS_7TeV_ALT.txt hzz4l_4eS_7TeV_ALT.txt hzz4l_4muS_7TeV_ALT.txt > hzz4l_4lS_7TeV_ALT.txt

combineCards.py hzz4l_2e2muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_4muS_8TeV.txt > hzz4l_4lS_8TeV.txt
combineCards.py hzz4l_2e2muS_8TeV_ALT.txt hzz4l_4eS_8TeV_ALT.txt hzz4l_4muS_8TeV_ALT.txt > hzz4l_4lS_8TeV_ALT.txt

#combineCards.py hzz4l_2e2muS_7TeV.txt hzz4l_4eS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_2e2muS_8TeV.txt hzz4l_4eS_8TeV.txt hzz4l_4muS_8TeV.txt > hzz4l_4lS_7and8TeV.txt
combineCards.py   hzz4l_2e2muS_7TeV.txt   hzz4l_2e2muS_8TeV.txt hzz4l_4muS_8TeV.txt hzz4l_4eS_7TeV.txt hzz4l_4muS_7TeV.txt hzz4l_4eS_8TeV.txt > hzz4l_4lS_7and8TeV.txt
combineCards.py hzz4l_2e2muS_7TeV_ALT.txt hzz4l_4eS_7TeV_ALT.txt hzz4l_4muS_7TeV_ALT.txt hzz4l_2e2muS_8TeV_ALT.txt hzz4l_4eS_8TeV_ALT.txt hzz4l_4muS_8TeV_ALT.txt > hzz4l_4lS_7and8TeV_ALT.txt
