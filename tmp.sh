#!/bin/sh


rm -f out/test_suite/GMD2017/CliffFineSandBays/*
cp in/test_suite/GMD2017/CliffFineSandBays/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/GMD2017/UndefendedCoastline/*
cp in/test_suite/GMD2017/UndefendedCoastline/cme.ini .
./cme
echo ===============================================================================

rm -f out/CSE/Typology/Cliff/*
cp in/CSE/0_Typology/Cliff/cme.ini .
./cme
echo ===============================================================================

rm -f out/CSE/Typology/Dune/*
cp in/CSE/0_Typology/Dune/cme.ini .
./cme
echo ===============================================================================

rm -f out/two_coasts/*
cp in/two_coasts/cme.ini .
./cme
echo ===============================================================================

mv cme.ini.OLD cme.ini 2> /dev/null
mv cme.yaml.OLD cme.yaml 2> /dev/null

# OLD ===============================================================================
# rm -f out/test_suite/Wilf/Scen001/*
# cp in/test_suite/Wilf/Scen001/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen002/*
# cp in/test_suite/Wilf/Scen002/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen003/*
# cp in/test_suite/Wilf/Scen003/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen004/*
# cp in/test_suite/Wilf/Scen004/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen005/*
# cp in/test_suite/Wilf/Scen005/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen006/*
# cp in/test_suite/Wilf/Scen006/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen007/*
# cp in/test_suite/Wilf/Scen007/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen008/*
# cp in/test_suite/Wilf/Scen008/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen009/*
# cp in/test_suite/Wilf/Scen009/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen010/*
# cp in/test_suite/Wilf/Scen010/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen011/*
# cp in/test_suite/Wilf/Scen011/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen012/*
# cp in/test_suite/Wilf/Scen012/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen013/*
# cp in/test_suite/Wilf/Scen013/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen014/*
# cp in/test_suite/Wilf/Scen014/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen015/*
# cp in/test_suite/Wilf/Scen015/cme.ini .
# ./cme
# echo ===============================================================================
#
# rm -f out/test_suite/Wilf/Scen016/*
# cp in/test_suite/Wilf/Scen016/cme.ini .
# ./cme
# echo ===============================================================================


