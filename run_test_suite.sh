#!/bin/sh

cp cme.ini cme.ini.OLD 2> /dev/null
cp cme.yaml cme.yaml.OLD 2> /dev/null

rm -f out/test_suite/minimal_cons_wave_angle_230/*
cp in/test_suite/minimal_cons_wave_angle_230/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_uncons_wave_angle_230/*
cp in/test_suite/minimal_uncons_wave_angle_230/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_wave_angle_270/*
cp in/test_suite/minimal_cons_wave_angle_270/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_wave_angle_310/*
cp in/test_suite/minimal_cons_wave_angle_310/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_215/*
cp in/test_suite/minimal_cons_with_intervention_wave_angle_215/cme.ini .
./cme
echo ===============================================================================
#
# rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_235/*
# cp in/test_suite/minimal_cons_with_intervention_wave_angle_235/cme.ini .
# ./cme
# echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_270/*
cp in/test_suite/minimal_cons_with_intervention_wave_angle_270/cme.ini .
./cme
echo ===============================================================================

# rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_285/*
# cp in/test_suite/minimal_cons_with_intervention_wave_angle_285/cme.ini .
# ./cme
# echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_305/*
cp in/test_suite/minimal_cons_with_intervention_wave_angle_305/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_two_wavestations/*
cp in/test_suite/minimal_cons_with_two_wavestations/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_sediment_input_point/*
cp in/test_suite/minimal_cons_with_sediment_input_point/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_sediment_input_line/*
cp in/test_suite/minimal_cons_with_sediment_input_line/cme.ini .
./cme
echo ===============================================================================

# rm -f out/test_suite/minimal_cons_with_sediment_input_block/*
# cp in/test_suite/minimal_cons_with_sediment_input_block/cme.ini .
# ./cme
# echo ===============================================================================

rm -f out/test_suite/minimal_cons_sand_and_coarse_wave_angle_230/*
cp in/test_suite/minimal_cons_sand_and_coarse_wave_angle_230/cme.ini .
./cme
echo ===============================================================================
rm -f out/test_suite/Happisburgh/*
cp in/test_suite/Happisburgh/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/Manuel_C003_0001/*
cp in/test_suite/Manuel_C003_0001/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/GMD2017/CliffFineBays/*
cp in/test_suite/GMD2017/CliffFineBays/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/GMD2017/Groyne/*
cp in/test_suite/GMD2017/Groyne/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/GMD2017/CliffFineSandBays/*
cp in/test_suite/GMD2017/CliffFineSandBays/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/GMD2017/UndefendedCoastline/*
cp in/test_suite/GMD2017/UndefendedCoastline/cme.ini .
./cme
echo ===============================================================================

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
echo ===============================================================================

rm -f out/CSE/Typology/Cliff/*
cp in/CSE/0_Typology/Cliff/cme.ini .
./cme
echo ===============================================================================

mv cme.ini.OLD cme.ini 2> /dev/null
mv cme.yaml.OLD cme.yaml 2> /dev/null
