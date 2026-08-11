#!/bin/sh

mv cme.ini cme.ini.OLD 2> /dev/null
mv cme.yaml cme.yaml.OLD 2> /dev/null

rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_215/*
cp in/test_suite/minimal_cons_with_intervention_wave_angle_215/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_270/*
cp in/test_suite/minimal_cons_with_intervention_wave_angle_270/cme.ini .
./cme
echo ===============================================================================

rm -f out/test_suite/minimal_cons_with_intervention_wave_angle_305/*
cp in/test_suite/minimal_cons_with_intervention_wave_angle_305/cme.ini .
./cme
echo ===============================================================================

