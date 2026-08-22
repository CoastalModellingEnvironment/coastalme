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
