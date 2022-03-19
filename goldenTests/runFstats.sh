#!/usr/bin/env bash

stack run xerxes -- fstats --statFile fstats.txt -d ../../published_data/2012_PattersonGenetics

# .---------------------------------.------------------------.-----------------------.---------------------.
# |            Statistic            |        Estimate        |        StdErr         |       Z score       |
# :=================================:========================:=======================:=====================:
# | F3(French,Italian_North,Mbuti)  | 0.25759035709059575    | 2.6932196401541274e-3 | 95.64402147158499   |
# | F3(French,Han,Mbuti)            | 0.21675082698270387    | 2.5382690832689257e-3 | 85.39316355835685   |
# | F3(Sardinian,Pima,French)       | -5.347734780554392e-3  | 3.9418584723910377e-4 | -13.566531670302671 |
# | F4(French,Russian,Han,Mbuti)    | -1.677777202063968e-3  | 9.14189301402711e-5   | -18.35262346091368  |
# | F4(Sardinian,French,Pima,Mbuti) | -1.4383893357278188e-3 | 1.1524553732797402e-4 | -12.481084899924126 |
# '---------------------------------'------------------------'-----------------------'---------------------'
