#!/usr/bin/env bash

stack run xerxes -- fstats --statFile fstats.txt -d ../../../community-archive/2012_PattersonGenetics -f testTable.txt --blockTableFile testBlockTable.txt

# .-----------.-----------.---------------.--------.-------.---------.----------------.--------------------.------------------.
# | Statistic |     a     |       b       |   c    |   d   | NrSites | Estimate_Total | Estimate_Jackknife | StdErr_Jackknife |
# :===========:===========:===============:========:=======:=========:================:====================:==================:
# | F3        | French    | Italian_North | Mbuti  |       | 593124  | 5.9698e-2      | 5.9698e-2          | 5.1423e-4        |
# | F3        | French    | Han           | Mbuti  |       | 593124  | 5.0233e-2      | 5.0233e-2          | 5.0324e-4        |
# | F3        | Sardinian | Pima          | French |       | 593124  | -1.2483e-3     | -1.2483e-3         | 9.2510e-5        |
# | F4        | French    | Russian       | Han    | Mbuti | 593124  | -1.6778e-3     | -1.6778e-3         | 9.1419e-5        |
# | F4        | Sardinian | French        | Pima   | Mbuti | 593124  | -1.4384e-3     | -1.4384e-3         | 1.1525e-4        |
# '-----------'-----------'---------------'--------'-------'---------'----------------'--------------------'------------------'

stack run xerxes -- fstats --statConfig config.yaml -d ../../../community-archive --maxSnps 100000 -f testTable2.txt --blockTableFile testBlockTable2.txt

