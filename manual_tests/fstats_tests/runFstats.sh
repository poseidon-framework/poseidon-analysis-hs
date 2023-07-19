#!/usr/bin/env bash

stack run xerxes -- fstats --statFile fstats.txt -d ../../../published_data/2012_PattersonGenetics -f testTable.txt --blockTableFile testBlockTable.txt

# .-----------.-----------.---------------.--------.-------.---------.------------.-----------.---------------------.
# | Statistic |     a     |       b       |   c    |   d   | NrSites |  Estimate  |  StdErr   |       Z score       |
# :===========:===========:===============:========:=======:=========:============:===========:=====================:
# | F3        | French    | Italian_North | Mbuti  |       | 593124  | 0.2576     | 2.6932e-3 | 95.64402147158148   |
# | F3        | French    | Han           | Mbuti  |       | 593124  | 0.2168     | 2.5383e-3 | 85.39316355835682   |
# | F3        | Sardinian | Pima          | French |       | 593124  | -5.3477e-3 | 3.9419e-4 | -13.566531670301005 |
# | F4        | French    | Russian       | Han    | Mbuti | 593124  | -1.6778e-3 | 9.1419e-5 | -18.35262346091248  |
# | F4        | Sardinian | French        | Pima   | Mbuti | 593124  | -1.4384e-3 | 1.1525e-4 | -12.481084899924868 |
# '-----------'-----------'---------------'--------'-------'---------'------------'-----------'---------------------'

# stack run xerxes -- fstats --statConfig config.yaml -d ../../../published_data --maxSnps 100000

