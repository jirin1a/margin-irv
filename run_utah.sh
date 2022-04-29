for i in data/profile_format/*txt; do echo "$i"; ./analyzeirv -ballots $i -logfile logs/log.$(echo $i|sed "s/.*\///g").txt -debug; done
