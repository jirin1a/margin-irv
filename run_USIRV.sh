#for i in USIRV/*txt; do echo "$i"; ./analyzeirv -ballots $i -logfile logs/log.$(echo $i|sed "s/.*\///g").txt -debug; done
for i in USIRV/*txt; do echo "$i"; ./analyzeirv -task 1 -ballots $i -logfile logs/task1/log.$(echo $i|sed "s/.*\///g").txt -debug; done
