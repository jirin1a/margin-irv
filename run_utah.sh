#for i in data/profile_format/*txt; do echo "$i"; ./analyzeirv -task 0 -ballots $i -logfile logs/task0/log.$(echo $i|sed "s/.*\///g").txt -debug; done
#for i in data/profile_format/*txt; do echo "$i"; ./analyzeirv -task 1 -ballots $i -logfile logs/task1/log.$(echo $i|sed "s/.*\///g").txt -debug; done
#for i in data/profile_format/*txt; do echo "$i"; ./analyzeirv -task 2 -ballots $i -logfile logs/task2/log.$(echo $i|sed "s/.*\///g").txt -debug; done
for i in data/profile_format/*txt; do echo "$i"; ./analyzeirv -task 3 -ballots $i -logfile logs/task3/log.$(echo $i|sed "s/.*\///g").txt -debug; done
