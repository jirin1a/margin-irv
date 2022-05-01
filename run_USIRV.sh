for i in USIRV/*txt; do echo "$i"; ./analyzeirv -ballots $i -logfile logs/log.$(echo $i|sed "s/.*\///g").txt -debug; done
