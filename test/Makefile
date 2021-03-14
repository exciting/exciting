EXESER = exciting_serial
EXEMPI = exciting_purempi
EXEMPISMP = exciting_mpismp

all :: 
	make serial

serial :: 
	python3 runtest.py -a run -e $(EXESER)
	python3 runtest.py -a clean

mpi :: 
	python3 runtest.py -a run -e $(EXEMPI)
	python3 runtest.py -a clean

mpiandsmp ::
	python3 runtest.py -a run -e $(EXEMPISMP)
	python3 runtest.py -a clean

reference :: 
	python3 runtest.py -a ref -e $(EXESER)
	python3 runtest.py -a clean

cleantests :: 
	python3 runtest.py -a clean

