NOW=`date  +%Y-%m-%d`
PWD=`pwd`



default: build
	-cd devel; $(MAKE) 
	-cd examples/part1; $(MAKE)
	-cd examples/part2; $(MAKE)
	-cd examples/part3; $(MAKE)
	-cd examples/part4; $(MAKE)

build:
	-cd src; $(MAKE) 

release: clean doxygen changelog
	-rm -r ~/scratch/numcxx/numcxx-$(NOW)
	cd ..; rsync -avu --exclude=QUARRY --exclude=.hg  numcxx/* ~/scratch/numcxx/numcxx-$(NOW)
	cd ~/scratch/numcxx/numcxx-$(NOW); $(MAKE) test; $(MAKE) clean
	cd ~/scratch/numcxx; tar czvf numcxx-$(NOW).tgz numcxx-$(NOW)
	cp -p ~/scratch/numcxx/numcxx-$(NOW).tgz $(HOME)/Wias/www/fuhrmann/sitesrc/blobs
	cp -p ~/scratch/numcxx/numcxx-$(NOW).tgz $(HOME)/Wias/www/fuhrmann/sitesrc/blobs/numcxx-latest.tgz
	hg tag -f RELEASE_$(NOW)

doxygen: 
	doxygen

clean:  
	-rm -r lib
	-rm -r html
	-rm INSTALL.pdf
	-cd src; $(MAKE) clean
	-cd devel; $(MAKE) clean
	-cd examples/part1; $(MAKE) clean
	-cd examples/part2; $(MAKE) clean
	-cd examples/part3; $(MAKE) clean
	-cd examples/part4; $(MAKE) clean
 

test:
	cd devel; $(MAKE) ; $(MAKE) test; echo devel ok
	cd examples/part1; WIR=$(PWD) $(MAKE) -e; $(MAKE) test; echo part1 ok
	cd examples/part2; WIR=$(PWD) $(MAKE) -e; $(MAKE) test; echo part2 ok
	cd examples/part3; WIR=$(PWD) $(MAKE) -e; $(MAKE) test; echo part3 ok
	cd examples/part4; WIR=$(PWD) $(MAKE) -e; $(MAKE) test; echo part4 ok

changelog:
	hg log --template "{date(date, '%Y-%m-%d-%H%M')}:\n {fill(desc,60,'                ','                 ')}\n\n" > CHANGELOG.txt
