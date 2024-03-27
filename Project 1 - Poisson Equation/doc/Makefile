all:
	xelatex document

release:
	mv document.pdf ..

clean:
	$(shell if [ -e document.bcf ];then rm document.bcf; fi)
	$(shell if [ -e document.log ];then rm document.log; fi)
	$(shell if [ -e document.out ];then rm document.out; fi)
	$(shell if [ -e document.toc ];then rm document.toc; fi)
	$(shell if [ -e document.run.xml ];then rm document.run.xml; fi)
	rm document.aux