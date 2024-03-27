all:
	xelatex report

release:
	mv report.pdf ..

clean:
	$(shell if [ -e report.bcf ];then rm report.bcf; fi)
	$(shell if [ -e report.log ];then rm report.log; fi)
	$(shell if [ -e report.out ];then rm report.out; fi)
	$(shell if [ -e report.toc ];then rm report.toc; fi)
	$(shell if [ -e report.run.xml ];then rm report.run.xml; fi)
	rm report.aux