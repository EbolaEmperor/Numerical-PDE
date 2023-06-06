all:
	xelatex slide
	xelatex slide

release:
	mv slide.pdf ..

clean:
	$(shell if [ -e slide.bcf ];then rm slide.bcf; fi)
	$(shell if [ -e slide.log ];then rm slide.log; fi)
	$(shell if [ -e slide.out ];then rm slide.out; fi)
	$(shell if [ -e slide.toc ];then rm slide.toc; fi)
	$(shell if [ -e slide.run.xml ];then rm slide.run.xml; fi)
	rm slide.aux
