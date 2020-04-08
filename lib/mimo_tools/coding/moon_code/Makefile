#
#	Make file for putting code over for class use
#

CLASSNUM  = 7670
CLASS = ee$(CLASSNUM)

SRCDIR = /home/tmoon/classes/$(CLASS)
BASEDIR = /home/tmoon/classes
HTMLDIR = $(BASEDIR)/$(CLASS)/web
INSTALLDIR = /opt/www/classes/ece/$(CLASSNUM)/provided
PROVDIR = /home/tmoon/eccbook/provided

install:
	makedirhtml.p
	scp -r $(PROVDIR)/* tmoon@sol:$(INSTALLDIR)

