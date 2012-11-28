# This makefile is used to update templates for internationalization ("i18n")
# of ag.

# INSTRUCTIONS
# ------------
#
# To update existing templates with strings from the latest python files, 
# just run 'make'.
#
# To create a translation template for a new language (e.g. de_DE.UTF8), run:
#     LANG=de_DE.UTF8 make init
# 

# List of files calling gettext.
TRANSLATIONFILES=../ag.py ../calculator.py ../colors.py ../constraints.py ../crystal.py ../debug.py ../dft.py ../energyforces.py ../execute.py ../graphene.py ../graphs.py ../gtkexcepthook.py ../gui.py ../minimize.py ../movie.py ../nanoparticle.py ../nanotube.py ../progress.py ../pybutton.py ../quickinfo.py ../render.py ../repeat.py ../rotate.py ../scaling.py ../settings.py ../setupwindow.py ../simulation.py ../status.py ../surfaceslab.py ../widgets.py

i18n: ag.pot update-po compile

# This will update the English template (ag.pot) with English messages from 
# the Python source code.
ag.pot: ${TRANSLATIONFILES}
	xgettext --add-comments --language=Python --keyword=_ --output=ag.pot --msgid-bugs-address=ase-developers@listserv.fysik.dtu.dk --copyright-holder=CAMD --package-name=ase ${TRANSLATIONFILES}

# This will create an empty translation file ($LANG.po, where $LANG is
# an environment variable) from the English template ag.pot.
# The po-file header will have to be edited afterwards.
init: ag.pot
	mkdir -p ${LANG}/LC_MESSAGES
	msginit --locale=${LANG} --input=ag.pot --output-file=${LANG}/LC_MESSAGES/ag.po

# This will merge new/changed English strings from the template ag.pot
# into all existing translations (*.po), maximally
# reusing the existing translations.
update-po: ag.pot
	for FILE in $(wildcard */LC_MESSAGES/ag.po); do \
	  msgmerge --previous --update $$FILE ag.pot ;\
	done

# This will compile all translations (*.po) into binary files in gettext's
# search directory.
compile: */LC_MESSAGES/ag.po
	for DIR in $(wildcard */LC_MESSAGES/); do \
	  msgfmt -cv --output-file=$$DIR/ag.mo $$DIR/ag.po ;\
	done
