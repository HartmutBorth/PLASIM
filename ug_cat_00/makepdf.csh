#----------------------------#
# Create PDF of User's Guide #
#----------------------------#

pdflatex ug_cat
bibtex   ug_cat
pdflatex ug_cat

