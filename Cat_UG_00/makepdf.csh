#----------------------------#
# Create PDF of User's Guide #
#----------------------------#

pdflatex cat_ug
bibtex   cat_ug
pdflatex cat_ug

