#!/bin/bash
pdflatex ASAS_product_chart.tex
convert -density 300 ASAS_product_chart.pdf -quality 100 ASAS_product_chart.png
rm ASAS_product_chart.aux ASAS_product_chart.log ASAS_product_chart.pdf
