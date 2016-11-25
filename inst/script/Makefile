DESeq2_pdf.Rmd: DESeq2.Rmd
	sed -e 's/rmarkdown::html_document:/BiocStyle::pdf_document2:/' -e '/highlight: pygments/d' -e 's/Analyzing RNA-seq data with DESeq2/Analyzing RNA-seq data with DESeq2 (PDF)/g' -e 's/This is the source document/This is a derived document, DO NOT EDIT/' -e '/DESeq2 package version: `r packageVersion("DESeq2")`/d' DESeq2.Rmd > DESeq2_pdf.Rmd
