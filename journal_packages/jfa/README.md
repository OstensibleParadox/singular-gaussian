# Journal of Functional Analysis submission package

Journal: Journal of Functional Analysis
Publisher: Elsevier
Official guide: https://www.sciencedirect.com/journal/journal-of-functional-analysis/publish/guide-for-authors
Official LaTeX instructions: https://www.elsevier.com/researcher/author/policies-and-guidelines/latex-instructions
Template URL: https://assets.ctfassets.net/o78em1y1w4i4/4MpsJHO0MOJ2xZuwGTAbOZ/7bc64af36477c5d6cfce335a1f872363/elsarticle.zip

Package files:
- manuscript.pdf: current anonymous compiled manuscript PDF, submitted as a standalone file.
- jfa_cover_letter.pdf: cover letter, submitted as a standalone file.
- jfa_anonymous_source.zip: flat LaTeX source package.

Source package notes:
- main.tex uses Elsevier elsarticle preprint format.
- main.tex sets journal to Journal of Functional Analysis.
- author line remains Anonymous Author(s).
- bibliography is supplied through main.bbl and paper1_ref.bib.
- figures are flattened into source root; no subfolders inside source zip.
- elsarticle.cls and elsarticle-num.bst are included.
- manuscript and cover letter PDFs are not included in the source zip.

Submission reminder:
- Elsevier says source files should be editable; PDF alone is not source.
- Elsevier says source files include tex, bib, figures, classes/styles not in TeX Live.
- Elsevier EM cannot process LaTeX files in folder structures; source zip is flat.
