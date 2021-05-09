## updating website
Use `rmarkdown` to render Rmd to md
```
rmarkdown::render("vignettes/getting_started.Rmd", rmarkdown::md_document(variant = "markdown_github"))
```
Then copy over image and md files to `docs/`