## updating website
Use `rmarkdown` to render Rmd to md
```
render("vignettes/getting_started.Rmd", md_document(variant = "markdown_github"))
```
Then copy over image and md files to `docs/`