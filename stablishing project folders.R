library(renv)

# Optional but recommended: initialize Git
usethis::use_git()

# Reproducible env
renv::init()

# Project structure
usethis::use_directory("R")
usethis::use_directory("scripts")
usethis::use_directory("data/raw")
usethis::use_directory("data/processed")
usethis::use_directory("models")
usethis::use_directory("results")
usethis::use_directory("figures")
usethis::use_directory("docs")

# Basic README
usethis::use_readme_md()

# Ignore typical files/dirs
usethis::use_git_ignore(c(".Rhistory", ".RData", ".Rproj.user", "data/raw/", "models/", "results/"))
